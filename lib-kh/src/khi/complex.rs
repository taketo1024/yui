//! The involutive Khovanov chain complex `CKhI(D, τ) = Cone(CKh(D) -^{Q(1+τ)}-> Q·CKh(D))`
//! for a strongly invertible link, where `Q² = 0` (Definition 2.2 of the reference).
//! Built here as the [`ChainMap::cone`] of `1 + τ`.
//!
//! Reference:
//! - T. Sano, "Involutive Khovanov homology and equivariant knots",
//!   Algebr. Geom. Topol. 25 (2025), 5059–5111.
//!   <https://doi.org/10.2140/agt.2025.25.5059>, <https://arxiv.org/abs/2404.08568>

use std::ops::{Index, RangeInclusive};
use std::sync::OnceLock;
use delegate::delegate;

use itertools::Itertools;
use yui_core::lc::Lc;
use yui_core::{EucRing, EucRingOps, IteratorExt, Ring, RingOps};
use yui_homology::{ChainComplex1, ChainMap, ToSeqString, ToTableString, GrMod1, GrMod2, Summand};
use yui_link::InvLink;

use crate::kh::{KhComplex, KhGen};
use crate::tng::builder::SymBuildConfig;
use crate::khi::KhIHomology;
use crate::khi::{KhIGen, KhIGenExt, from_cone_gen, to_cone_gen};
use crate::util::Bigraded;

pub type KhIChain<R> = Lc<KhIGen, R>;

pub type KhIComplexSummand<R> = Summand<KhIGen, R>;

#[derive(Clone)]
pub struct KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    inner: ChainComplex1<KhIGen, R>,
    canon_cycles: Vec<KhIChain<R>>,
    deg_shift: (isize, isize),
    cache_bigr: OnceLock<GrMod2<KhIGen, R>>,
}

impl<R> KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn new(l: &InvLink, h: &R, t: &R, reduced: bool) -> Self {
        Self::new_partial(l, h, t, reduced, None)
    }

    // restricts the build to `h_range`; since `KhI_i = C_i ⊕ C_{i-1}`, the cone needs `C` over `[a-1, b]`.
    pub fn new_partial(l: &InvLink, h: &R, t: &R, reduced: bool, h_range: Option<RangeInclusive<isize>>) -> Self {
        Self::new_with_config(l, h, t, reduced, SymBuildConfig { h_range, ..Default::default() })
    }

    // `config.h_range` is the desired cone range `[a, b]`; the cone needs `C` over `[a-1, b]`,
    // so the build range is shifted down by one while the rest of `config` is kept.
    pub fn new_with_config(l: &InvLink, h: &R, t: &R, reduced: bool, config: SymBuildConfig) -> Self {
        use crate::tng::builder::SymTngBuilder;

        if config.cone_cob {
            return Self::cone_cob_complex(l, h, t, reduced, config);
        }

        let config = SymBuildConfig {
            h_range: config.h_range.map(|r| KhComplex::<R>::clamp_h_range(l.inner(), reduced, r)),
            ..config
        };
        let h_range = config.h_range.clone();
        let build_config = SymBuildConfig {
            h_range: h_range.as_ref().map(|r| (*r.start() - 1) ..= *r.end()),
            ..config
        };
        let b = SymTngBuilder::from_inv_link(l, h, t, reduced).with_config(build_config).run();
        let tau_map = b.tau_map();

        let b = b.into_inner();
        let canon_cycles = b.eval_elements();
        let complex = b.into_tng_complex().into_raw_complex();
        let c = KhComplex::from_raw_complex(l.inner(), h, t, reduced, complex, canon_cycles);

        match h_range {
            Some(range) => Self::cone_of(c, tau_map, range),
            None => Self::from_kh_complex(c, tau_map),
        }
    }

    // Cobordism-level cone: ConeBuilder yields the coned TngComplex + canon classes directly; extract
    // KhIGen by stripping the cone bit (both for the complex and the canon cycles).
    fn cone_cob_complex(l: &InvLink, h: &R, t: &R, reduced: bool, config: SymBuildConfig) -> Self {
        use crate::tng::builder::ConeBuilder;

        let config = SymBuildConfig {
            h_range: config.h_range.map(|r| KhComplex::<R>::clamp_h_range(l.inner(), reduced, r)),
            ..config
        };
        let h_range = config.h_range.clone();
        let build_config = SymBuildConfig {
            h_range: h_range.as_ref().map(|r| (*r.start() - 1) ..= *r.end()),
            ..config
        };

        let cone = ConeBuilder::from_inv_link(l, h, t, reduced).with_config(build_config).run();
        let canon_cycles = cone.eval_elements().into_iter()
            .map(|z| z.map_keys(|x| from_cone_gen(&x)))
            .collect_vec();

        let raw = cone.into_tng_complex().into_raw_complex();
        let inner = raw.map_keys(from_cone_gen, to_cone_gen);
        let inner = match h_range {
            Some(range) => inner.truncated(range),
            None => inner,
        };

        let deg_shift = KhComplex::<R>::deg_shift_for(l.inner(), reduced);
        Self::new_impl(inner, canon_cycles, deg_shift)
    }

    pub fn new_no_simplify(l: &InvLink, h: &R, t: &R, reduced: bool) -> Self {
        assert_eq!(R::one() + R::one(), R::zero(), "char(R) != 2");
        assert!(!reduced || (l.base_pt().is_some() && t.is_zero()));

        let c = KhComplex::new_no_simplify(l.inner(), h, t, reduced);
        Self::from_kh_complex(c, crate::khi::tau::tau_map(l))
    }

    pub(crate) fn from_kh_complex<F>(c: KhComplex<R>, map: F) -> Self
    where F: Fn(&KhGen) -> KhGen + Send + Sync + 'static {
        // the cone extends one degree above the complex.
        let h_range = c.h_range();
        let h_range = *h_range.start() ..= (h_range.end() + 1);
        Self::cone_of(c, map, h_range)
    }

    pub(crate) fn cone_of<F>(c: KhComplex<R>, map: F, h_range: RangeInclusive<isize>) -> Self
    where F: Fn(&KhGen) -> KhGen + Send + Sync + 'static {
        let deg_shift = c.deg_shift();

        let canon_cycles = c.canon_cycles().iter().flat_map(|z| {
            let bz = z.clone().map_keys(KhIGen::from_left);
            let qz = z.clone().map_keys(KhIGen::from_right);
            [bz, qz]
        }).sorted_by_key(|z| z.keys().map(|x| x.rel_h_deg()).min().unwrap_or(0)).collect_vec();

        // KhI is the mapping cone of (1 + τ) : KhComplex → KhComplex.
        let one_plus_tau = ChainMap::new(c.inner(), c.inner(), 0, move |_, z| {
            z.clone() + z.apply(|x| Lc::from(map(x)))
        });
        let inner = one_plus_tau.cone(h_range, false);

        KhIComplex::new_impl(inner, canon_cycles, deg_shift)
    }

    pub(crate) fn new_impl(inner: ChainComplex1<KhIGen, R>, canon_cycles: Vec<KhIChain<R>>, deg_shift: (isize, isize)) -> Self {
        Self { inner, canon_cycles, deg_shift, cache_bigr: OnceLock::new() }
    }

    pub fn inner(&self) -> &ChainComplex1<KhIGen, R> {
        &self.inner
    }

    delegate! {
        to self.inner {
            pub fn support(&self) -> impl Iterator<Item = &isize> + '_;
            pub fn is_supported(&self, i: isize) -> bool;
            pub fn d_deg(&self) -> isize;
            pub fn d(&self, i: isize, z: &KhIChain<R>) -> KhIChain<R>;
            pub fn describe_d(&self) -> String;
            pub fn describe_d_at(&self, i: isize) -> String;
        }
    }

    pub fn deg_shift(&self) -> (isize, isize) {
        self.deg_shift
    }

    pub fn h_deg_of(&self, x: &KhIGen) -> isize {
        self.deg_shift.0 + x.rel_h_deg()
    }

    pub fn q_deg_of(&self, x: &KhIGen) -> isize {
        self.deg_shift.1 + x.rel_q_deg()
    }

    pub fn h_deg_of_chain(&self, z: &KhIChain<R>) -> isize {
        z.keys().map(|x| self.h_deg_of(x)).min().unwrap_or(0)
    }

    pub fn q_deg_of_chain(&self, z: &KhIChain<R>) -> isize {
        z.keys().map(|x| self.q_deg_of(x)).min().unwrap_or(0)
    }

    pub fn h_range(&self) -> RangeInclusive<isize> {
        self.support().copied().range().unwrap_or(0..=-1)
    }

    pub fn q_range(&self) -> RangeInclusive<isize> {
        self.support().flat_map(|&i|
            self[i].raw_generators().iter().map(|x| self.q_deg_of(x))
        ).range().unwrap_or(0..=-1)
    }

    pub fn canon_cycles(&self) -> &[KhIChain<R>] { 
        &self.canon_cycles
    }

    pub fn truncated(&self, range: RangeInclusive<isize>) -> Self {
        Self::new_impl(
            self.inner.truncated(range),
            self.canon_cycles.clone(),
            self.deg_shift,
        )
    }

    pub fn homology(&self) -> KhIHomology<R>
    where R: EucRing, for<'x> &'x R: EucRingOps<R> {
        KhIHomology::from(self)
    }

    fn cached_bigraded(&self) -> &GrMod2<KhIGen, R> {
        self.cache_bigr.get_or_init(|| self.bigraded())
    }
}

impl<R> Bigraded<KhIGen, R> for KhIComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn base(&self) -> &GrMod1<KhIGen, R> { self.inner.summands() }
    fn decomp_key(&self, z: &KhIChain<R>) -> isize { self.q_deg_of_chain(z) }
}

impl<R> Index<isize> for KhIComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = KhIComplexSummand<R>;

    delegate! {
        to self.inner {
            fn index(&self, index: isize) -> &Self::Output;
        }
    }
}

impl<R> Index<(isize, isize)> for KhIComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = KhIComplexSummand<R>;

    fn index(&self, index: (isize, isize)) -> &Self::Output {
        &self.cached_bigraded()[index]
    }
}

impl<R> ToSeqString<isize> for KhIComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    delegate! {
        to self.inner { 
            fn label(&self) -> String;
            fn indices(&self) -> Vec<isize>;
            fn entry_at(&self, i: &isize) -> String;
        }
    }
}

impl<R> ToTableString<isize> for KhIComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn labels(&self) -> (String, String) { 
        ("i".to_string(), "j".to_string())
    }

    fn indices(&self) -> (Vec<isize>, Vec<isize>) { 
        (self.h_range().collect(), self.q_range().step_by(2).collect())
    }

    fn entry_at(&self, i: &isize, j: &isize) -> String {
        if self[(*i, *j)].is_zero() {
            ".".to_string()
        } else {
            self[(*i, *j)].to_string()
        }
    }
}

#[cfg(test)]
mod tests {
    use yui_core::poly::Poly;
    use yui_core::num::FF2;
    use num_traits::{Zero, One};
    use super::*;

    // canon-cycle tests are invariant under v1/v2 algorithm choice
    macro_rules! canon_tests {
        ($build:expr) => {
            #[test]
            fn canon_fbn() {
                let l = InvLink::test_data("3_1");

                type R = FF2;
                let (h, t) = (R::one(), R::zero());
                let c = $build(&l, &h, &t, false);

                let zs = c.canon_cycles.clone();

                assert_eq!(zs.len(), 4);
                assert!(zs[0].homogeneous_value(|x| c.h_deg_of(x)) == Some(0));
                assert!(zs[1].homogeneous_value(|x| c.h_deg_of(x)) == Some(0));
                assert!(zs[2].homogeneous_value(|x| c.h_deg_of(x)) == Some(1));
                assert!(zs[3].homogeneous_value(|x| c.h_deg_of(x)) == Some(1));

                for (i, z) in zs.iter().enumerate() {
                    let i = (i / 2) as isize;
                    assert!(c.d(i, z).is_zero());
                }
            }

            #[test]
            fn canon_fbn_red() {
                let l = InvLink::test_data("3_1");

                type R = FF2;
                let (h, t) = (R::one(), R::zero());
                let c = $build(&l, &h, &t, true);

                let zs = c.canon_cycles.clone();

                assert_eq!(zs.len(), 2);
                assert!(zs[0].homogeneous_value(|x| c.h_deg_of(x)) == Some(0));
                assert!(zs[1].homogeneous_value(|x| c.h_deg_of(x)) == Some(1));

                for (i, z) in zs.iter().enumerate() {
                    let i = i as isize;
                    assert!(c.d(i, z).is_zero());
                }
            }

            #[test]
            fn canon_bn() {
                let l = InvLink::test_data("3_1");

                type R = FF2;
                type P = Poly<'H', R>;
                let (h, t) = (P::variable(), P::zero());
                let c = $build(&l, &h, &t, false);

                let zs = c.canon_cycles.clone();

                assert_eq!(zs.len(), 4);
                assert!(zs[0].homogeneous_value(|x| c.h_deg_of(x)) == Some(0));
                assert!(zs[1].homogeneous_value(|x| c.h_deg_of(x)) == Some(0));
                assert!(zs[2].homogeneous_value(|x| c.h_deg_of(x)) == Some(1));
                assert!(zs[3].homogeneous_value(|x| c.h_deg_of(x)) == Some(1));

                for (i, z) in zs.iter().enumerate() {
                    let i = (i / 2) as isize;
                    assert!(c.d(i, z).is_zero());
                }
            }

            #[test]
            fn canon_bn_red() {
                let l = InvLink::test_data("3_1");

                type R = FF2;
                type P = Poly<'H', R>;
                let (h, t) = (P::variable(), P::zero());
                let c = $build(&l, &h, &t, true);

                let zs = c.canon_cycles.clone();

                assert_eq!(zs.len(), 2);
                assert!(zs[0].homogeneous_value(|x| c.h_deg_of(x)) == Some(0));
                assert!(zs[1].homogeneous_value(|x| c.h_deg_of(x)) == Some(1));

                for (i, z) in zs.iter().enumerate() {
                    let i = i as isize;
                    assert!(c.d(i, z).is_zero());
                }
            }
        };
    }

    mod v2 {
        use super::*;
        canon_tests!(KhIComplex::new);

        #[test]
        fn complex_kh() {
            let l = InvLink::test_data("3_1");

            type R = FF2;
            let (h, t) = (R::zero(), R::zero());
            let c = KhIComplex::new(&l, &h, &t, false);

            assert_eq!(c[0].rank(), 2);
            assert_eq!(c[1].rank(), 2);
            assert_eq!(c[2].rank(), 2);
            assert_eq!(c[3].rank(), 4);
            assert_eq!(c[4].rank(), 2);

            c.inner().check_d_all();
        }

        #[test]
        fn complex_fbn() {
            let l = InvLink::test_data("3_1");

            type R = FF2;
            let (h, t) = (R::one(), R::zero());
            let c = KhIComplex::new(&l, &h, &t, false);

            assert_eq!(c[0].rank(), 2);
            assert_eq!(c[1].rank(), 2);
            assert_eq!(c[2].rank(), 0);
            assert_eq!(c[3].rank(), 0);
            assert_eq!(c[4].rank(), 0);

            c.inner().check_d_all();
        }

        #[test]
        fn complex_bn() {
            let l = InvLink::test_data("3_1");

            type R = FF2;
            type P = Poly<'H', R>;
            let (h, t) = (P::variable(), P::zero());

            let c = KhIComplex::new(&l, &h, &t, false);

            assert_eq!(c[0].rank(), 2);
            assert_eq!(c[1].rank(), 2);
            assert_eq!(c[2].rank(), 2);
            assert_eq!(c[3].rank(), 4);
            assert_eq!(c[4].rank(), 2);

            c.inner().check_d_all();
        }

        #[test]
        fn complex_red() {
            let l = InvLink::test_data("3_1");

            type R = FF2;
            let (h, t) = (R::zero(), R::zero());
            let c = KhIComplex::new(&l, &h, &t, true);

            assert_eq!(c[0].rank(), 1);
            assert_eq!(c[1].rank(), 1);
            assert_eq!(c[2].rank(), 1);
            assert_eq!(c[3].rank(), 2);
            assert_eq!(c[4].rank(), 1);

            c.inner().check_d_all();
        }

        #[test]
        fn complex_kh_bigr() {
            let l = InvLink::test_data("3_1");

            type R = FF2;
            let (h, t) = (R::zero(), R::zero());
            let c = KhIComplex::new(&l, &h, &t, false);

            assert_eq!(c[(0, 1)].rank(), 1);
            assert_eq!(c[(0, 3)].rank(), 1);
            assert_eq!(c[(1, 1)].rank(), 1);
            assert_eq!(c[(1, 3)].rank(), 1);
            assert_eq!(c[(2, 5)].rank(), 1);
            assert_eq!(c[(2, 7)].rank(), 1);
            assert_eq!(c[(3, 5)].rank(), 1);
            assert_eq!(c[(3, 7)].rank(), 2);
            assert_eq!(c[(3, 9)].rank(), 1);
            assert_eq!(c[(4, 7)].rank(), 1);
            assert_eq!(c[(4, 9)].rank(), 1);
        }

        #[test]
        fn complex_kh_red_bigr() {
            let l = InvLink::test_data("3_1");

            type R = FF2;
            let (h, t) = (R::zero(), R::zero());
            let c = KhIComplex::new(&l, &h, &t, true);

            assert_eq!(c[(0, 2)].rank(), 1);
            assert_eq!(c[(1, 2)].rank(), 1);
            assert_eq!(c[(2, 6)].rank(), 1);
            assert_eq!(c[(3, 6)].rank(), 1);
            assert_eq!(c[(3, 8)].rank(), 1);
            assert_eq!(c[(4, 8)].rank(), 1);
        }

        // the cobordism-level cone must give the same bigraded homology as the matrix cone.
        #[test]
        fn cone_cob_matches_matrix() {
            use yui_homology::isize2;

            type R = FF2;
            let (h, t) = (R::zero(), R::zero());

            let nonzero = |m: &GrMod2<KhIGen, R>| -> Vec<(isize2, usize)> {
                m.support().map(|&k| (k, m[k].rank())).filter(|(_, r)| *r > 0).sorted().collect()
            };

            for name in ["3_1", "4_1", "6_3"] {
                for reduced in [false, true] {
                    let l = InvLink::test_data(name);
                    let matrix = KhIComplex::new(&l, &h, &t, reduced).homology().bigraded();
                    let config = SymBuildConfig { cone_cob: true, ..Default::default() };
                    let cone = KhIComplex::new_with_config(&l, &h, &t, reduced, config).homology().bigraded();
                    assert_eq!(nonzero(&matrix), nonzero(&cone), "{name} reduced={reduced}");
                }
            }
        }
    }

    mod v1 {
        use super::*;
        canon_tests!(KhIComplex::new_no_simplify);

        #[test]
        fn complex_kh() {
            let l = InvLink::test_data("3_1");

            type R = FF2;
            let (h, t) = (R::zero(), R::zero());
            let c = KhIComplex::new_no_simplify(&l, &h, &t, false);

            assert_eq!(c[0].rank(), 4);
            assert_eq!(c[1].rank(), 10);
            assert_eq!(c[2].rank(), 18);
            assert_eq!(c[3].rank(), 20);
            assert_eq!(c[4].rank(), 8);

            c.inner().check_d_all();
        }

        #[test]
        fn complex_fbn() {
            let l = InvLink::test_data("3_1");

            type R = FF2;
            let (h, t) = (R::one(), R::zero());
            let c = KhIComplex::new_no_simplify(&l, &h, &t, false);

            assert_eq!(c[0].rank(), 4);
            assert_eq!(c[1].rank(), 10);
            assert_eq!(c[2].rank(), 18);
            assert_eq!(c[3].rank(), 20);
            assert_eq!(c[4].rank(), 8);

            c.inner().check_d_all();
        }

        #[test]
        fn complex_bn() {
            let l = InvLink::test_data("3_1");

            type R = FF2;
            type P = Poly<'H', R>;
            let (h, t) = (P::variable(), P::zero());

            let c = KhIComplex::new_no_simplify(&l, &h, &t, false);

            assert_eq!(c[0].rank(), 4);
            assert_eq!(c[1].rank(), 10);
            assert_eq!(c[2].rank(), 18);
            assert_eq!(c[3].rank(), 20);
            assert_eq!(c[4].rank(), 8);

            c.inner().check_d_all();
        }

        #[test]
        fn complex_red() {
            let l = InvLink::test_data("3_1");

            type R = FF2;
            let (h, t) = (R::zero(), R::zero());
            let c = KhIComplex::new_no_simplify(&l, &h, &t, true);

            assert_eq!(c[0].rank(), 2);
            assert_eq!(c[1].rank(), 5);
            assert_eq!(c[2].rank(), 9);
            assert_eq!(c[3].rank(), 10);
            assert_eq!(c[4].rank(), 4);

            c.inner().check_d_all();
        }
    }
}