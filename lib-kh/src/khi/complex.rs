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
use yui_core::abst::{EucRing, EucRingOps, Ring, RingOps};
use yui_core::ext::{empty_range, IteratorExt};
use yui_homology::{ChainComplex1, ChainMap, ToSeqString, ToTableString, GrMod1, GrMod2, Summand};
use yui_link::InvLink;

use crate::kh::{KhComplex, KhGen};
use crate::tng::builder::{SymBuildConfig, assert_supported_symmetry};
use crate::khi::KhIHomology;
use crate::khi::{KhIGen, KhIGenExt};
use crate::util::Bigraded;
use yui_core::util::tex::TeX;
use yui_homology::tex::{ToTexSeq, ToTexTable};

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

    /// The default KhI construction: the cobordism-level cone (`ConeBuilder`) yields the coned complex
    /// + canon classes directly, and `into_raw_complex` converts once at the boundary (matrix-backed).
    /// The equivalent matrix-level cone is kept for reference as `new_with_config_v1`.
    pub fn new_with_config(l: &InvLink, h: &R, t: &R, reduced: bool, config: SymBuildConfig) -> Self {
        assert_supported_symmetry(l);

        let (h_range, build_config) = Self::cone_build_config(l, reduced, config);
        Self::build_cone(l, h, t, reduced, build_config, h_range)
    }

    // The ssi (`V2`) entry: `config.h_range` is used literally for the build (the driver pre-widens
    // it), and only `raw_range` is converted — the solves never touch the other degrees.
    pub(crate) fn new_windowed(l: &InvLink, h: &R, t: &R, reduced: bool, config: SymBuildConfig, raw_range: RangeInclusive<isize>) -> Self {
        Self::build_cone(l, h, t, reduced, config, Some(raw_range))
    }

    fn build_cone(l: &InvLink, h: &R, t: &R, reduced: bool, build_config: SymBuildConfig, raw_range: Option<RangeInclusive<isize>>) -> Self {
        use crate::tng::builder::ConeBuilder;
        assert_eq!(R::one() + R::one(), R::zero(), "char(R) != 2"); // the cobordism cone is char-2 only

        let cone = ConeBuilder::from_inv_link(l, h, t, reduced).with_config(build_config).run();

        // sort by h-degree (all `B` then all `Q`) to match the matrix cone's canon-cycle order.
        let canon_cycles = cone.eval_khi_elements().into_iter()
            .sorted_by_key(|z| z.keys().map(|x| x.rel_h_deg()).min().unwrap_or(0))
            .collect_vec();

        let inner = cone.into_raw_complex(raw_range);

        let deg_shift = KhComplex::<R>::deg_shift_for(l.inner(), reduced);
        Self::new_impl(inner, canon_cycles, deg_shift)
    }

    pub fn new_no_simplify(l: &InvLink, h: &R, t: &R, reduced: bool) -> Self {
        assert_eq!(R::one() + R::one(), R::zero(), "char(R) != 2");
        assert_supported_symmetry(l);
        assert!(
            !reduced || (l.base_pt().is_some_and(|e| l.is_on_axis(e)) && t.is_zero()),
            "reduced requires t = 0 and a base point on the axis"
        );

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
        self.support().copied().range().unwrap_or_else(empty_range)
    }

    pub fn q_range(&self) -> RangeInclusive<isize> {
        self.support().flat_map(|&i|
            self[i].raw_generators().iter().map(|x| self.q_deg_of(x))
        ).range().unwrap_or_else(empty_range)
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

    // `config.h_range` is the desired cone range `[a, b]`, clamped; since `KhI_i = C_i ⊕ C_{i-1}`,
    // the build gets `[a-1, b]` while the rest of `config` is kept.
    fn cone_build_config(l: &InvLink, reduced: bool, config: SymBuildConfig) -> (Option<RangeInclusive<isize>>, SymBuildConfig) {
        let h_range = config.h_range.clone().map(|r| KhComplex::<R>::clamp_h_range(l.inner(), reduced, r));
        let build_config = SymBuildConfig {
            h_range: h_range.as_ref().map(|r| (*r.start() - 1) ..= *r.end()),
            ..config
        };
        (h_range, build_config)
    }
}

// Reference: the "honest" matrix-level cone — build the sym `KhComplex`, then cone `(1+τ)` at the
// matrix level (`cone_of`). Superseded by the cobordism cone in `new_with_config`; kept for
// cross-checks (see the `cone_canon_ssi_matches_matrix` test).
#[allow(dead_code)]
impl<R> KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn new_with_config_v1(l: &InvLink, h: &R, t: &R, reduced: bool, config: SymBuildConfig) -> Self {
        use crate::tng::builder::SymTngBuilder;

        let (h_range, build_config) = Self::cone_build_config(l, reduced, config);
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

impl<R> ToTexSeq<isize> for KhIComplex<R>
where R: Ring + TeX, for<'x> &'x R: RingOps<R> {
    fn tex_entry_at(&self, i: &isize) -> String {
        if self[*i].is_zero() {
            ".".to_string()
        } else {
            self[*i].tex_string()
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

impl<R> ToTexTable<isize> for KhIComplex<R>
where R: Ring + TeX, for<'x> &'x R: RingOps<R> {
    fn tex_entry_at(&self, i: &isize, j: &isize) -> String {
        if self[(*i, *j)].is_zero() {
            ".".to_string()
        } else {
            self[(*i, *j)].tex_string()
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

        // the default cobordism cone (`new`) must give the same bigraded homology as the honest
        // matrix cone (`new_with_config_v1`).
        #[test]
        fn cone_matches_v1() {
            use yui_homology::isize2;

            type R = FF2;
            let (h, t) = (R::zero(), R::zero());

            let nonzero = |m: &GrMod2<KhIGen, R>| -> Vec<(isize2, usize)> {
                m.support().map(|&k| (k, m[k].rank())).filter(|(_, r)| *r > 0).sorted().collect()
            };

            for name in ["3_1", "4_1", "6_3"] {
                for reduced in [false, true] {
                    let l = InvLink::test_data(name);
                    let v1 = KhIComplex::new_with_config_v1(&l, &h, &t, reduced, SymBuildConfig::default()).homology().bigraded();
                    let v2 = KhIComplex::new(&l, &h, &t, reduced).homology().bigraded();
                    assert_eq!(nonzero(&v1), nonzero(&v2), "{name} reduced={reduced}");
                }
            }
        }

        // The sym+cone q-filter: over Khovanov (d preserves q), a q-window keeps exactly the
        // in-window generators, and bigraded KhI homology at each kept (i, q) is unchanged.
        // Runs the single-pass build (final-merge prune + finalize deloop filter) and a chunked
        // build (sym chunk-merges too).
        #[test]
        fn q_filter_matches_full() {
            use crate::tng::builder::CutOption;
            type R = FF2;
            let l = InvLink::test_data("6_3");
            let (h, t) = (R::zero(), R::zero());

            let full = KhIComplex::new(&l, &h, &t, false);
            let (h_range, q_range) = (full.h_range(), full.q_range());
            let full_h = full.homology();

            let (lo, hi) = (*q_range.start() + 2, *q_range.end() - 2);

            for cut in [CutOption::None, CutOption::Auto(3)] {
                let config = SymBuildConfig { q_range: Some(lo ..= hi), cut: cut.clone(), ..Default::default() };
                let win = KhIComplex::new_with_config(&l, &h, &t, false, config);

                for i in win.h_range() {
                    for x in win[i].raw_generators() {
                        assert!((lo ..= hi).contains(&win.q_deg_of(x)), "gen out of window: ({i}, {}), cut={cut:?}", win.q_deg_of(x));
                    }
                }

                let win_h = win.homology();
                for i in h_range.clone() {
                    for q in q_range.clone().step_by(2) {
                        let expected = if (lo ..= hi).contains(&q) { full_h[(i, q)].rank() } else { 0 };
                        assert_eq!(win_h[(i, q)].rank(), expected, "rank ({i}, {q}), cut={cut:?}");
                    }
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

        // an off-axis base point leaves the reduced complex without the τ-images of its
        // generators, which `vectorize` then drops without a word.
        fn based_off_axis() -> InvLink {
            let l = InvLink::test_data("3_1");
            let e_map: Vec<_> = l.inner().edges().into_iter().map(|e| (e, l.inv_edge(e))).collect();
            let off = l.inner().edges().into_iter().find(|&e| !l.is_on_axis(e)).unwrap();
            InvLink::new(l.inner().clone().with_base_pt(off), e_map)
        }

        #[test]
        #[should_panic(expected = "base point on the axis")]
        fn reduced_rejects_off_axis_base_pt() {
            let (h, t) = (FF2::zero(), FF2::zero());
            let _ = KhIComplex::new_no_simplify(&based_off_axis(), &h, &t, true);
        }

        #[test]
        #[should_panic(expected = "base point on the axis")]
        fn reduced_rejects_off_axis_base_pt_cone() {
            let (h, t) = (FF2::zero(), FF2::zero());
            let _ = KhIComplex::new(&based_off_axis(), &h, &t, true);
        }
    }
}