use std::ops::{Index, RangeInclusive};
use std::sync::OnceLock;
use delegate::delegate;

use itertools::{Either, Itertools};
use yui_core::lc::Lc;
use yui_core::{EucRing, EucRingOps, IteratorExt, Ring, RingOps};
use yui_homology::{ChainComplex1, ToSeqString, ToTableString, GrMod1, GrMod2, Summand};
use yui_link::InvLink;

use crate::kh::{KhChain, KhComplex, KhGen};
use crate::khi::KhIHomology;
use crate::khi::{KhIGen, KhIGenExt};
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
        use crate::khi::internal::v2::builder::SymTngBuilder;

        SymTngBuilder::build_khi_complex(l, h, t, reduced)
    }

    pub fn new_no_simplify(l: &InvLink, h: &R, t: &R, reduced: bool) -> Self {
        assert_eq!(R::one() + R::one(), R::zero(), "char(R) != 2");
        assert!(!reduced || (l.base_pt().is_some() && t.is_zero()));

        let c = KhComplex::new_no_simplify(l.inner(), h, t, reduced);
        Self::from_kh_complex(c, crate::khi::tau::tau_map(l))
    }

    pub fn from_kh_complex<'a, F>(c: KhComplex<R>, map: F) -> Self
    where F: Fn(&KhGen) -> KhGen + Send + Sync + 'static {
        let deg_shift = c.deg_shift();
        let h_range = c.h_range();
        let h_range = *h_range.start() ..= (h_range.end() + 1);

        let canon_cycles = c.canon_cycles().iter().flat_map(|z| { 
            let bz = z.clone().map_keys(KhIGen::from_left);
            let qz = z.clone().map_keys(KhIGen::from_right);
            [bz, qz]
        }).sorted_by_key(|z| z.keys().map(|x| x.rel_h_deg()).min().unwrap_or(0)).collect_vec();

        // TODO use mapping cone

        let summands = GrMod1::generate(h_range, |i| { 
            let b_gens = c[i].raw_generators().iter().map(|x| KhIGen::from_left(*x));
            let q_gens = c[i - 1].raw_generators().iter().map(|x| KhIGen::from_right(*x));
            Summand::from_raw_generators(Iterator::chain(b_gens, q_gens))
        });

        let d = move |i: isize, x: &KhIGen| -> KhIChain<R> {
            match x.inner() {
                Either::Left(x) => {
                    let z = KhChain::from(*x);
                    let dx = c.d(i, &z).map_keys(KhIGen::from_left);
                    let qx = KhIChain::from(KhIGen::from_right(*x));
                    let qtx = {
                        let tx = map(x);
                        KhIChain::from(KhIGen::from_right(tx))
                    };
                    dx + qx + qtx
                },
                Either::Right(x) => {
                    let z = KhChain::from(*x);
                    c.d(i, &z).map_keys(KhIGen::from_right)
                }
            }
        };

        let inner = ChainComplex1::new(summands, 1, move |i, z| { 
            z.apply(|x| d(i, x))
        });

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