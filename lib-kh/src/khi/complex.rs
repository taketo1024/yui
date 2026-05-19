use std::ops::{Index, RangeInclusive};
use std::sync::OnceLock;
use delegate::delegate;

use itertools::Itertools;
use yui_core::lc::Lc;
use yui_core::{EucRing, EucRingOps, Ring, RingOps};
use yui_homology::{ChainComplex, ChainComplexTrait, DisplaySeq, DisplayTable, Grid1, Grid2, GridIter, GridTrait, Summand, SummandTrait};
use yui_link::InvLink;
use yui_matrix::sparse::SpMat;

use crate::kh::{KhChain, KhChainExt, KhComplex, KhChainGen};
use crate::khi::KhIHomology;
use crate::khi::KhIGen;
use crate::misc::{make_gen_grid, range_of};

pub type KhIChain<R> = Lc<KhIGen, R>;

impl<R> KhChainExt for KhIChain<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn h_deg(&self) -> isize {
        self.gens().map(|x| x.h_deg()).min().unwrap_or(0)
    }
    
    fn q_deg(&self) -> isize {
        self.gens().map(|x| x.q_deg()).min().unwrap_or(0)
    }
}

pub type KhIComplexSummand<R> = Summand<KhIGen, R>;

#[derive(Clone)]
pub struct KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    inner: ChainComplex<KhIGen, R>,
    canon_cycles: Vec<KhIChain<R>>,
    deg_shift: (isize, isize),
    gen_grid: OnceLock<Grid2<KhIComplexSummand<R>>>,
}

impl<R> KhIComplex<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn new(l: &InvLink, h: &R, t: &R, reduced: bool) -> Self { 
        use crate::khi::internal::v2::builder::SymTngBuilder;

        SymTngBuilder::build_khi_complex(l, h, t, reduced)
    }

    pub fn new_no_simplify(l: &InvLink, h: &R, t: &R, reduced: bool) -> Self { 
        use crate::khi::internal::v1::cube::KhICube;
        use crate::kh::KhComplex;

        assert_eq!(R::one() + R::one(), R::zero(), "char(R) != 2");
        assert!(!reduced || (l.base_pt().is_some() && t.is_zero()));

        let deg_shift = KhComplex::deg_shift_for(l.inner(), reduced);

        // TODO use mapping cone

        let cube = KhICube::new(l, h, t, reduced, deg_shift);
        let inner = cube.into_complex();

        let canon_cycles = if l.base_pt().is_some() && l.is_knot() {
            let p = l.base_pt().unwrap();
            let zs = KhComplex::make_canon_cycles(l.inner(), p, &R::zero(), h, reduced, deg_shift);
            Iterator::chain(
                zs.iter().map(|z| z.clone().map_gens(|x| KhIGen::B(x))),
                zs.iter().map(|z| z.clone().map_gens(|x| KhIGen::Q(x)))
            ).collect()
        } else { 
            vec![]
        };

        Self::new_impl(inner, canon_cycles, deg_shift)
    }

    pub fn from_kh_complex<'a, F>(c: KhComplex<R>, map: F) -> Self
    where F: Fn(&KhChainGen) -> KhChainGen + Send + Sync + 'static {
        let deg_shift = c.deg_shift();
        let h_range = c.h_range();
        let h_range = *h_range.start() ..= (h_range.end() + 1);

        let canon_cycles = c.canon_cycles().iter().flat_map(|z| { 
            let bz = z.clone().map_gens(|x| KhIGen::B(x));
            let qz = z.clone().map_gens(|x| KhIGen::Q(x));
            [bz, qz]
        }).sorted_by_key(|z| z.h_deg()).collect_vec();

        // TODO use mapping cone

        let summands = Grid1::generate(h_range, |i| { 
            let b_gens = c[i].raw_generators().iter().map(|x| KhIGen::B(*x));
            let q_gens = c[i - 1].raw_generators().iter().map(|x| KhIGen::Q(*x));
            Summand::from_raw_generators(Iterator::chain(b_gens, q_gens))
        });

        let d = move |i: isize, x: &KhIGen| -> KhIChain<R> { 
            match x { 
                KhIGen::B(x) => {
                    let z = KhChain::from(*x);
                    let dx = c.d(i, &z).map_gens(|y| KhIGen::B(y));
                    let qx = KhIChain::from(KhIGen::Q(*x));
                    let qtx = {
                        let tx = map(x);
                        KhIChain::from(KhIGen::Q(tx))
                    };
                    dx + qx + qtx
                },
                KhIGen::Q(x) => {
                    let z = KhChain::from(*x);
                    c.d(i, &z).map_gens(|y| KhIGen::Q(y))
                }
            }
        };

        let inner = ChainComplex::new(summands, 1, move |i, z| { 
            z.apply(|x| d(i, x))
        });

        KhIComplex::new_impl(inner, canon_cycles, deg_shift)
    }

    pub(crate) fn new_impl(inner: ChainComplex<KhIGen, R>, canon_cycles: Vec<KhIChain<R>>, deg_shift: (isize, isize)) -> Self {
        Self { inner, canon_cycles, deg_shift, gen_grid: OnceLock::new() }
    }

    pub fn inner(&self) -> &ChainComplex<KhIGen, R> {
        &self.inner
    }

    pub fn h_range(&self) -> RangeInclusive<isize> {
        range_of(self.support().copied())
    }

    pub fn q_range(&self) -> RangeInclusive<isize> {
        range_of(self.support().flat_map(|&i|
            self[i].raw_generators().iter().map(|x| x.q_deg())
        ))
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

    fn gen_grid(&self) -> &Grid2<KhIComplexSummand<R>> {
        self.gen_grid.get_or_init(|| make_gen_grid(self.inner.summands()))
    }

    pub fn homology(&self) -> KhIHomology<R>
    where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
        KhIHomology::from(self)
    }
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

    delegate! {
        to self.gen_grid() {
            fn index(&self, index: (isize, isize)) -> &Self::Output;
        }
    }
}

impl<R> GridTrait<isize> for KhIComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Item = KhIComplexSummand<R>;
    type Support<'a> = GridIter<'a, isize, Self::Item> where Self: 'a, R: 'a;

    delegate! {
        to self.inner {
            fn support(&self) -> Self::Support<'_>;
            fn is_supported(&self, i: isize) -> bool;
            fn get(&self, i: isize) -> &Self::Item;
            fn get_default(&self) -> &Self::Item;
        }
    }
}

impl<R> ChainComplexTrait<isize> for KhIComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;
    type Element = KhIChain<R>;

    delegate! { 
        to self.inner { 
            fn rank(&self, i: isize) -> usize;
            fn d_deg(&self) -> isize;
            fn d(&self, i: isize, z: &Self::Element) -> Self::Element;
            fn d_matrix(&self, i: isize) -> SpMat<R>;
        }
    }
}

impl<R> DisplaySeq<isize> for KhIComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    delegate! {
        to self.inner { 
            fn display_label(&self) -> String;
            fn display_indices(&self) -> Vec<isize>;
            fn display_at(&self, i: &isize) -> String;
        }
    }
}

impl<R> DisplayTable<isize> for KhIComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn display_labels(&self) -> (String, String) { 
        ("i".to_string(), "j".to_string())
    }

    fn display_indices(&self) -> (Vec<isize>, Vec<isize>) { 
        (self.h_range().into_iter().collect(), self.q_range().step_by(2).collect())
    }

    fn display_at(&self, i: &isize, j: &isize) -> String {
        if self[(*i, *j)].is_zero() {
            ".".to_string()
        } else {
            self[(*i, *j)].to_string()
        }
    }
}

#[cfg(test)]
mod tests {
    use yui_core::poly::HPoly;
    use yui_core::num::FF2;
    use num_traits::{Zero, One};
    use yui_homology::{ChainComplexTrait, SummandTrait};
    use super::*;

    #[test]
    fn complex_kh() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        let (h, t) = (R::zero(), R::zero());
        let c = KhIComplex::new(&l, &h, &t, false);

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 2);
        assert_eq!(c[2].rank(), 2);
        assert_eq!(c[3].rank(), 4);
        assert_eq!(c[4].rank(), 2);
            
        c.check_d_all();
    }

    #[test]
    fn complex_fbn() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        let (h, t) = (R::one(), R::zero());
        let c = KhIComplex::new(&l, &h, &t, false);

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 2);
        assert_eq!(c[2].rank(), 0);
        assert_eq!(c[3].rank(), 0);
        assert_eq!(c[4].rank(), 0);
        
        c.check_d_all();
    }

    #[test]
    fn complex_bn() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        type P = HPoly<'H', R>;
        let (h, t) = (P::variable(), P::zero());

        let c = KhIComplex::new(&l, &h, &t, false);

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 2);
        assert_eq!(c[2].rank(), 2);
        assert_eq!(c[3].rank(), 4);
        assert_eq!(c[4].rank(), 2);
        
        c.check_d_all();
    }

    #[test]
    fn complex_red() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        let (h, t) = (R::zero(), R::zero());
        let c = KhIComplex::new(&l, &h, &t, true);

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 1);
        assert_eq!(c[2].rank(), 1);
        assert_eq!(c[3].rank(), 2);
        assert_eq!(c[4].rank(), 1);
        
        c.check_d_all();
    }

    #[test]
    fn complex_kh_bigr() { 
        let l = InvLink::load("3_1").unwrap();

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
        let l = InvLink::load("3_1").unwrap();

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

    #[test]
    fn canon_fbn() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        let (h, t) = (R::one(), R::zero());
        let c = KhIComplex::new(&l, &h, &t, false);

        let zs = c.canon_cycles.clone();

        assert_eq!(zs.len(), 4);
        assert!(zs[0].gens().all(|x| x.h_deg() == 0));
        assert!(zs[1].gens().all(|x| x.h_deg() == 0));
        assert!(zs[2].gens().all(|x| x.h_deg() == 1));
        assert!(zs[3].gens().all(|x| x.h_deg() == 1));

        for (i, z) in zs.iter().enumerate() { 
            let i = (i / 2) as isize;
            assert!(c.d(i, z).is_zero());
        }
    }

    #[test]
    fn canon_fbn_red() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        let (h, t) = (R::one(), R::zero());
        let c = KhIComplex::new(&l, &h, &t, true);

        let zs = c.canon_cycles.clone();

        assert_eq!(zs.len(), 2);
        assert!(zs[0].gens().all(|x| x.h_deg() == 0));
        assert!(zs[1].gens().all(|x| x.h_deg() == 1));

        for (i, z) in zs.iter().enumerate() { 
            let i = i as isize;
            assert!(c.d(i, z).is_zero());
        }
    }

    #[test]
    fn canon_bn() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        type P = HPoly<'H', R>;
        let (h, t) = (P::variable(), P::zero());
        let c = KhIComplex::new(&l, &h, &t, false);

        let zs = c.canon_cycles.clone();

        assert_eq!(zs.len(), 4);
        assert!(zs[0].gens().all(|x| x.h_deg() == 0));
        assert!(zs[1].gens().all(|x| x.h_deg() == 0));
        assert!(zs[2].gens().all(|x| x.h_deg() == 1));
        assert!(zs[3].gens().all(|x| x.h_deg() == 1));

        for (i, z) in zs.iter().enumerate() { 
            let i = (i / 2) as isize;
            assert!(c.d(i, z).is_zero());
        }
    }

    #[test]
    fn canon_bn_red() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        type P = HPoly<'H', R>;
        let (h, t) = (P::variable(), P::zero());
        let c = KhIComplex::new(&l, &h, &t, true);
        
        let zs = c.canon_cycles.clone();

        assert_eq!(zs.len(), 2);
        assert!(zs[0].gens().all(|x| x.h_deg() == 0));
        assert!(zs[1].gens().all(|x| x.h_deg() == 1));

        for (i, z) in zs.iter().enumerate() { 
            let i = i as isize;
            assert!(c.d(i, z).is_zero());
        }
    }
}

#[cfg(test)]
mod tests_v1 {
    use yui_core::poly::HPoly;
    use yui_core::num::FF2;
    use num_traits::{Zero, One};
    use yui_homology::{ChainComplexTrait, SummandTrait};
    use super::*;

    #[test]
    fn complex_kh() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        let (h, t) = (R::zero(), R::zero());
        let c = KhIComplex::new_no_simplify(&l, &h, &t, false);

        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 10);
        assert_eq!(c[2].rank(), 18);
        assert_eq!(c[3].rank(), 20);
        assert_eq!(c[4].rank(), 8);
            
        c.check_d_all();
    }

    #[test]
    fn complex_fbn() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        let (h, t) = (R::one(), R::zero());
        let c = KhIComplex::new_no_simplify(&l, &h, &t, false);

        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 10);
        assert_eq!(c[2].rank(), 18);
        assert_eq!(c[3].rank(), 20);
        assert_eq!(c[4].rank(), 8);
        
        c.check_d_all();
    }

    #[test]
    fn complex_bn() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        type P = HPoly<'H', R>;
        let (h, t) = (P::variable(), P::zero());

        let c = KhIComplex::new_no_simplify(&l, &h, &t, false);

        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 10);
        assert_eq!(c[2].rank(), 18);
        assert_eq!(c[3].rank(), 20);
        assert_eq!(c[4].rank(), 8);
        
        c.check_d_all();
    }

    #[test]
    fn complex_red() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        let (h, t) = (R::zero(), R::zero());
        let c = KhIComplex::new_no_simplify(&l, &h, &t, true);

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 5);
        assert_eq!(c[2].rank(), 9);
        assert_eq!(c[3].rank(), 10);
        assert_eq!(c[4].rank(), 4);
        
        c.check_d_all();
    }

    #[test]
    fn canon_fbn() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        let (h, t) = (R::one(), R::zero());
        let c = KhIComplex::new_no_simplify(&l, &h, &t, false);

        let zs = c.canon_cycles.clone();

        assert_eq!(zs.len(), 4);
        assert!(zs[0].gens().all(|x| x.h_deg() == 0));
        assert!(zs[1].gens().all(|x| x.h_deg() == 0));
        assert!(zs[2].gens().all(|x| x.h_deg() == 1));
        assert!(zs[3].gens().all(|x| x.h_deg() == 1));

        for (i, z) in zs.iter().enumerate() { 
            let i = (i / 2) as isize;
            assert!(c.d(i, z).is_zero());
        }
    }

    #[test]
    fn canon_fbn_red() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        let (h, t) = (R::one(), R::zero());
        let c = KhIComplex::new_no_simplify(&l, &h, &t, true);

        let zs = c.canon_cycles.clone();

        assert_eq!(zs.len(), 2);
        assert!(zs[0].gens().all(|x| x.h_deg() == 0));
        assert!(zs[1].gens().all(|x| x.h_deg() == 1));

        for (i, z) in zs.iter().enumerate() { 
            let i = i as isize;
            assert!(c.d(i, z).is_zero());
        }
    }

    #[test]
    fn canon_bn() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        type P = HPoly<'H', R>;
        let (h, t) = (P::variable(), P::zero());
        let c = KhIComplex::new_no_simplify(&l, &h, &t, false);

        let zs = c.canon_cycles.clone();

        assert_eq!(zs.len(), 4);
        assert!(zs[0].gens().all(|x| x.h_deg() == 0));
        assert!(zs[1].gens().all(|x| x.h_deg() == 0));
        assert!(zs[2].gens().all(|x| x.h_deg() == 1));
        assert!(zs[3].gens().all(|x| x.h_deg() == 1));

        for (i, z) in zs.iter().enumerate() { 
            let i = (i / 2) as isize;
            assert!(c.d(i, z).is_zero());
        }
    }

    #[test]
    fn canon_bn_red() { 
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        type P = HPoly<'H', R>;
        let (h, t) = (P::variable(), P::zero());
        let c = KhIComplex::new_no_simplify(&l, &h, &t, true);
        
        let zs = c.canon_cycles.clone();

        assert_eq!(zs.len(), 2);
        assert!(zs[0].gens().all(|x| x.h_deg() == 0));
        assert!(zs[1].gens().all(|x| x.h_deg() == 1));

        for (i, z) in zs.iter().enumerate() { 
            let i = i as isize;
            assert!(c.d(i, z).is_zero());
        }
    }
}