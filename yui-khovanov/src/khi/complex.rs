use std::ops::{Index, RangeInclusive};
use cartesian::cartesian;
use delegate::delegate;

use itertools::Itertools;
use yui::lc::Lc;
use yui::{EucRing, EucRingOps, Ring, RingOps};
use yui_homology::{isize2, ChainComplexTrait, Grid1, Grid2, GridTrait, ChainComplex, ChainComplex2, Summand};
use yui_link::InvLink;
use yui_matrix::sparse::SpMat;

use crate::kh::{KhChain, KhChainExt, KhComplex, KhGen};
use crate::khi::KhIHomology;
use crate::khi::KhIGen;
use crate::misc::range_of;

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
    deg_shift: (isize, isize)
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

        let deg_shift = KhComplex::deg_shift_for(l.link(), reduced);
        let cube = KhICube::new(l, h, t, reduced, deg_shift);
        let inner = cube.into_complex();

        let canon_cycles = if l.base_pt().is_some() && l.link().is_knot() {
            let p = l.base_pt().unwrap();
            let zs = KhComplex::make_canon_cycles(l.link(), p, &R::zero(), h, reduced, deg_shift);
            Iterator::chain(
                zs.iter().map(|z| z.map_gens(|x| KhIGen::B(*x))),
                zs.iter().map(|z| z.map_gens(|x| KhIGen::Q(*x)))
            ).collect()
        } else { 
            vec![]
        };

        Self::new_impl(inner, canon_cycles, deg_shift)
    }

    pub fn from_kh_complex<'a, F>(c: KhComplex<R>, map: F) -> Self
    where F: Fn(&KhGen) -> KhGen + Send + Sync + 'static {
        let deg_shift = c.deg_shift();
        let h_range = c.h_range();
        let h_range = *h_range.start() ..= (h_range.end() + 1);

        let canon_cycles = c.canon_cycles().iter().flat_map(|z| { 
            let bz = z.map_gens(|x| KhIGen::B(*x));
            let qz = z.map_gens(|x| KhIGen::Q(*x));
            [bz, qz]
        }).sorted_by_key(|z| z.h_deg()).collect_vec();

        let summands = Grid1::generate(h_range, |i| { 
            let b_gens = c[i].raw_gens().iter().map(|x| KhIGen::B(*x));
            let q_gens = c[i - 1].raw_gens().iter().map(|x| KhIGen::Q(*x));
            Summand::from_raw_gens(Iterator::chain(b_gens, q_gens))
        });

        let d = move |i: isize, x: &KhIGen| -> KhIChain<R> { 
            match x { 
                KhIGen::B(x) => {
                    let z = KhChain::from(*x);
                    let dx = c.d(i, &z).map_gens(|y| KhIGen::B(*y));
                    let qx = KhIChain::from(KhIGen::Q(*x));
                    let qtx = {
                        let tx = map(x);
                        KhIChain::from(KhIGen::Q(tx))
                    };
                    dx + qx + qtx
                },
                KhIGen::Q(x) => {
                    let z = KhChain::from(*x);
                    c.d(i, &z).map_gens(|y| KhIGen::Q(*y))
                }
            }
        };

        let inner = ChainComplex::new(summands, 1, move |i, z| { 
            z.apply(|x| d(i, x))
        });

        KhIComplex::new_impl(inner, canon_cycles, deg_shift)
    }

    pub(crate) fn new_impl(inner: ChainComplex<KhIGen, R>, canon_cycles: Vec<KhIChain<R>>, deg_shift: (isize, isize)) -> Self { 
        Self { inner, canon_cycles, deg_shift }
    }

    pub fn inner(&self) -> &ChainComplex<KhIGen, R> {
        &self.inner
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        range_of(self.support())
    }

    pub fn q_range(&self) -> RangeInclusive<isize> {
        range_of(self.support().flat_map(|i| 
            self[i].raw_gens().iter().map(|x| x.q_deg())
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

    pub fn gen_grid(&self) -> Grid2<Summand<KhIGen, R>> {
        let h_range = self.h_range();
        let q_range = self.q_range().step_by(2);
        let support = cartesian!(h_range, q_range.clone()).map(|(i, j)| 
            isize2(i, j)
        );

        Grid2::generate(support, |idx| { 
            let isize2(i, j) = idx;
            let gens = self[i].raw_gens().iter().filter(|x| { 
                x.q_deg() == j
            }).cloned();
            Summand::from_raw_gens(gens)
        })
    }

    pub fn into_bigraded(self) -> ChainComplex2<KhIGen, R> {
        // TODO assert h == 0
        let summands = self.gen_grid();

        ChainComplex2::new(summands, isize2(1, 0), move |idx, x| { 
            let i = idx.0;
            let x = KhIChain::from(x.clone());
            let dx = self.d(i, &x);
            dx.into_iter().collect()
        })
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

impl<R> GridTrait<isize> for KhIComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Support = std::vec::IntoIter<isize>;
    type Item = KhIComplexSummand<R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Support;
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

#[cfg(test)]
mod tests {
    use yui::poly::HPoly;
    use yui::FF2;
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
        let c = KhIComplex::new(&l, &h, &t, false).into_bigraded();

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

        c.check_d_all();
    }

    #[test]
    fn complex_kh_red_bigr() {
        let l = InvLink::load("3_1").unwrap();

        type R = FF2;
        let (h, t) = (R::zero(), R::zero());
        let c = KhIComplex::new(&l, &h, &t, true).into_bigraded();

        assert_eq!(c[(0, 2)].rank(), 1);
        assert_eq!(c[(1, 2)].rank(), 1);
        assert_eq!(c[(2, 6)].rank(), 1);
        assert_eq!(c[(3, 6)].rank(), 1);
        assert_eq!(c[(3, 8)].rank(), 1);
        assert_eq!(c[(4, 8)].rank(), 1);
        
        c.check_d_all();
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
    use yui::poly::HPoly;
    use yui::FF2;
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