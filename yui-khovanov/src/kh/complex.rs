use std::ops::{RangeInclusive, Index};
use cartesian::cartesian;

use delegate::delegate;
use yui::lc::Lc;
use yui::{Ring, RingOps, EucRing, EucRingOps};
use yui_link::Link;
use yui_homology::{isize2, ChainComplexTrait, Grid2, GridTrait, ChainComplex, ChainComplex2, Summand};
use yui_matrix::sparse::SpMat;

use crate::kh::{KhGen, KhHomology, KhHomologyBigraded};
use crate::misc::range_of;

use super::KhAlgStr;

pub type KhChain<R> = Lc<KhGen, R>;
pub trait KhChainExt { 
    fn h_deg(&self) -> isize;
    fn q_deg(&self) -> isize;
}

impl<R> KhChainExt for KhChain<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn h_deg(&self) -> isize {
        self.gens().map(|x| x.h_deg()).min().unwrap_or(0)
    }
    
    fn q_deg(&self) -> isize {
        self.gens().map(|x| x.q_deg()).min().unwrap_or(0)
    }
}

pub type KhComplexSummand<R> = Summand<KhGen, R>;

#[derive(Clone)]
pub struct KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    inner: ChainComplex<KhGen, R>,
    str: KhAlgStr<R>,
    deg_shift: (isize, isize),
    reduced: bool,
    canon_cycles: Vec<KhChain<R>>,
}

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(link: &Link, h: &R, t: &R, reduced: bool) -> Self {
        cfg_if::cfg_if! { 
        if #[cfg(feature = "old")] { 
            Self::new_v1(link, h, t, reduced)
        } else { 
            Self::new_v2(link, h, t, reduced)
        }}
    }

    #[cfg(not(feature = "old"))]
    fn new_v2(l: &Link, h: &R, t: &R, reduced: bool) -> Self { 
        use crate::kh::internal::v2::builder::TngComplexBuilder;

        assert!(!reduced || (!l.is_empty() && t.is_zero()));

        TngComplexBuilder::build_kh_complex(l, h, t, reduced)
    }

    #[cfg(feature = "old")]
    fn new_v1(l: &Link, h: &R, t: &R, reduced: bool) -> Self { 
        use crate::kh::internal::v1::cube::KhCube;

        assert!(!reduced || (!l.is_empty() && t.is_zero()));

        let red_e = reduced.then(|| l.first_edge().unwrap());
        let deg_shift = Self::deg_shift_for(l, reduced);
        
        let cube = KhCube::new(l, h, t, red_e, deg_shift);
        let str = cube.str().clone();
        let complex = cube.into_complex();

        let canon_cycles = if t.is_zero() && l.is_knot() {
            let p = l.first_edge().unwrap();
            Self::make_canon_cycles(l, p, &R::zero(), h, reduced, deg_shift)
        } else { 
            vec![]
        };

        KhComplex::new_impl(complex, str, deg_shift, reduced, canon_cycles)
    }

    pub(crate) fn new_impl(inner: ChainComplex<KhGen, R>, str: KhAlgStr<R>, deg_shift: (isize, isize), reduced: bool, canon_cycles: Vec<KhChain<R>>) -> Self { 
        KhComplex { inner, str, deg_shift, reduced, canon_cycles }
    }

    pub fn str(&self) -> &KhAlgStr<R> { 
        &self.str
    }

    pub fn deg_shift(&self) -> (isize, isize) { 
        self.deg_shift
    }

    pub fn is_reduced(&self) -> bool { 
        self.reduced
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        range_of(self.support())
    }

    pub fn q_range(&self) -> RangeInclusive<isize> {
        range_of(self.support().flat_map(|i| 
            self[i].raw_gens().iter().map(|x| x.q_deg())
        ))
    }

    pub fn canon_cycles(&self) -> &Vec<KhChain<R>> { 
        &self.canon_cycles
    }

    pub fn inner(&self) -> &ChainComplex<KhGen, R> {
        &self.inner
    }

    pub fn truncated(&self, range: RangeInclusive<isize>) -> Self {
        Self::new_impl(
            self.inner.truncated(range), 
            self.str.clone(), 
            self.deg_shift, 
            self.reduced, 
            self.canon_cycles.clone()
        )
    }

    pub fn gen_grid(&self) -> Grid2<Summand<KhGen, R>> { 
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

    pub fn into_bigraded(self) -> KhComplexBigraded<R> {
        assert!(self.str.h().is_zero());
        assert!(self.str.t().is_zero());

        let str = self.str.clone();
        let deg_shift = self.deg_shift;
        let reduced = self.reduced;
        let canon_cycles = self.canon_cycles.clone();
        let summands = self.gen_grid();

        let inner = ChainComplex2::new(summands, isize2(1, 0), move |idx, x| { 
            let i = idx.0;
            let x = KhChain::from(x.clone());
            let dx = self.d(i, &x);
            dx.into_iter().collect()
        });

        KhComplexBigraded::new_impl(inner, str, deg_shift, reduced, canon_cycles)
    }

    pub fn deg_shift_for(l: &Link, reduced: bool) -> (isize, isize) {
        let (n_pos, n_neg) = l.signed_crossing_nums();
        let (n_pos, n_neg) = (n_pos as isize, n_neg as isize);
        let h = -n_neg;
        let q = n_pos - 2 * n_neg;
        let e = if reduced { 1 } else { 0 };
        (h, q + e)
    }

    delegate! { 
        to self.inner { 
            pub fn d(&self, i: isize, z: &KhChain<R>) -> KhChain<R>;
        }
    }
}

impl<R> Index<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = KhComplexSummand<R>;

    delegate! { 
        to self.inner {
            fn index(&self, index: isize) -> &Self::Output;
        }
    }
}

impl<R> GridTrait<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Support = std::vec::IntoIter<isize>;
    type Item = KhComplexSummand<R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Support;
            fn is_supported(&self, i: isize) -> bool;
            fn get(&self, i: isize) -> &Self::Item;
            fn get_default(&self) -> &Self::Item;
        }
    }
}

impl<R> ChainComplexTrait<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;
    type Element = KhChain<R>;

    delegate! { 
        to self.inner { 
            fn rank(&self, i: isize) -> usize;
            fn d_deg(&self) -> isize;
            fn d(&self, i: isize, z: &Self::Element) -> Self::Element;
            fn d_matrix(&self, i: isize) -> SpMat<R>;
        }
    }
}

impl<R> KhComplex<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology(&self) -> KhHomology<R> {
        self.into()
    }
}

#[derive(Clone)]
pub struct KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    inner: ChainComplex2<KhGen, R>,
    str: KhAlgStr<R>,
    deg_shift: (isize, isize),
    reduced: bool,
    canon_cycles: Vec<KhChain<R>>,
}

impl<R> KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(l: &Link, h: &R, t: &R, reduced: bool) -> Self { 
        assert!(h.is_zero() && t.is_zero());

        let c = KhComplex::new(&l, h, t, reduced);
        c.into_bigraded()
    }

    fn new_impl(
        inner: ChainComplex2<KhGen, R>,
        str: KhAlgStr<R>,
        deg_shift: (isize, isize),
        reduced: bool,
        canon_cycles: Vec<KhChain<R>>
    ) -> Self {
        Self { inner, str, deg_shift, reduced, canon_cycles }
    }

    pub fn str(&self) -> &KhAlgStr<R> { 
        &self.str
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        range_of(self.support().map(|idx| idx.0))
    }

    pub fn q_range(&self) -> RangeInclusive<isize> {
        range_of(self.support().map(|idx| idx.1))
    }

    pub fn deg_shift(&self) -> (isize, isize) { 
        self.deg_shift
    }

    pub fn is_reduced(&self) -> bool { 
        self.reduced
    }

    pub fn canon_cycles(&self) -> &Vec<KhChain<R>> { 
        &self.canon_cycles
    }

    pub fn inner(&self) -> &ChainComplex2<KhGen, R> {
        &self.inner
    }
}

impl<R> Index<(isize, isize)> for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = KhComplexSummand<R>;

    delegate! { 
        to self.inner { 
            fn index(&self, index: (isize, isize)) -> &Self::Output;
        }
    }
}

impl<R> GridTrait<isize2> for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Support = std::vec::IntoIter<isize2>;
    type Item = KhComplexSummand<R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Support;
            fn is_supported(&self, i: isize2) -> bool;
            fn get(&self, i: isize2) -> &Self::Item;
            fn get_default(&self) -> &Self::Item;
        }
    }
}

impl<R> ChainComplexTrait<isize2> for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;
    type Element = KhChain<R>;

    delegate! { 
        to self.inner { 
            fn rank(&self, i: isize2) -> usize;
            fn d_deg(&self) -> isize2;
            fn d(&self, i: isize2, z: &Self::Element) -> Self::Element;
            fn d_matrix(&self, i: isize2) -> SpMat<Self::R>;
        }
    }
}

impl<R> KhComplexBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology(&self) -> KhHomologyBigraded<R> {
        self.into()
    }
}

#[cfg(test)]
mod tests {
    use yui_homology::{ChainComplexTrait, SummandTrait};
    use yui_link::Link;

    use crate::kh::KhComplexBigraded;

    use super::KhComplex;

    #[test]
    fn ckh_trefoil() {
        let l = Link::trefoil();
        let c = KhComplex::new(&l, &0, &0, false);

        assert_eq!(c.h_range(), -3..=0);

        cfg_if::cfg_if! { 
        if #[cfg(feature = "old")] { 
            assert_eq!(c[-3].rank(), 8);
            assert_eq!(c[-2].rank(), 12);
            assert_eq!(c[-1].rank(), 6);
            assert_eq!(c[ 0].rank(), 4);    
        } else { 
            assert_eq!(c[-3].rank(), 2);
            assert_eq!(c[-2].rank(), 2);
            assert_eq!(c[-1].rank(), 0);
            assert_eq!(c[ 0].rank(), 2);
        }}  

        c.check_d_all();
    }

    #[test]
    fn ckh_trefoil_red() {
        let l = Link::trefoil();
        let c = KhComplex::new(&l, &0, &0, true);

        assert_eq!(c.h_range(), -3..=0);

        cfg_if::cfg_if! { 
        if #[cfg(feature = "old")] { 
            assert_eq!(c[-3].rank(), 4);
            assert_eq!(c[-2].rank(), 6);
            assert_eq!(c[-1].rank(), 3);
            assert_eq!(c[ 0].rank(), 2);
        } else {
            assert_eq!(c[-3].rank(), 1);
            assert_eq!(c[-2].rank(), 1);
            assert_eq!(c[-1].rank(), 0);
            assert_eq!(c[ 0].rank(), 1);
        }}

        c.check_d_all();
    }

    #[test]
    fn ckh_bigr_trefoil() {
        let l = Link::trefoil();
        let c = KhComplexBigraded::new(&l, &0, &0, false);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c.q_range(), -9..=-1);

        cfg_if::cfg_if! { 
        if #[cfg(feature = "old")] { 
            // TODO
        } else { 
            assert_eq!(c[(-3, -9)].rank(), 1);
            assert_eq!(c[(-3, -7)].rank(), 1);
            assert_eq!(c[(-2, -7)].rank(), 1);
            assert_eq!(c[(-2, -5)].rank(), 1);
            assert_eq!(c[(0, -3)].rank(), 1);
            assert_eq!(c[(0, -1)].rank(), 1);
        }}

        c.check_d_all();
    }

    #[test]
    fn ckh_bigr_red_trefoil() {
        let l = Link::trefoil();
        let c = KhComplexBigraded::new(&l, &0, &0, true);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c.q_range(), -8..=-2);

        cfg_if::cfg_if! { 
        if #[cfg(feature = "old")] { 
            // TODO
        } else { 
            assert_eq!(c[(-3, -8)].rank(), 1);
            assert_eq!(c[(-2, -6)].rank(), 1);
            assert_eq!(c[(0, -2)].rank(), 1);
        }}

        c.check_d_all();
    }
}
