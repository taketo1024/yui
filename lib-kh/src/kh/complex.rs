use std::ops::{RangeInclusive, Index};
use cartesian::cartesian;

use delegate::delegate;
use yui_core::{Ring, RingOps, EucRing, EucRingOps};
use yui_link::Link;
use yui_homology::{isize2, ChainComplexTrait, Grid2, GridTrait, ChainComplex, Summand};
use yui_matrix::sparse::SpMat;

use crate::kh::r#gen::KhChain;
use crate::kh::{KhChainGen, KhHomology};
use crate::misc::range_of;

use super::KhAlg;

pub type KhComplexSummand<R> = Summand<KhChainGen, R>;

#[derive(Clone)]
pub struct KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    inner: ChainComplex<KhChainGen, R>,
    str: KhAlg<R>,
    deg_shift: (isize, isize),
    reduced: bool,
    canon_cycles: Vec<KhChain<R>>,
}

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        use crate::kh::internal::v2::builder::TngComplexBuilder;

        assert!(!reduced || (!l.is_empty() && t.is_zero()));

        TngComplexBuilder::build_kh_complex(l, h, t, reduced)
    }

    pub fn new_no_simplify(l: &Link, h: &R, t: &R, reduced: bool) -> Self { 
        use crate::kh::internal::v1::cube::KhCube;

        assert!(!reduced || (!l.is_empty() && t.is_zero()));

        let red_e = reduced.then(|| l.min_edge().unwrap());
        let deg_shift = Self::deg_shift_for(l, reduced);
        
        let cube = KhCube::new(l, h, t, red_e, deg_shift);
        let str = cube.str().clone();
        let complex = cube.into_complex();

        let canon_cycles = if t.is_zero() && l.is_knot() {
            let p = l.min_edge().unwrap();
            Self::make_canon_cycles(l, p, &R::zero(), h, reduced, deg_shift)
        } else { 
            vec![]
        };

        KhComplex::new_impl(complex, str, deg_shift, reduced, canon_cycles)
    }

    pub(crate) fn new_impl(inner: ChainComplex<KhChainGen, R>, str: KhAlg<R>, deg_shift: (isize, isize), reduced: bool, canon_cycles: Vec<KhChain<R>>) -> Self { 
        KhComplex { inner, str, deg_shift, reduced, canon_cycles }
    }

    pub fn str(&self) -> &KhAlg<R> { 
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

    pub fn inner(&self) -> &ChainComplex<KhChainGen, R> {
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

    pub fn gen_grid(&self) -> Grid2<Summand<KhChainGen, R>> { 
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

    pub fn deg_shift_for(l: &Link, reduced: bool) -> (isize, isize) {
        let (n_pos, n_neg) = l.count_signed_crossings();
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

#[cfg(test)]
mod tests {
    use yui_homology::{ChainComplexTrait, SummandTrait};
    use yui_link::Link;

    use super::KhComplex;

    #[test]
    fn ckh_trefoil() {
        let l = Link::trefoil();
        let c = KhComplex::new(&l, &0, &0, false);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c.q_range(), -9..=-1);

        assert_eq!(c[-3].rank(), 2);
        assert_eq!(c[-2].rank(), 2);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);

        c.check_d_all();
    }

    #[test]
    fn ckh_trefoil_red() {
        let l = Link::trefoil();
        let c = KhComplex::new(&l, &0, &0, true);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c.q_range(), -8..=-2);

        assert_eq!(c[-3].rank(), 1);
        assert_eq!(c[-2].rank(), 1);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 1);

        c.check_d_all();
    }

    #[test]
    fn gen_grid() {
        let l = Link::trefoil();
        let c = KhComplex::new(&l, &0, &0, false).gen_grid();

        assert_eq!(c[(-3, -9)].rank(), 1);
        assert_eq!(c[(-3, -7)].rank(), 1);
        assert_eq!(c[(-2, -7)].rank(), 1);
        assert_eq!(c[(-2, -5)].rank(), 1);
        assert_eq!(c[(0, -3)].rank(), 1);
        assert_eq!(c[(0, -1)].rank(), 1);
    }

    #[test]
    fn gen_grid_red() {
        let l = Link::trefoil();
        let c = KhComplex::new(&l, &0, &0, true).gen_grid();

        assert_eq!(c[(-3, -8)].rank(), 1);
        assert_eq!(c[(-2, -6)].rank(), 1);
        assert_eq!(c[(0, -2)].rank(), 1);
    }
}

#[cfg(test)]
mod tests_v1 {
    use yui_homology::{ChainComplexTrait, SummandTrait};
    use yui_link::Link;

    use super::KhComplex;

    #[test]
    fn ckh_trefoil() {
        let l = Link::trefoil();
        let c = KhComplex::new_no_simplify(&l, &0, &0, false);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c[-3].rank(), 8);
        assert_eq!(c[-2].rank(), 12);
        assert_eq!(c[-1].rank(), 6);
        assert_eq!(c[ 0].rank(), 4);    

        c.check_d_all();
    }

    #[test]
    fn ckh_trefoil_red() {
        let l = Link::trefoil();
        let c = KhComplex::new_no_simplify(&l, &0, &0, true);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c[-3].rank(), 4);
        assert_eq!(c[-2].rank(), 6);
        assert_eq!(c[-1].rank(), 3);
        assert_eq!(c[ 0].rank(), 2);

        c.check_d_all();
    }
}
