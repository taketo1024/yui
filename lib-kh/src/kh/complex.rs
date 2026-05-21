use std::ops::{RangeInclusive, Index};
use std::sync::OnceLock;

use delegate::delegate;
use yui_core::{IteratorExt, Ring, RingOps, EucRing, EucRingOps};
use yui_link::Link;
use yui_homology::{ChainComplex, ToSeqString, ToTableString, Grid2, GridIter, GridTrait, Summand};

use crate::kh::chain::KhChain;
use crate::kh::internal::v1::cube::KhCube;
use crate::kh::{KhState, KhHomology};
use crate::misc::make_gen_grid;

use super::KhAlg;

pub type KhComplexSummand<R> = Summand<KhState, R>;

// TODO: Make KhComplexTrait, and split impl into KhComplexV1 and V2. 

#[derive(Clone)]
pub struct KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    inner: ChainComplex<KhState, R>,
    str: KhAlg<R>,
    cube: KhCube<R>,
    deg_shift: (isize, isize),
    reduced: bool,
    canon_cycles: Vec<KhChain<R>>,
    gen_grid: OnceLock<Grid2<KhComplexSummand<R>>>,
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
        let complex = cube.clone().into_complex();

        let canon_cycles = if t.is_zero() && l.is_knot() {
            let p = l.min_edge().unwrap();
            Self::make_canon_cycles(l, p, &R::zero(), h, reduced, deg_shift)
        } else { 
            vec![]
        };

        KhComplex::new_impl(complex, str, cube, deg_shift, reduced, canon_cycles)
    }

    pub(crate) fn new_impl(inner: ChainComplex<KhState, R>, str: KhAlg<R>, cube: KhCube<R>, deg_shift: (isize, isize), reduced: bool, canon_cycles: Vec<KhChain<R>>) -> Self {
        KhComplex { inner, str, cube, deg_shift, reduced, canon_cycles, gen_grid: OnceLock::new() }
    }

    pub fn str(&self) -> &KhAlg<R> { 
        &self.str
    }

    pub fn cube(&self) -> &KhCube<R> {
        &self.cube
    }

    pub fn deg_shift(&self) -> (isize, isize) { 
        self.deg_shift
    }

    pub fn is_reduced(&self) -> bool { 
        self.reduced
    }

    pub fn h_range(&self) -> RangeInclusive<isize> {
        self.support().copied().range().unwrap_or(0..=-1)
    }

    pub fn q_range(&self) -> RangeInclusive<isize> {
        self.support().flat_map(|&i|
            self[i].raw_generators().iter().map(|x| x.q_deg())
        ).range().unwrap_or(0..=-1)
    }

    pub fn canon_cycles(&self) -> &Vec<KhChain<R>> { 
        &self.canon_cycles
    }

    pub fn inner(&self) -> &ChainComplex<KhState, R> {
        &self.inner
    }

    fn gen_grid(&self) -> &Grid2<KhComplexSummand<R>> {
        self.gen_grid.get_or_init(|| make_gen_grid(self.inner.summands()))
    }

    pub fn deg_shift_for(l: &Link, reduced: bool) -> (isize, isize) {
        let (n_pos, n_neg) = l.n_signed_crossings();
        let (n_pos, n_neg) = (n_pos as isize, n_neg as isize);
        let h = -n_neg;
        let q = n_pos - 2 * n_neg;
        let e = if reduced { 1 } else { 0 };
        (h, q + e)
    }

    delegate! {
        to self.inner {
            pub fn d_deg(&self) -> isize;
            pub fn d(&self, i: isize, z: &KhChain<R>) -> KhChain<R>;
            pub fn describe_d(&self) -> String;
            pub fn describe_d_at(&self, i: isize) -> String;
        }
    }
}

impl<R> KhComplex<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology(&self) -> KhHomology<R> {
        self.into()
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

impl<R> Index<(isize, isize)> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = KhComplexSummand<R>;

    delegate! { 
        to self.gen_grid() {
            fn index(&self, index: (isize, isize)) -> &Self::Output;
        }
    }
}

impl<R> GridTrait<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Item = KhComplexSummand<R>;
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

impl<R> ToSeqString<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    delegate! {
        to self.inner { 
            fn label(&self) -> String;
            fn indices(&self) -> Vec<isize>;
            fn entry_at(&self, i: &isize) -> String;
        }
    }
}

impl<R> ToTableString<isize> for KhComplex<R>
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
        use yui_link::Link;

    use super::KhComplex;

    #[test]
    fn ckh_trefoil() {
        let l = Link::test_data("3_1");
        let c = KhComplex::new(&l, &0, &0, false);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c.q_range(), -9..=-1);

        assert_eq!(c[-3].rank(), 2);
        assert_eq!(c[-2].rank(), 2);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);

        c.inner().check_d_all();
    }

    #[test]
    fn ckh_trefoil_red() {
        let l = Link::test_data("3_1");
        let c = KhComplex::new(&l, &0, &0, true);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c.q_range(), -8..=-2);

        assert_eq!(c[-3].rank(), 1);
        assert_eq!(c[-2].rank(), 1);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 1);

        c.inner().check_d_all();
    }

    #[test]
    fn ckh_trefoil_bigr() {
        let l = Link::test_data("3_1");
        let c = KhComplex::new(&l, &0, &0, false);

        assert_eq!(c[(-3, -9)].rank(), 1);
        assert_eq!(c[(-3, -7)].rank(), 1);
        assert_eq!(c[(-2, -7)].rank(), 1);
        assert_eq!(c[(-2, -5)].rank(), 1);
        assert_eq!(c[(0, -3)].rank(), 1);
        assert_eq!(c[(0, -1)].rank(), 1);
    }

    #[test]
    fn ckh_trefoil_bigr_red() {
        let l = Link::test_data("3_1");
        let c = KhComplex::new(&l, &0, &0, true);

        assert_eq!(c[(-3, -8)].rank(), 1);
        assert_eq!(c[(-2, -6)].rank(), 1);
        assert_eq!(c[(0, -2)].rank(), 1);
    }
}

#[cfg(test)]
mod tests_v1 {
        use yui_link::Link;

    use super::KhComplex;

    #[test]
    fn ckh_trefoil() {
        let l = Link::test_data("3_1");
        let c = KhComplex::new_no_simplify(&l, &0, &0, false);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c[-3].rank(), 8);
        assert_eq!(c[-2].rank(), 12);
        assert_eq!(c[-1].rank(), 6);
        assert_eq!(c[ 0].rank(), 4);    

        c.inner().check_d_all();
    }

    #[test]
    fn ckh_trefoil_red() {
        let l = Link::test_data("3_1");
        let c = KhComplex::new_no_simplify(&l, &0, &0, true);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c[-3].rank(), 4);
        assert_eq!(c[-2].rank(), 6);
        assert_eq!(c[-1].rank(), 3);
        assert_eq!(c[ 0].rank(), 2);

        c.inner().check_d_all();
    }
}
