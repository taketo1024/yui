use std::ops::Index;

use delegate::delegate;
use yui_core::{EucRing, EucRingOps, Ring, RingOps, Deg, isize2, isize3};
use yui_matrix::sparse::{SpVec, Trans};

use crate::utils::HomologyCalc;
use crate::{GridBase, GridIter, DisplayForGrid, ChainComplexTrait, RModStr};

use super::grid::GridTrait;

pub type Homology<R>  = HomologyBase<isize,  R>;
pub type Homology2<R> = HomologyBase<isize2, R>;
pub type Homology3<R> = HomologyBase<isize3, R>;

#[derive(Debug, Clone)]
pub struct HomologySummand<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    rank: usize,
    tors: Vec<R>,
    trans: Option<Trans<R>>
}

impl<R> HomologySummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(rank: usize, tors: Vec<R>, trans: Option<Trans<R>>) -> Self { 
        Self { rank, tors, trans }
    }

    pub fn zero() -> Self { 
        Self::new(0, vec![], None)
    }

    pub fn ngens(&self) -> usize { 
        self.rank + self.tors.len()
    }

    pub fn gen_vec(&self, i: usize) -> Option<SpVec<R>> {
        assert!(i < self.ngens());
        self.trans.as_ref().map(|t| t.backward_mat().col_vec(i))
    }

    pub fn vec_from_cpx(&self, v: &SpVec<R>) -> Option<SpVec<R>> { 
        self.trans.as_ref().map(|t| t.forward(v))
    }

    pub fn vec_to_cpx(&self, v: &SpVec<R>) -> Option<SpVec<R>> { 
        self.trans.as_ref().map(|t| t.backward(v))
    }
}

impl<R> Default for HomologySummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn default() -> Self {
        Self::zero()
    }
}

impl<R> RModStr for HomologySummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn rank(&self) -> usize {
        self.rank
    }

    fn tors(&self) -> &Vec<Self::R> {
        &self.tors
    }
}

impl<R> DisplayForGrid for HomologySummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn display_for_grid(&self) -> String {
        self.math_symbol()
    }
}

pub struct HomologyBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    summands: GridBase<I, HomologySummand<R>>
}

impl<I, R> GridTrait<I> for HomologyBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type Itr = GridIter<I>;
    type E = HomologySummand<R>;
    
    delegate! { 
        to self.summands { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: I) -> bool;
            fn get(&self, i: I) -> &Self::E;
        }
    }
}

impl<I, R> HomologyBase<I, R>
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn compute_from<C>(complex: &C, with_trans: bool) -> Self
    where 
        C: ChainComplexTrait<I, R = R>,
        C::E: RModStr<R = R>
    {
        let summands = GridBase::new(
            complex.support(), 
            |i| {
                let i0 = i - complex.d_deg();
                let d0 = complex.d_matrix(i0);
                let d1 = complex.d_matrix(i );
                HomologyCalc::calculate(d0, d1, with_trans)
            }
        );
        Self { summands }
    }
}

impl<I, R> Index<I> for HomologyBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = HomologySummand<R>;
    fn index(&self, i: I) -> &Self::Output {
        self.get(i)
    }
}

impl<R> Index<(isize, isize)> for Homology2<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = HomologySummand<R>;
    fn index(&self, i: (isize, isize)) -> &Self::Output {
        self.get(i.into())
    }
}

impl<R> Index<(isize, isize, isize)> for Homology3<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = HomologySummand<R>;
    fn index(&self, i: (isize, isize, isize)) -> &Self::Output {
        self.get(i.into())
    }
}

#[cfg(test)]
mod tests { 
    use yui_matrix::sparse::SpVec;
    use crate::{ChainComplex, RModStr};

    #[test]
    fn zero() { 
        let c = ChainComplex::<i32>::zero();
        let h = c.homology(false);
        
        assert!(h[0].is_zero());
    }

    #[test]
    fn single() { 
        let c = ChainComplex::<i32>::one();
        let h = c.homology(false);
        
        assert_eq!(h[0].rank(), 1);
        assert!( h[0].is_free());
        assert!(!h[0].is_zero());
    }

    #[test]
    fn one_to_one() { 
        let c = ChainComplex::<i32>::one_one(1);
        let h = c.homology(false);

        assert!(h[0].is_zero());
        assert!(h[1].is_zero());
    }

    #[test]
    fn two_to_one() { 
        let c = ChainComplex::<i32>::two_one(1, -1);
        let h = c.homology(false);

        assert!(h[0].is_zero());
        assert_eq!(h[1].rank(), 1);
        assert!(h[1].is_free());
    }

    #[test]
    fn one_to_two() { 
        let c = ChainComplex::<i32>::one_two(1, -1);
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());
        assert!(h[1].is_zero());
    }

    #[test]
    fn torsion() { 
        let c = ChainComplex::<i32>::one_one(2);
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[0].tors(), &vec![2]);
        assert!(!h[0].is_free());
    }

    #[test]
    fn d3() {
        let c = ChainComplex::<i32>::d3();
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[2].is_free(), true);

        assert_eq!(h[3].rank(), 0);
        assert_eq!(h[3].is_free(), true);
    }

    #[test]
    fn s2() {
        let c = ChainComplex::<i32>::s2();
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn t2() {
        let c = ChainComplex::<i32>::t2();
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn rp2() {
        let c = ChainComplex::<i32>::rp2();
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].tors(), &vec![2]);
        assert_eq!(h[1].is_free(), false);

        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn s2_gens() {
        let c = ChainComplex::<i32>::s2();
        let h = c.homology(true);

        let h2 = &h[2];
        assert_eq!(h2.ngens(), 1);

        let z = h2.gen_vec(0).unwrap();
        assert!(!z.is_zero());
        assert!(c.d(2, &z).is_zero());
        assert_eq!(h2.vec_from_cpx(&z), Some(SpVec::from(vec![1])));
    }

    #[test]
    fn t2_gens() {
        let c = ChainComplex::<i32>::t2();
        let h = c.homology(true);

        let h2 = &h[2];
        assert_eq!(h2.ngens(), 1);

        let z = h2.gen_vec(0).unwrap();
        assert!(!z.is_zero());
        assert!(c.d(2, &z).is_zero());
        assert_eq!(h2.vec_from_cpx(&z), Some(SpVec::from(vec![1])));
        assert_eq!(h2.ngens(), 1);

        let h1 = &h[1];
        assert_eq!(h1.ngens(), 2);

        let a = h1.gen_vec(0).unwrap();
        let b = h1.gen_vec(1).unwrap();

        assert!(!a.is_zero());
        assert!(!b.is_zero());
        assert!(c.d(1, &a).is_zero());
        assert!(c.d(1, &b).is_zero());
        assert_eq!(h1.vec_from_cpx(&a), Some(SpVec::from(vec![1, 0])));
        assert_eq!(h1.vec_from_cpx(&b), Some(SpVec::from(vec![0, 1])));
    }

    #[test]
    fn rp2_gens() {
        let c = ChainComplex::<i32>::rp2();
        let h = c.homology(true);

        let h1 = &h[1];
        assert_eq!(h1.ngens(), 1);

        let z = h1.gen_vec(0).unwrap();

        assert!(!z.is_zero());
        assert!(c.d(1, &z).is_zero());
        assert_eq!(h1.vec_from_cpx(&z), Some(SpVec::from(vec![1])));
    }
}