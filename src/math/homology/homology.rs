use std::fmt::Display;
use std::ops::Index;

use crate::math::traits::{Ring, RingOps, EucRing, EucRingOps};
use super::base::{GradedRModStr, RModStr, GenericRModStr, RModGrid, AdditiveIndex, AdditiveIndexRange};
use super::complex::ChainComplex;
use super::utils::homology_calc::HomologyCalc;

pub trait HomologyComputable<S>: ChainComplex
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>,
    S: RModStr<R = Self::R>
{
    fn homology_at(&self, i: Self::Index) -> S;
}

impl<R, C> HomologyComputable<GenericRModStr<R>> for C
where 
    R: EucRing, for<'x> &'x R: EucRingOps<R>,
    C: ChainComplex<R = R>,
    C::Output: RModStr<R = R>
{ 
    fn homology_at(&self, k: C::Index) -> GenericRModStr<Self::R> {
        let d1 = self.d_matrix(k - self.d_degree());
        let d2 = self.d_matrix(k);
        HomologyCalc::calculate(d1, d2)
    }
}

pub trait Homology: GradedRModStr
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>
{
    fn is_zero(&self) -> bool {
        self.range().all(|i| self[i].is_zero())
    }

    fn is_free(&self) -> bool {
        self.range().all(|i| self[i].is_free())
    }

    fn fmt_default(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in self.range() { 
            write!(f, "H[{}]: {}\n", i, self[i])?
        }
        Ok(())
    }
}

pub struct GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    grid: RModGrid<GenericRModStr<R>, I>
}

impl<R, C> From<C> for GenericHomology<R, C::IndexRange>
where
    R: Ring, for<'x> &'x R: RingOps<R>,
    C: HomologyComputable<GenericRModStr<R>, R = R>,
    C::IndexRange: Clone,
    C::Output: RModStr<R = C::R>
{
    fn from(c: C) -> Self {
        let range = c.range();
        let grid = RModGrid::new(range, |i| {
            let h_i = c.homology_at(i);
            if !h_i.is_zero() { Some(h_i) } else { None }
        });
        Self { grid }
    }
}

impl<R, I> GradedRModStr for GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    type R = R;
    type Index = I::Item;
    type IndexRange = I;

    fn in_range(&self, k: Self::Index) -> bool {
        self.grid.in_range(k)
    }

    fn range(&self) -> Self::IndexRange {
        self.grid.range()
    }
}

impl<R, I> Index<I::Item> for GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    type Output = GenericRModStr<R>;

    fn index(&self, k: I::Item) -> &Self::Output {
        &self.grid[k]
    }
}

impl<R, I> Homology for GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
}

impl<R, I> Display for GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt_default(f)
    }
}

#[cfg(test)]
mod tests { 
    use crate::math::matrix::sparse::*;
    use super::*;
    use super::super::complex::tests::*;

    #[test]
    fn cancel_pair() { 
        let c = TestChainComplex::<i32>::descending(
            vec![ SpMat::from_vec((1, 1), vec![1]) ],
        );

        let h = GenericHomology::from(c);
        
        assert_eq!(h[0].rank(), 0);
        assert!(h[0].is_free());
    }

    #[test]
    fn torsion() { 
        let c = TestChainComplex::<i32>::descending( 
            vec![ SpMat::from_vec((1, 1), vec![2]) ],
        );

        let h = GenericHomology::from(c);
        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[0].tors(), &vec![2]);
    }

    #[test]
    fn d3() {
        let c = TestChainComplex::<i32>::d3();
        let h = GenericHomology::from(c);

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
        let c = TestChainComplex::<i32>::s2();
        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn t2() {
        let c = TestChainComplex::<i32>::t2();
        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn rp2() {
        let c = TestChainComplex::<i32>::rp2();
        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].tors(), &vec![2]);
        assert_eq!(h[1].is_free(), false);

        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn d3_dual() {
        let c = TestChainComplex::<i32>::d3().dual();
        let h = GenericHomology::from(c);

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
    fn s2_dual() {
        let c = TestChainComplex::<i32>::s2().dual();
        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn t2_dual() {
        let c = TestChainComplex::<i32>::t2().dual();
        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn rp2_dual() {
        let c = TestChainComplex::<i32>::rp2().dual();
        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[2].tors(), &vec![2]);
        assert_eq!(h[2].is_free(), false);
    }
}