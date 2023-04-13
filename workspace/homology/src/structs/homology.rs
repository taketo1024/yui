use std::fmt::Display;
use std::ops::Index;

use yui_core::{Ring, RingOps};
use crate::{GridIdx, GridItr, RModStr, RModGrid, GenericRModStr, GenericRModGrid, Homology, HomologyComputable, Grid};

pub struct GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr,
    I::Item: GridIdx
{
    grid: GenericRModGrid<GenericRModStr<R>, I>
}

impl<R, C> From<C> for GenericHomology<R, C::IdxIter>
where
    R: Ring, for<'x> &'x R: RingOps<R>,
    C: HomologyComputable<GenericRModStr<R>, R = R>,
    C::IdxIter: Clone,
    C::Output: RModStr<R = C::R>
{
    fn from(c: C) -> Self {
        let range = c.indices();
        let grid = GenericRModGrid::new(range, |i| {
            let h_i = c.homology_at(i);
            if !h_i.is_zero() { Some(h_i) } else { None }
        });
        Self { grid }
    }
}

impl<R, I> Grid for GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr,
    I::Item: GridIdx
{
    type Idx = I::Item;
    type IdxIter = I;

    fn contains_idx(&self, k: Self::Idx) -> bool {
        self.grid.contains_idx(k)
    }

    fn indices(&self) -> Self::IdxIter {
        self.grid.indices()
    }
}

impl<R, I> Index<I::Item> for GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr,
    I::Item: GridIdx
{
    type Output = GenericRModStr<R>;

    fn index(&self, k: I::Item) -> &Self::Output {
        &self.grid[k]
    }
}

impl<R, I> Homology for GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr,
    I::Item: GridIdx
{
}

impl<R, I> Display for GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr,
    I::Item: GridIdx
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt_default(f, "H")
    }
}

#[cfg(test)]
mod tests { 
    use super::*;
    use yui_matrix::sparse::*;
    use crate::test::TestChainComplex;

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