use std::collections::HashMap;
use std::ops::Index;
use yui_matrix::sparse::*;
use yui_core::{Ring, RingOps};

use crate::{RModStr, RModGrid, GenericRModGrid, GenericRModStr, ChainComplex};
use crate::utils::reducer2::ChainReducer;

pub struct Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>,
    C::Output: RModStr<R = C::R>
{ 
    original: C,
    grid: GenericRModGrid<GenericRModStr<C::R>, C::IdxIter>,
    d_matrices: HashMap<C::Idx, SpMat<C::R>>,
}

impl<C> From<C> for Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>,
    C::Output: RModStr<R = C::R>
{
    fn from(c: C) -> Self {
        let mut red = ChainReducer::new(&c);
        red.process();

        let d_matrices = c.indices().map(|i| 
            (i, red.take_matrix(i))
        ).collect::<HashMap<_, _>>();

        let grid = GenericRModGrid::new(c.indices(), |i| { 
            let n = d_matrices[&i].cols();
            if n > 0 { 
                let s = GenericRModStr::new(n, vec![]);
                Some(s)
            } else { 
                None
            }
        });

        Self { original: c, grid, d_matrices }
    }
}

impl<C> Index<C::Idx> for Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>,
    C::Output: RModStr<R = C::R>
{
    type Output = GenericRModStr<C::R>;

    fn index(&self, index: C::Idx) -> &Self::Output {
        &self.grid[index]
    }
}

impl<C> RModGrid for Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>,
    C::Output: RModStr<R = C::R>
{
    type R = C::R;
    type Idx = C::Idx;
    type IdxIter = C::IdxIter;

    fn contains_idx(&self, k: Self::Idx) -> bool {
        self.original.contains_idx(k)
    }

    fn indices(&self) -> Self::IdxIter {
        self.original.indices()
    }
}

impl<C> ChainComplex for Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>,
    C::Output: RModStr<R = C::R>
{
    fn d_degree(&self) -> Self::Idx {
        self.original.d_degree()
    }

    fn d_matrix(&self, k: Self::Idx) -> SpMat<Self::R> {
        if self.contains_idx(k) { 
            self.d_matrices[&k].clone()
        } else {
            let m = self[k + self.d_degree()].rank();
            let n = self[k].rank();
            SpMat::zero((m, n))
        }
    }
}

#[cfg(test)]
mod tests { 
    use super::*;
    use crate::{RModStr, GenericHomology};
    use crate::test::TestChainComplex;

    #[test]
    fn cancel_pair() { 
        let c = TestChainComplex::<i32>::descending(
            vec![ SpMat::from_vec((1, 1), vec![1]) ],
        );
        let c = Reduced::from(c);

        assert_eq!(c[0].rank(), 0);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn torsion() { 
        let c = TestChainComplex::<i32>::descending( 
            vec![ SpMat::from_vec((1, 1), vec![2]) ],
        );
        let c = Reduced::from(c);

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 1);
    }

    #[test]
    fn d3() {
        let c = TestChainComplex::<i32>::d3();
        let c = Reduced::from(c);

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 0);
        assert_eq!(c[2].rank(), 0);
        assert_eq!(c[3].rank(), 0);

        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[3].rank(), 0);
    }

    #[test]
    fn s2() {
        let c = TestChainComplex::<i32>::s2();
        let c = Reduced::from(c);

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 0);
        assert_eq!(c[2].rank(), 1);

        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 1);
    }

    #[test]
    fn t2() {
        let c = TestChainComplex::<i32>::t2();
        let c = Reduced::from(c);

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 2);
        assert_eq!(c[2].rank(), 1);

        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[2].rank(), 1);
    }

    #[test]
    fn rp2() {
        let c = TestChainComplex::<i32>::rp2();
        let c = Reduced::from(c);

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 1);
        assert_eq!(c[2].rank(), 1);

        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].tors(), &vec![2]);
        assert_eq!(h[2].rank(), 0);
    }

    #[test]
    fn d3_dual() {
        let c = TestChainComplex::<i32>::d3().dual();
        let c = Reduced::from(c);

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 0);
        assert_eq!(c[2].rank(), 0);
        assert_eq!(c[3].rank(), 0);

        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[3].rank(), 0);
    }

    #[test]
    fn s2_dual() {
        let c = TestChainComplex::<i32>::s2().dual();
        let c = Reduced::from(c);

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 0);
        assert_eq!(c[2].rank(), 1);

        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 1);
    }

    #[test]
    fn t2_dual() {
        let c = TestChainComplex::<i32>::t2().dual();
        let c = Reduced::from(c);

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 2);
        assert_eq!(c[2].rank(), 1);

        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[2].rank(), 1);
    }

    #[test]
    fn rp2_dual() {
        let c = TestChainComplex::<i32>::rp2().dual();
        let c = Reduced::from(c);

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 1);
        assert_eq!(c[2].rank(), 1);

        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[2].tors(), &vec![2]);
    }
    
}