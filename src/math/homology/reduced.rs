use std::collections::HashMap;
use std::ops::Index;
use crate::math::matrix::sparse::*;
use crate::math::traits::{Ring, RingOps};
use super::base::{GradedRModStr, RModStr, RModGrid, GenericRModStr};
use super::complex::ChainComplex;
use super::utils::reducer::ChainReducer;

pub struct Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>,
    C::Output: RModStr<R = C::R>
{ 
    original: C,
    grid: RModGrid<GenericRModStr<C::R>, C::IndexRange>,
    d_matrices: HashMap<C::Index, SpMat<C::R>>,
}

impl<C> From<C> for Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>,
    C::Output: RModStr<R = C::R>
{
    fn from(c: C) -> Self {
        let mut d_matrices = HashMap::new();

        for k in c.range() {
            let deg = c.d_degree();
            let (k0, k1, k2) = (k - deg, k, k + deg);

            let a0 = d_matrices.remove(&k0).unwrap_or(c.d_matrix(k0)); // prev
            let a1 = d_matrices.remove(&k1).unwrap_or(c.d_matrix(k1)); // target
            let a2 = c.d_matrix(k2); // next

            let (b0, b1, b2) = ChainReducer::reduce(a0, a1, a2);
            
            d_matrices.insert(k0, b0);
            d_matrices.insert(k1, b1);
            d_matrices.insert(k2, b2);
        }

        let grid = RModGrid::new(c.range(), |i| { 
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

impl<C> Index<C::Index> for Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>,
    C::Output: RModStr<R = C::R>
{
    type Output = GenericRModStr<C::R>;

    fn index(&self, index: C::Index) -> &Self::Output {
        &self.grid[index]
    }
}

impl<C> GradedRModStr for Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>,
    C::Output: RModStr<R = C::R>
{
    type R = C::R;
    type Index = C::Index;
    type IndexRange = C::IndexRange;

    fn in_range(&self, k: Self::Index) -> bool {
        self.original.in_range(k)
    }

    fn range(&self) -> Self::IndexRange {
        self.original.range()
    }
}

impl<C> ChainComplex for Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>,
    C::Output: RModStr<R = C::R>
{
    fn d_degree(&self) -> Self::Index {
        self.original.d_degree()
    }

    fn d_matrix(&self, k: Self::Index) -> SpMat<Self::R> {
        if self.in_range(k) { 
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
    use crate::math::homology::base::RModStr;
    use crate::math::homology::homology::GenericHomology;
    use super::*;
    use super::super::complex::tests::*;

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