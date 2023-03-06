use std::collections::HashMap;
use std::fmt::Display;
use std::iter::Rev;
use std::ops::{Index, RangeInclusive};
use yui_matrix::sparse::*;
use yui_core::{RingOps, Ring};
use crate::utils::ChainReducer;
use crate::{GridItr, GridIdx, RModStr, RModGrid, ChainComplex, GenericRModStr, GenericRModGrid};

pub struct GenericChainComplex<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr,
    I::Item: GridIdx
{
    grid: GenericRModGrid<GenericRModStr<R>, I>,
    d_degree: I::Item,
    d_matrices: HashMap<I::Item, SpMat<R>>
}

impl<R, I> GenericChainComplex<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr,
    I::Item: GridIdx
{
    pub fn new(range: I, d_degree: I::Item, d_matrices: HashMap<I::Item, SpMat<R>>) -> Self {
        let grid = GenericRModGrid::new(range, |i| {
            if let Some(d) = d_matrices.get(&i) {
                let n = d.cols();
                Some(GenericRModStr::new(n, vec![]))
            } else if let Some(d) = d_matrices.get(&(i - d_degree)) {
                let n = d.rows();
                Some(GenericRModStr::new(n, vec![]))
            } else {
                None
            }
        });

        Self { grid, d_degree, d_matrices }
    }

    pub fn generate<F>(range: I, d_degree: I::Item, mut d_matrices: F) -> Self 
    where F: FnMut(I::Item) -> Option<SpMat<R>> {
        let d_matrices = range.clone().flat_map(|i| 
            if let Some(d_i) = d_matrices(i) { 
                Some((i, d_i))
            } else { 
                None
            }
        ).collect();

        Self::new(range, d_degree, d_matrices)
    }

    pub fn simplify(&self) -> GenericChainComplex<R, I> { 
        let mut red = ChainReducer::new(self);
        red.process();

        Self::generate(
            self.indices(),
            self.d_degree(),
            |i| Some(red.take_matrix(i))
        )
    }
}

impl<R, I> GenericChainComplex<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr + DoubleEndedIterator,
    I::Item: GridIdx
{
    pub fn dual(&self) -> GenericChainComplex<R, Rev<I>> {
        let range = self.indices().rev();
        let d_degree = -self.d_degree();
        let d_matrices = self.d_matrices.iter().map(|(&i, d)| 
            (i - d_degree, d.transpose().to_owned())
        ).collect();
        GenericChainComplex::new(range, d_degree, d_matrices)
    }
}

impl<R, I> Display for GenericChainComplex<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr,
    I::Item: GridIdx
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt_default(f, "C")
    }
}

impl<R, I> Index<I::Item> for GenericChainComplex<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr,
    I::Item: GridIdx
{
    type Output = GenericRModStr<R>;

    fn index(&self, index: I::Item) -> &Self::Output {
        &self.grid[index]
    }
}

impl<R, I> RModGrid for GenericChainComplex<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr,
    I::Item: GridIdx
{
    type R = R;
    type Idx = I::Item;
    type IdxIter = I;
    
    fn contains_idx(&self, k: Self::Idx) -> bool {
        self.grid.contains_idx(k)
    }

    fn indices(&self) -> Self::IdxIter {
        self.grid.indices()
    }
}

impl<R, I> ChainComplex for GenericChainComplex<R, I> 
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr,
    I::Item: GridIdx
{
    fn d_degree(&self) -> Self::Idx {
        self.d_degree
    }

    fn d_matrix(&self, k: Self::Idx) -> SpMat<R> {
        if let Some(d) = self.d_matrices.get(&k) { 
            d.clone()
        } else {
            let m = self[k + self.d_degree].rank();
            let n = self[k].rank();
            SpMat::zero((m, n))
        }
    }
}

impl<R> GenericChainComplex<R, Rev<RangeInclusive<isize>>> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn descending(d_matrices: Vec<SpMat<R>>) -> Self {
        let n = d_matrices.len() as isize;
        let range = (0..=n).rev();
        let d_degree = -1;
        let d_matrices = d_matrices.into_iter().enumerate().map(|(i, d)| 
            ((i + 1) as isize, d)
        ).collect();
        Self::new(range, d_degree, d_matrices)
    }
}

impl<R> GenericChainComplex<R, RangeInclusive<isize>>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn ascending(d_matrices: Vec<SpMat<R>>) -> Self {
        let n = d_matrices.len() as isize;
        let range = 0..=n;
        let d_degree = 1;
        let d_matrices = d_matrices.into_iter().enumerate().map(|(i, d)| 
            (i as isize, d)
        ).collect();
        GenericChainComplex::new(range, d_degree, d_matrices)
    }
}

#[cfg(test)]
mod tests { 
    use super::*;
    use crate::test::{TestChainComplex, ChainComplexValidation};

    #[test]
    fn zero_complex() {
        let c = TestChainComplex::<i32>::zero();
        assert_eq!(c.d_degree(), -1);
        assert_eq!(c[0].rank(), 0);
    }

    #[test]
    fn d3() {
        let c = TestChainComplex::<i32>::d3();

        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 4);
        assert_eq!(c[3].rank(), 1);

        c.check_d_all();
    }

    #[test]
    fn s2() {
        let c = TestChainComplex::<i32>::s2();

        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 4);

        c.check_d_all();
    }

    #[test]
    fn t2() {
        let c = TestChainComplex::<i32>::t2();

        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c[0].rank(), 9);
        assert_eq!(c[1].rank(), 27);
        assert_eq!(c[2].rank(), 18);

        c.check_d_all();
    }

    #[test]
    fn rp2() {
        let c = TestChainComplex::<i32>::rp2();

        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c[0].rank(), 6);
        assert_eq!(c[1].rank(), 15);
        assert_eq!(c[2].rank(), 10);

        c.check_d_all();
    }

    #[test]
    fn d3_retr() {
        let c = TestChainComplex::<i32>::d3().simplify();

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 0);
        assert_eq!(c[2].rank(), 0);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();
    }

    #[test]
    fn s2_retr() {
        let c = TestChainComplex::<i32>::s2().simplify();

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 0);
        assert_eq!(c[2].rank(), 1);

        c.check_d_all();
    }

    #[test]
    fn t2_retr() {
        let c = TestChainComplex::<i32>::t2().simplify();

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 2);
        assert_eq!(c[2].rank(), 1);

        c.check_d_all();
    }

    #[test]
    fn rp2_retr() {
        let c = TestChainComplex::<i32>::rp2().simplify();

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 1);
        assert_eq!(c[2].rank(), 1);

        c.check_d_all();
    }

    #[test]
    fn d3_dual() {
        let c = TestChainComplex::<i32>::d3().dual();

        assert_eq!(c.d_degree(), 1);
        
        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 4);
        assert_eq!(c[3].rank(), 1);

        c.check_d_all();
    }

    #[test]
    fn s2_dual() {
        let c = TestChainComplex::<i32>::s2().dual();

        assert_eq!(c.d_degree(), 1);
        
        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 4);

        c.check_d_all();
    }

    #[test]
    fn t2_dual() {
        let c = TestChainComplex::<i32>::t2().dual();

        assert_eq!(c.d_degree(), 1);
        
        assert_eq!(c[0].rank(), 9);
        assert_eq!(c[1].rank(), 27);
        assert_eq!(c[2].rank(), 18);

        c.check_d_all();
    }

    #[test]
    fn rp2_dual() {
        let c = TestChainComplex::<i32>::rp2().dual();

        assert_eq!(c.d_degree(), 1);
        
        assert_eq!(c[0].rank(), 6);
        assert_eq!(c[1].rank(), 15);
        assert_eq!(c[2].rank(), 10);

        c.check_d_all();
    }

}