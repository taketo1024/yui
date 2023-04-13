use std::collections::HashMap;
use std::fmt::Display;
use std::iter::Rev;
use std::ops::Index;
use yui_matrix::sparse::*;
use yui_core::{RingOps, Ring, EucRingOps, EucRing};
use crate::utils::{ChainReducer, HomologyCalc};
use crate::{GridItr, GridIdx, RModStr, RModGrid, ChainComplex, GenericRModStr, GenericRModGrid, HomologyComputable, GenericHomology, Grid};

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
    pub fn new<Iter>(range: I, d_degree: I::Item, d_matrices: Iter) -> Self
    where Iter: Iterator<Item = (I::Item, SpMat<R>)> {
        let d_matrices = d_matrices.collect::<HashMap<_, _>>();
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
        );

        Self::new(range, d_degree, d_matrices)
    }

    pub fn simplify(&self) -> GenericChainComplex<R, I> { 
        self.simplify_with(vec![]).0
    }

    pub fn simplify_with(&self, vecs: Vec<(I::Item, SpVec<R>)>) -> (GenericChainComplex<R, I>, Vec<(I::Item, SpVec<R>)>) { 
        let mut red = ChainReducer::new(self);
        for (i, v) in vecs { 
            red.set_vec(i, v);
        }
        red.process();

        let vecs = red.take_all_vecs();

        let c = Self::generate(
            self.indices(),
            self.d_degree(),
            |i| Some(red.take_matrix(i))
        );

        (c, vecs)
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
        );
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

impl<R, I> Grid for GenericChainComplex<R, I>
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

impl<R, I> HomologyComputable for GenericChainComplex<R, I>
where 
    R: EucRing, for<'x> &'x R: EucRingOps<R>,
    I: GridItr,
    I::Item: GridIdx
{
    type Homology = GenericHomology<R, I>;
    type HomologySummand = GenericRModStr<R>;

    fn homology_at(&self, i: Self::Idx) -> Self::HomologySummand {
        let d1 = self.d_matrix(i - self.d_degree());
        let d2 = self.d_matrix(i);
        HomologyCalc::calculate(d1, d2)
    }

    fn homology(&self) -> Self::Homology {
        let range = self.indices();
        let grid = GenericRModGrid::new(range, |i| {
            let h_i = self.homology_at(i);
            if !h_i.is_zero() { Some(h_i) } else { None }
        });
        Self::Homology::new(grid)
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