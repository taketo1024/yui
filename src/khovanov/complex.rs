use std::ops::{RangeInclusive, Index};
use std::vec::IntoIter;

use crate::math::homology::base::{GradedRModStr, RModGrid, Idx2};
use crate::math::homology::free::FreeRModStr;
use crate::math::traits::{Ring, RingOps};
use crate::math::matrix::CsMatElem;
use crate::math::homology::complex::ChainComplex;
use crate::links::Link;
use super::cube::{KhEnhState, KhCube};

pub struct KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    cube: KhCube<R>,
    grid: RModGrid<R, FreeRModStr<R, KhEnhState>, isize, RangeInclusive<isize>>
}

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(l: &Link) -> Self { 
        Self::new_ht(l, R::zero(), R::zero())
    }

    pub fn new_ht(l: &Link, h: R, t: R) -> Self { 
        let cube = KhCube::new_ht(l, h, t);
        let grid = RModGrid::new(cube.h_range(), |i| {
            let gens = cube.generators(i);
            FreeRModStr::new(gens)
        });
        KhComplex { cube, grid }
    }
}

impl<R> Index<isize> for KhComplex<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> { 
    type Output = FreeRModStr<R, KhEnhState>;
    
    fn index(&self, index: isize) -> &Self::Output {
        &self.grid[index]
    }
}

impl<R> GradedRModStr for KhComplex<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> { 
    type R = R;
    type Index = isize;
    type IndexRange = RangeInclusive<isize>;

    fn in_range(&self, k: Self::Index) -> bool {
        self.grid.in_range(k)
    }

    fn range(&self) -> Self::IndexRange {
        self.grid.range()
    }
}

impl<R> ChainComplex for KhComplex<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> { 
    fn d_degree(&self) -> Self::Index {
        1
    }

    fn d_matrix(&self, k: Self::Index) -> sprs::CsMat<Self::R> {
        use crate::math::homology::free::make_matrix;
        let source = self.grid[k].generators();
        let target = self.grid[k + 1].generators();
        make_matrix(source, target, |x| self.cube.differentiate(x))
    }
}

pub struct KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    cube: KhCube<R>,
    grid: RModGrid<R, FreeRModStr<R, KhEnhState>, Idx2, IntoIter<Idx2>>
}

impl<R> KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(l: &Link) -> Self { 
        let cube = KhCube::new(l);
        let mut grid1 = RModGrid::new(cube.h_range(), |i| {
            let gens = cube.generators(i);
            FreeRModStr::new(gens)
        });

        let h_range = cube.h_range();
        let q_range = cube.q_range();
        
        let range = Idx2::iterate(
            Idx2(*h_range.start(), *q_range.start()), 
            Idx2(*h_range.end(),   *q_range.end()), 
            (1, 2)
        );

        let j0 = cube.shift().1;
        let grid = RModGrid::new(range, |idx| {
            let (i, j) = idx.as_tuple();
            grid1[i].extract(|x| x.q_deg() == j - j0)
        });

        debug_assert!(grid1.is_zero());

        Self { cube, grid }
    }
}

impl<R> Index<Idx2> for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = FreeRModStr<R, KhEnhState>;

    fn index(&self, index: Idx2) -> &Self::Output {
        &self.grid[index]
    }
}

impl<R> GradedRModStr for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;
    type Index = Idx2;
    type IndexRange = IntoIter<Idx2>;

    fn in_range(&self, k: Self::Index) -> bool {
        self.grid.in_range(k)
    }

    fn range(&self) -> Self::IndexRange {
        self.grid.range()
    }
}

impl<R> ChainComplex for KhComplexBigraded<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> {
    fn d_degree(&self) -> Self::Index {
        Idx2(1, 0)
    }

    fn d_matrix(&self, idx: Self::Index) -> sprs::CsMat<Self::R> {
        use crate::math::homology::free::make_matrix;
        let source = self.grid[idx].generators();
        let target = self.grid[idx + self.d_degree()].generators();
        make_matrix(source, target, |x| self.cube.differentiate(x))
    }
}

#[cfg(test)]
mod tests {
    use crate::links::Link;
    use crate::math::homology::base::GradedRModStr;
    use crate::math::homology::complex::ChainComplex;
    use super::KhComplex;

    #[test]
    fn kh_empty() {
        let l = Link::empty();
        let c = KhComplex::<i32>::new(&l);

        assert_eq!(c.range(), 0..=0);
        c.check_d_all();
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let c = KhComplex::<i32>::new(&l);

        assert_eq!(c.range(), 0..=0);
        c.check_d_all();
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let c = KhComplex::<i32>::new(&l);

        assert_eq!(c.range(), -3..=0);
        assert_eq!(c[-3].generators().len(), 8);
        assert_eq!(c[-2].generators().len(), 12);
        assert_eq!(c[-1].generators().len(), 6);
        assert_eq!(c[ 0].generators().len(), 4);

        c.check_d_all();
    }

    #[test]
    fn kh_figure8() {
        let l = Link::figure8();
        let c = KhComplex::<i32>::new(&l);

        assert_eq!(c.range(), -2..=2);
        assert_eq!(c[-2].generators().len(), 8);
        assert_eq!(c[-1].generators().len(), 16);
        assert_eq!(c[ 0].generators().len(), 18);
        assert_eq!(c[ 1].generators().len(), 16);
        assert_eq!(c[ 2].generators().len(), 8);

        c.check_d_all();
    }
}