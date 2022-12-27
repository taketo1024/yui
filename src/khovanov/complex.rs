use std::ops::{RangeInclusive, Index};

use crate::math::homology::base::{GradedRModStr, RModGrid};
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
        Self::from_cube(cube)
    }

    pub fn from_cube(cube: KhCube<R>) -> Self { 
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
    fn rank(&self, _k: Self::Index) -> usize {
        todo!("remove")
    }

    fn d_degree(&self) -> Self::Index {
        1
    }

    fn d_matrix(&self, k: Self::Index) -> sprs::CsMat<Self::R> {
        let source = self.grid[k].generators();
        let target = self.grid[k + 1].generators();
        crate::math::homology::free::make_matrix(source, target, |x| self.cube.differentiate(x))
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