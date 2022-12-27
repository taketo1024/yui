use std::ops::RangeInclusive;

use crate::math::homology::base::Graded;
use crate::math::homology::free::FreeChainComplex;
use crate::math::traits::{Ring, RingOps};
use crate::math::matrix::CsMatElem;
use crate::math::homology::complex::ChainComplex;
use crate::links::Link;
use super::cube::{KhEnhState, KhCube};

pub struct KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    cube: KhCube<R>
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
        KhComplex { cube }
    }
}

impl<R> Graded for KhComplex<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> { 
    type Index = isize;
    type IndexRange = RangeInclusive<isize>;

    fn in_range(&self, k: Self::Index) -> bool {
        self.range().contains(&k)
    }

    fn range(&self) -> Self::IndexRange {
        self.cube.h_range()
    }
}

impl<R> ChainComplex for KhComplex<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> { 
    type R = R;

    fn rank(&self, k: Self::Index) -> usize {
        self.cube.generators(k).len()
    }

    fn d_degree(&self) -> Self::Index {
        1
    }

    fn d_matrix(&self, k: Self::Index) -> sprs::CsMat<Self::R> {
        let source = self.cube.generators(k);
        let target = self.cube.generators(k + 1);
        crate::math::homology::free::make_matrix(&source, &target, |x| self.cube.differentiate(x))
    }
}

impl<R> FreeChainComplex for KhComplex<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> { 
    type Generator = KhEnhState;
        
    fn generators(&self, k: Self::Index) -> Vec<&Self::Generator> {
        self.cube.generators(k)
    }

    fn differentiate(&self, _k: Self::Index, x:&Self::Generator) -> Vec<(Self::Generator, Self::R)> {
        self.cube.differentiate(x)
    }    
}

#[cfg(test)]
mod tests {
    use crate::links::Link;
    use crate::math::homology::base::Graded;
    use crate::math::homology::complex::ChainComplex;
    use crate::math::homology::free::FreeChainComplex;
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
        assert_eq!(c.generators(-3).len(), 8);
        assert_eq!(c.generators(-2).len(), 12);
        assert_eq!(c.generators(-1).len(), 6);
        assert_eq!(c.generators( 0).len(), 4);

        c.check_d_all();
    }

    #[test]
    fn kh_figure8() {
        let l = Link::figure8();
        let c = KhComplex::<i32>::new(&l);

        assert_eq!(c.range(), -2..=2);
        assert_eq!(c.generators(-2).len(), 8);
        assert_eq!(c.generators(-1).len(), 16);
        assert_eq!(c.generators( 0).len(), 18);
        assert_eq!(c.generators( 1).len(), 16);
        assert_eq!(c.generators( 2).len(), 8);

        c.check_d_all();
    }
}