use std::ops::RangeInclusive;

use crate::math::traits::{Ring, RingOps};
use crate::math::matrix::CsMatElem;
use crate::math::homology::complex::{ChainComplex, Graded};
use crate::links::Link;
use super::algebra::KhAlgStr;
use super::cube::{KhEnhState, KhCube};

pub struct KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    cube: KhCube<R>,
    shift: (isize, isize)
}

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(l: &Link) -> Self { 
        Self::new_ht(l, R::zero(), R::zero())
    }

    pub fn new_ht(l: &Link, h: R, t: R) -> Self { 
        let str = KhAlgStr::new(h, t);
        let cube = KhCube::new(l, str);
        let (n_pos, n_neg) = l.signed_crossing_nums();
        let (n_pos, n_neg) = (n_pos as isize, n_neg as isize);
        let shift = (-n_neg, n_pos - 2 * n_neg);
        KhComplex { cube, shift }
    }

    pub fn shift(&self) -> (isize, isize) { 
        self.shift
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
        let n = self.cube.dim() as isize;
        let s = self.shift.0;
        s ..= s + n
    }
}

impl<R> ChainComplex for KhComplex<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> { 
    type R = R;
    type Generator = KhEnhState;

    fn d_degree(&self) -> Self::Index {
        1
    }

    fn generators(&self, k: Self::Index) -> Vec<&Self::Generator> {
        let s = self.shift.0;
        let k = (k - s) as usize;
        self.cube.generators(k)
    }

    fn differentiate(&self, _k: Self::Index, x:&Self::Generator) -> Vec<(Self::Generator, Self::R)> {
        self.cube.differentiate(x)
    }

    fn d_matrix(&self, k: Self::Index) -> sprs::CsMat<Self::R> {
        self.impl_d_matrix_from_differentiate(k)
    }
}

#[cfg(test)]
mod tests {
    use crate::links::Link;
    use crate::math::homology::complex::{ChainComplex, Graded};
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