use std::ops::RangeInclusive;

use crate::math::homology::chain_complex::ChainComplex;
use crate::math::homology::homology::{HomologyComputable, SimpleHomology};
use crate::math::matrix::CsMatElem;
use crate::math::traits::{Ring, RingOps, EucRing, EucRingOps};
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

impl<R> ChainComplex for KhComplex<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> { 
    type R = R;
    type Generator = KhEnhState;

    fn d_degree(&self) -> isize {
        1
    }

    fn hdeg_range(&self) -> RangeInclusive<isize> {
        let n = self.cube.dim() as isize;
        let s = self.shift.0;
        s ..= s + n
    }

    fn generators(&self, k: isize) -> Vec<&Self::Generator> {
        let s = self.shift.0;
        let k = (k - s) as usize;
        self.cube.generators(k)
    }

    fn differentiate(&self, _k: isize, x:&Self::Generator) -> Vec<(Self::Generator, Self::R)> {
        self.cube.differentiate(x)
    }
}

impl<R> HomologyComputable for KhComplex<R>
where R: EucRing + CsMatElem, for<'x> &'x R: EucRingOps<R> { 
    type Homology = SimpleHomology<R>;
    fn homology(&self) -> Self::Homology {
        SimpleHomology::from(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::{links::Link, khovanov::algebra::KhAlgStr, math::homology::chain_complex::ChainComplex};
    use super::KhComplex;

    #[test]
    fn kh_empty() {
        let l = Link::empty();
        let c = KhComplex::<i32>::new(&l);

        assert_eq!(c.hdeg_range(), 0..=0);
        c.check_d_all();
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let c = KhComplex::<i32>::new(&l);

        assert_eq!(c.hdeg_range(), 0..=0);
        c.check_d_all();
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let c = KhComplex::<i32>::new(&l);

        assert_eq!(c.hdeg_range(), -3..=0);
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

        assert_eq!(c.hdeg_range(), -2..=2);
        assert_eq!(c.generators(-2).len(), 8);
        assert_eq!(c.generators(-1).len(), 16);
        assert_eq!(c.generators( 0).len(), 18);
        assert_eq!(c.generators( 1).len(), 16);
        assert_eq!(c.generators( 2).len(), 8);

        c.check_d_all();
    }
}