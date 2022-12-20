use std::fmt::Display;
use std::ops::{RangeInclusive, Index};

use crate::math::homology::homology::{Homology, HomologySummand};
use crate::math::homology::homology::GenericHomology;
use crate::math::homology::reduce::Reduced;
use crate::math::matrix::CsMatElem;
use crate::math::traits::{EucRing, EucRingOps};
use crate::links::Link;
use super::complex::KhComplex;

pub struct KhHomology<R> 
where R: EucRing + CsMatElem, for<'x> &'x R: EucRingOps<R> { 
    homology: GenericHomology<R>,
    shift: (isize, isize)
}

impl<R> KhHomology<R> 
where R: EucRing + CsMatElem, for<'x> &'x R: EucRingOps<R> { 
    pub fn new(l: &Link) -> Self {
        Self::new_ht(l, R::zero(), R::zero())
    }

    pub fn new_ht(l: &Link, h: R, t: R) -> Self {
        Self::_new(l, h, t, true)
    }

    fn _new(l: &Link, h: R, t: R, reduction: bool) -> Self {
        let complex = KhComplex::new_ht(l, h, t);
        let shift = complex.shift();
        let homology = if reduction { 
            let reduced = Reduced::new(complex);
            GenericHomology::from(reduced)
        } else { 
            GenericHomology::from(complex)
        };
        Self { homology, shift }
    }

    pub fn shift(&self) -> (isize, isize) { 
        self.shift
    }
}

impl<R> Index<isize> for KhHomology<R> 
where 
    R: EucRing + CsMatElem,
    for<'x> &'x R: EucRingOps<R>
{
    type Output = HomologySummand<R>;

    fn index(&self, index: isize) -> &Self::Output {
        self.homology.index(index)
    }
}

impl<R> Homology for KhHomology<R> 
where 
    R: EucRing + CsMatElem,
    for<'x> &'x R: EucRingOps<R>
{
    type R = R;

    fn range(&self) -> RangeInclusive<isize> {
        self.homology.range()
    }
}

impl<R> Display for KhHomology<R> 
where 
    R: EucRing + CsMatElem,
    for<'x> &'x R: EucRingOps<R>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.homology.fmt(f)
    }
}

#[cfg(test)]
mod tests {
    use crate::links::Link;
    use crate::math::homology::homology::Homology;
    use super::KhHomology;
    
    #[test]
    fn kh_empty() {
        let l = Link::empty();
        let h = KhHomology::<i32>::_new(&l, 0, 0, false);

        assert_eq!(h.range(), 0..=0);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let h = KhHomology::<i32>::_new(&l, 0, 0, false);

        assert_eq!(h.range(), 0..=0);
        
        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[0].is_free(), true);
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let h = KhHomology::<i32>::_new(&l, 0, 0, false);

        assert_eq!(h.range(), -3..=0);

        assert_eq!(h[-3].rank(), 1);
        assert_eq!(h[-3].is_free(), true);

        assert_eq!(h[-2].rank(), 1);
        assert_eq!(h[-2].tors(), &vec![2]);

        assert_eq!(h[-1].is_zero(), true);

        assert_eq!(h[ 0].rank(), 2);
        assert_eq!(h[ 0].is_free(), true);
    }

    #[test]
    fn kh_trefoil_mirror() {
        let l = Link::trefoil().mirror();
        let h = KhHomology::<i32>::_new(&l, 0, 0, false);

        assert_eq!(h.range(), 0..=3);

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].is_zero(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);

        assert_eq!(h[3].rank(), 1);
        assert_eq!(h[3].tors(), &vec![2]);
    }

    #[test]
    fn kh_figure8() {
        let l = Link::figure8();
        let h = KhHomology::<i32>::_new(&l, 0, 0, false);

        assert_eq!(h.range(), -2..=2);

        assert_eq!(h[-2].rank(), 1);
        assert_eq!(h[-2].is_free(), true);

        assert_eq!(h[-1].rank(), 1);
        assert_eq!(h[-1].tors(), &vec![2]);

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 1);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].tors(), &vec![2]);
    }
}