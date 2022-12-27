use std::collections::HashMap;
use std::fmt::Display;
use std::ops::{RangeInclusive, Index};

use crate::math::homology::base::{Graded, Summand};
use crate::math::homology::homology::{Homology, HomologySummand};
use crate::math::homology::homology::GenericHomology;
use crate::math::homology::reduce::Reduced;
use crate::math::matrix::CsMatElem;
use crate::math::traits::{EucRing, EucRingOps};
use crate::links::Link;
use super::complex::KhComplex;
use super::cube::KhCube;

pub struct KhHomology<R> 
where 
    R: EucRing + CsMatElem, 
    for<'x> &'x R: EucRingOps<R> 
{ 
    homology: GenericHomology<Reduced<KhComplex<R>>>,
}

impl<R> KhHomology<R> 
where 
    R: EucRing + CsMatElem, 
    for<'x> &'x R: EucRingOps<R> 
{ 
    pub fn new(l: &Link) -> Self {
        Self::new_ht(l, R::zero(), R::zero())
    }

    pub fn new_ht(l: &Link, h: R, t: R) -> Self {
        let complex = KhComplex::new_ht(l, h, t);
        Self::from(complex)
    }
}

impl<R> From<KhComplex<R>> for KhHomology<R>
where 
    R: EucRing + CsMatElem, 
    for<'x> &'x R: EucRingOps<R> 
{
    fn from(c: KhComplex<R>) -> Self {
        let reduced = Reduced::from(c);
        let homology = GenericHomology::from(reduced);
        Self { homology }
    }
}

impl<R> Graded for KhHomology<R>
where 
    R: EucRing + CsMatElem, 
    for<'x> &'x R: EucRingOps<R> 
{ 
    type Index = isize;
    type IndexRange = RangeInclusive<isize>;

    fn in_range(&self, k: Self::Index) -> bool {
        self.homology.in_range(k)
    }

    fn range(&self) -> Self::IndexRange {
        self.homology.range()
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

pub struct BigradedKhHomology<R> 
where 
    R: EucRing + CsMatElem, 
    for<'x> &'x R: EucRingOps<R> 
{ 
    decomposed: HashMap<isize, KhHomology<R>>,
    h_range: RangeInclusive<isize>,
    q_range: RangeInclusive<isize>,
    zero: HomologySummand<R>
}

impl<R> BigradedKhHomology<R>
where 
    R: EucRing + CsMatElem, 
    for<'x> &'x R: EucRingOps<R> 
{ 
    pub fn new(l: &Link) -> Self {
        let cube = KhCube::<R>::new(l);

        let h_range = cube.h_range();
        let q_range = cube.q_range();
        let j0 = cube.shift().1;

        let decomposed = q_range.clone().step_by(2).map(|j| {
            let subcube = cube.filter(|x| x.q_deg() == j - j0);
            let c = KhComplex::from_cube(subcube);
            (j, KhHomology::from(c))
        }).collect();

        let zero = HomologySummand::zero();

        Self { decomposed, h_range, q_range, zero }
    }
}

impl<R> Index<[isize; 2]> for BigradedKhHomology<R>
where 
    R: EucRing + CsMatElem, 
    for<'x> &'x R: EucRingOps<R> 
{ 
    type Output = HomologySummand<R>;
    fn index(&self, index: [isize; 2]) -> &Self::Output {
        let [i, j] = index;
        if let Some(h) = self.decomposed.get(&j) { 
            &h[i]
        } else {
            &self.zero
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{links::Link, math::homology::base::{Graded, Summand}};
    use super::KhHomology;
    
    #[test]
    fn kh_empty() {
        let l = Link::empty();
        let h = KhHomology::<i32>::new(&l);

        assert_eq!(h.range(), 0..=0);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let h = KhHomology::<i32>::new(&l);

        assert_eq!(h.range(), 0..=0);
        
        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[0].is_free(), true);
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let h = KhHomology::<i32>::new(&l);

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
        let h = KhHomology::<i32>::new(&l);

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
        let h = KhHomology::<i32>::new(&l);

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