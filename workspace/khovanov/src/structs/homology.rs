use std::fmt::Display;
use std::ops::{RangeInclusive, Index};

use yui_homology::{Idx2, Idx2Iter, ChainComplex, Grid, GenericRModStr, Homology, GenericHomology, fmt::FmtList};
use yui_core::{EucRing, EucRingOps, Ring, RingOps};
use yui_link::Link;
use super::complex::{KhComplex, KhComplexBigraded};

pub type KhHomologySummand<R> = GenericRModStr<R>;

pub struct KhHomology<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    homology: GenericHomology<R, RangeInclusive<isize>>
}

impl<R> KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        let complex = KhComplex::new(l, h, t, reduced);
        Self::from(&complex)
    }
}

impl<R> From<&KhComplex<R>> for KhHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    fn from(c: &KhComplex<R>) -> Self {
        let c = c.as_generic().simplify();
        let homology = GenericHomology::from(c);
        Self { homology }
    }
}

impl<R> Index<isize> for KhHomology<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = KhHomologySummand<R>;

    fn index(&self, index: isize) -> &Self::Output {
        self.homology.index(index)
    }
}

impl<R> Grid for KhHomology<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Idx = isize;
    type IdxIter = RangeInclusive<isize>;
    type Output = KhHomologySummand<R>;

    fn contains_idx(&self, k: Self::Idx) -> bool {
        self.homology.contains_idx(k)
    }

    fn indices(&self) -> Self::IdxIter {
        self.homology.indices()
    }

    fn get(&self, i: Self::Idx) -> Option<&Self::Output> {
        self.homology.get(i)
    }
}

impl<R> Homology for KhHomology<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {}

impl<R> Display for KhHomology<R> 
where R: Ring, for<'x> &'x R: EucRingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.homology.fmt(f)
    }
}

pub struct KhHomologyBigraded<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    homology: GenericHomology<R, Idx2Iter>
}

impl<R> KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(l: Link, reduced: bool) -> Self {
        let c = KhComplexBigraded::new(l, reduced);
        Self::from(&c)
    }
}

impl<R> From<&KhComplexBigraded<R>> for KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    fn from(c: &KhComplexBigraded<R>) -> Self {
        let c = c.as_generic().simplify();
        let homology = GenericHomology::from(c);
        Self { homology }
    }
}

impl<R> Index<[isize; 2]> for KhHomologyBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = GenericRModStr<R>;

    fn index(&self, index: [isize; 2]) -> &Self::Output {
        let index = Idx2::from(index);
        &self.homology[index]
    }
}

impl<R> Grid for KhHomologyBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Idx = Idx2;
    type IdxIter = Idx2Iter;
    type Output = GenericRModStr<R>;

    fn contains_idx(&self, k: Self::Idx) -> bool {
        self.homology.contains_idx(k)
    }

    fn indices(&self) -> Self::IdxIter {
        self.homology.indices()
    }

    fn get(&self, i: Self::Idx) -> Option<&Self::Output> {
        self.homology.get(i)
    }
}

impl<R> Homology for KhHomologyBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {}

impl<R> Display for KhHomologyBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt_list(f)
    }
}

#[cfg(test)]
mod tests {
    use yui_link::Link;
    use yui_homology::{RModStr, Grid};
    use super::*;
    
    #[test]
    fn kh_empty() {
        let l = Link::empty();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.indices(), 0..=0);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.indices(), 0..=0);
        
        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[0].is_free(), true);
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.indices(), -3..=0);

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
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.indices(), 0..=3);

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
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.indices(), -2..=2);

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

    #[test]
    fn kh_empty_bigr() {
        let l = Link::empty();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[[0,0]].rank(), 1);
        assert_eq!(h[[0,0]].is_free(), true);
    }

    #[test]
    fn kh_unknot_bigr() {
        let l = Link::unknot();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[[0,-1]].rank(), 1);
        assert_eq!(h[[0,-1]].is_free(), true);
        assert_eq!(h[[0, 1]].rank(), 1);
        assert_eq!(h[[0, 1]].is_free(), true);
    }

    #[test]
    fn kh_unknot_bigr_red() {
        let l = Link::unknot();
        let h = KhHomologyBigraded::<i32>::new(l, true);

        assert_eq!(h[[0, 0]].rank(), 1);
        assert_eq!(h[[0, 0]].is_free(), true);
    }

    #[test]
    fn kh_trefoil_bigr() {
        let l = Link::trefoil();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[[-3,-9]].rank(), 1);
        assert_eq!(h[[-3,-9]].is_free(), true);
        assert_eq!(h[[-2,-7]].rank(), 0);
        assert_eq!(h[[-2,-7]].tors(), &vec![2]);
        assert_eq!(h[[-2,-5]].rank(), 1);
        assert_eq!(h[[-2,-5]].is_free(), true);
        assert_eq!(h[[ 0,-3]].rank(), 1);
        assert_eq!(h[[ 0,-3]].is_free(), true);
        assert_eq!(h[[ 0,-1]].rank(), 1);
        assert_eq!(h[[ 0,-1]].is_free(), true);
    }

    #[test]
    fn kh_trefoil_mirror_bigr() {
        let l = Link::trefoil().mirror();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[[0, 1]].rank(), 1);
        assert_eq!(h[[0, 1]].is_free(), true);
        assert_eq!(h[[0, 3]].rank(), 1);
        assert_eq!(h[[0, 3]].is_free(), true);
        assert_eq!(h[[2, 5]].rank(), 1);
        assert_eq!(h[[2, 5]].is_free(), true);
        assert_eq!(h[[3, 7]].rank(), 0);
        assert_eq!(h[[3, 7]].tors(), &vec![2]);
        assert_eq!(h[[3, 9]].rank(), 1);
        assert_eq!(h[[3, 9]].is_free(), true);
    }

    #[test]
    fn kh_trefoil_bigr_red() {
        let l = Link::trefoil();
        let h = KhHomologyBigraded::<i32>::new(l, true);

        assert_eq!(h[[-3,-8]].rank(), 1);
        assert_eq!(h[[-3,-8]].is_free(), true);
        assert_eq!(h[[-2,-6]].rank(), 1);
        assert_eq!(h[[-2,-6]].is_free(), true);
        assert_eq!(h[[ 0,-2]].rank(), 1);
        assert_eq!(h[[ 0,-2]].is_free(), true);
    }

    #[test]
    fn kh_figure8_bigr() {
        let l = Link::figure8();
        let h = KhHomologyBigraded::<i32>::new(l, false);

        assert_eq!(h[[-2,-5]].rank(), 1);
        assert_eq!(h[[-2,-5]].is_free(), true);
        assert_eq!(h[[-1,-3]].rank(), 0);
        assert_eq!(h[[-1,-3]].tors(), &vec![2]);
        assert_eq!(h[[-1,-1]].rank(), 1);
        assert_eq!(h[[-1,-1]].is_free(), true);
        assert_eq!(h[[ 0,-1]].rank(), 1);
        assert_eq!(h[[ 0,-1]].is_free(), true);
        assert_eq!(h[[ 0, 1]].rank(), 1);
        assert_eq!(h[[ 0, 1]].is_free(), true);
        assert_eq!(h[[ 1, 1]].rank(), 1);
        assert_eq!(h[[ 1, 1]].is_free(), true);
        assert_eq!(h[[ 2, 3]].rank(), 0);
        assert_eq!(h[[ 2, 3]].tors(), &vec![2]);
        assert_eq!(h[[ 2, 5]].rank(), 1);
        assert_eq!(h[[ 2, 5]].is_free(), true);
   }

    #[test]
    fn kh_figure8_bigr_red() {
        let l = Link::figure8();
        let h = KhHomologyBigraded::<i32>::new(l, true);

        assert_eq!(h[[-2,-4]].rank(), 1);
        assert_eq!(h[[-2,-4]].is_free(), true);
        assert_eq!(h[[-1,-2]].rank(), 1);
        assert_eq!(h[[-1,-2]].is_free(), true);
        assert_eq!(h[[ 0, 0]].rank(), 1);
        assert_eq!(h[[ 0, 0]].is_free(), true);
        assert_eq!(h[[ 1, 2]].rank(), 1);
        assert_eq!(h[[ 1, 2]].is_free(), true);
        assert_eq!(h[[ 2, 4]].rank(), 1);
        assert_eq!(h[[ 2, 4]].is_free(), true);
   }
}