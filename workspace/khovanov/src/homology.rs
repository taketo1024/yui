use std::fmt::Display;
use std::ops::{RangeInclusive, Index};

use yui_homology::{Idx2, Idx2Range};
use yui_homology::{GradedRModStr, GenericRModStr, Reduced, Homology, GenericHomology};
use yui_core::{EucRing, EucRingOps, Ring, RingOps};
use yui_link::Link;
use super::complex::{KhComplex, KhComplexBigraded};

pub struct KhHomology<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    homology: GenericHomology<R, RangeInclusive<isize>>
}

impl<R> KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(l: Link, h: R, t: R, reduced: bool) -> Self {
        let complex = KhComplex::new(l, h, t, reduced);
        Self::from(complex)
    }

    pub fn unreduced(l: Link) -> Self { 
        Self::new(l, R::zero(), R::zero(), false)
    }

    pub fn unreduced_ht(l: Link, h: R, t: R) -> Self { 
        Self::new(l, h, t, false)
    }

    pub fn reduced(l: Link) -> Self { 
        Self::new(l, R::zero(), R::zero(), true)
    }

    pub fn reduced_ht(l: Link, h: R, t: R) -> Self { 
        Self::new(l, h, t, true)
    }
}

impl<R> From<KhComplex<R>> for KhHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    fn from(c: KhComplex<R>) -> Self {
        let reduced = Reduced::from(c);
        let homology = GenericHomology::from(reduced);
        Self { homology }
    }
}

impl<R> Index<isize> for KhHomology<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = GenericRModStr<R>;

    fn index(&self, index: isize) -> &Self::Output {
        self.homology.index(index)
    }
}

impl<R> GradedRModStr for KhHomology<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;
    type Index = isize;
    type IndexRange = RangeInclusive<isize>;

    fn in_range(&self, k: Self::Index) -> bool {
        self.homology.in_range(k)
    }

    fn range(&self) -> Self::IndexRange {
        self.homology.range()
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
    homology: GenericHomology<R, Idx2Range>
}

impl<R> KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(l: Link, reduced: bool) -> Self {
        let c = KhComplexBigraded::new(l, reduced);
        Self::from(c)
    }

    pub fn unreduced(l: Link) -> Self {
        Self::new(l, false)
    }

    pub fn reduced(l: Link) -> Self {
        Self::new(l, true)
    }
}

impl<R> From<KhComplexBigraded<R>> for KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    fn from(c: KhComplexBigraded<R>) -> Self {
        let reduced = Reduced::from(c);
        let homology = GenericHomology::from(reduced);
        Self { homology }
    }
}

impl<R> Index<Idx2> for KhHomologyBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = GenericRModStr<R>;

    fn index(&self, index: Idx2) -> &Self::Output {
        &self.homology[index]
    }
}

impl<R> GradedRModStr for KhHomologyBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;
    type Index = Idx2;
    type IndexRange = Idx2Range;

    fn in_range(&self, k: Self::Index) -> bool {
        self.homology.in_range(k)
    }

    fn range(&self) -> Self::IndexRange {
        self.homology.range()
    }
}

impl<R> Homology for KhHomologyBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {}

impl<R> Display for KhHomologyBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt_default(f)
    }
}

#[cfg(test)]
mod tests {
    use yui_link::Link;
    use yui_homology::{RModStr, GradedRModStr};
    use super::KhHomology;
    
    #[test]
    fn kh_empty() {
        let l = Link::empty();
        let h = KhHomology::<i32>::unreduced(l);

        assert_eq!(h.range(), 0..=0);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let h = KhHomology::<i32>::unreduced(l);

        assert_eq!(h.range(), 0..=0);
        
        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[0].is_free(), true);
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let h = KhHomology::<i32>::unreduced(l);

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
        let h = KhHomology::<i32>::unreduced(l);

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
        let h = KhHomology::<i32>::unreduced(l);

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