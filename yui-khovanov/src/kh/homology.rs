use std::ops::{RangeInclusive, Index};
use delegate::delegate;
use cartesian::cartesian;

use yui_homology::{isize2, Grid2, GridTrait, Homology, Homology2, Summand};
use yui::{EucRing, EucRingOps};
use yui_link::Link;

use crate::kh::KhGen;
use crate::misc::{collect_gen_info, range_of};

use super::{KhChain, KhComplex, KhComplexBigraded};

#[derive(Clone)]
pub struct KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    inner: Homology<KhGen, R>,
    ht: (R, R),
    deg_shift: (isize, isize),
    reduced: bool,
    canon_cycles: Vec<KhChain<R>>
}

impl<R> KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        let c = KhComplex::new(&l, h, t, reduced); 
        Self::from(&c)
    }
    
    pub(crate) fn new_impl(inner: Homology<KhGen, R>, ht: (R, R), deg_shift: (isize, isize), reduced: bool, canon_cycles: Vec<KhChain<R>>) -> Self { 
        Self { inner, ht, deg_shift, reduced, canon_cycles }
    }

    pub fn ht(&self) -> &(R, R) { 
        &self.ht
    }

    pub fn deg_shift(&self) -> (isize, isize) { 
        self.deg_shift
    }

    pub fn is_reduced(&self) -> bool { 
        self.reduced
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        range_of(self.support())
    }

    pub fn canon_cycles(&self) -> &Vec<KhChain<R>> { 
        &self.canon_cycles
    }

    pub fn inner(&self) -> &Homology<KhGen, R> { 
        &self.inner
    }

    pub fn truncated(&self, range: RangeInclusive<isize>) -> Self {
        Self::new_impl(
            self.inner.truncated(range), 
            self.ht.clone(), 
            self.deg_shift, 
            self.reduced, 
            self.canon_cycles.clone()
        )
    }

    // TODO take self instead of &self
    pub fn into_bigraded(&self) -> KhHomologyBigraded<R> { 
        // TODO: check (h, t)

        let ht = self.ht.clone();
        let deg_shift = self.deg_shift;
        let reduced = self.reduced;
        let canon_cycles = self.canon_cycles.clone();

        let table = collect_gen_info(self.inner());
        let h_range = range_of(table.keys().map(|i| i.0));
        let q_range = range_of(table.keys().map(|i| i.1)).step_by(2);
        let support = cartesian!(h_range, q_range.clone()).map(|(i, j)| 
            isize2(i, j)
        );

        let inner = Grid2::generate(support, move |idx| { 
            let i = idx.0;
            let Some(e) = table.get(&idx) else { 
                return Summand::zero()
            };
            
            let (rank, tors, indices) = e;
            let gens = self[i].raw_gens().clone(); 
            let trans = self[i].trans().sub(indices);
            Summand::new(gens, *rank, tors.clone(), trans)
        });

        KhHomologyBigraded::new_impl(inner, ht, deg_shift, reduced, canon_cycles)
    }
}

impl<R> From<&KhComplex<R>> for KhHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    fn from(c: &KhComplex<R>) -> Self {
        KhHomology::new_impl(
            c.inner().reduced().homology(), 
            c.ht().clone(), 
            c.deg_shift(), 
            c.is_reduced(),
            c.canon_cycles().clone()
        )
    }
}

impl<R> GridTrait<isize> for KhHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Support = std::vec::IntoIter<isize>;
    type Item = Summand<KhGen, R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Support;
            fn is_supported(&self, i: isize) -> bool;
            fn get(&self, i: isize) -> &Self::Item;
            fn get_default(&self) -> &Self::Item;
        }
    }
}

impl<R> Index<isize> for KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = Summand<KhGen, R>;

    delegate! { 
        to self.inner { 
            fn index(&self, index: isize) -> &Self::Output;
        }
    }
}

#[derive(Clone)]
pub struct KhHomologyBigraded<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    inner: Homology2<KhGen, R>,
    ht: (R, R),
    deg_shift: (isize, isize),
    reduced: bool,
    canon_cycles: Vec<KhChain<R>>
}

impl<R> KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        let h = KhHomology::new(l, h, t, reduced);
        h.into_bigraded()
    }

    fn new_impl(inner: Homology2<KhGen, R>, ht: (R, R), deg_shift: (isize, isize), reduced: bool, canon_cycles: Vec<KhChain<R>>) -> Self { 
        Self { inner, ht, deg_shift, reduced, canon_cycles }
    }

    pub fn ht(&self) -> &(R, R) { 
        &self.ht
    }

    pub fn deg_shift(&self) -> (isize, isize) { 
        self.deg_shift
    }

    pub fn is_reduced(&self) -> bool { 
        self.reduced
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        range_of(self.support().map(|i| i.0))
    }

    pub fn q_range(&self) -> RangeInclusive<isize> { 
        range_of(self.support().map(|i| i.1))
    }

    pub fn canon_cycles(&self) -> &Vec<KhChain<R>> { 
        &self.canon_cycles
    }

    pub fn inner(&self) -> &Homology2<KhGen, R> { 
        &self.inner
    }
}

impl<R> From<&KhComplexBigraded<R>> for KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    fn from(c: &KhComplexBigraded<R>) -> Self {
        KhHomologyBigraded::new_impl(
            c.inner().reduced().homology(), 
            c.ht().clone(), 
            c.deg_shift(), 
            c.is_reduced(),
            c.canon_cycles().clone()
        )
    }
}

impl<R> GridTrait<isize2> for KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Support = std::vec::IntoIter<isize2>;
    type Item = Summand<KhGen, R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Support;
            fn is_supported(&self, i: isize2) -> bool;
            fn get(&self, i: isize2) -> &Self::Item;
            fn get_default(&self) -> &Self::Item;
        }
    }
}

impl<R> Index<(isize, isize)> for KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = Summand<KhGen, R>;

    delegate! { 
        to self.inner { 
            fn index(&self, index: (isize, isize)) -> &Self::Output;
        }
    }
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use yui::poly::HPoly;
    use yui::FF2;
    use yui_homology::SummandTrait;
    use yui_link::Link;
    use super::*;
    
    #[test]
    fn kh_empty() {
        let l = Link::empty();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);
        
        assert_eq!(h[0].rank(), 2);
        assert!(h[0].is_free());
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), -3..=0);

        assert_eq!(h[-3].rank(), 1);
        assert!(h[-3].is_free());

        assert_eq!(h[-2].rank(), 1);
        assert_eq!(h[-2].tors(), &vec![2]);

        assert!(h[-1].is_zero());

        assert_eq!(h[ 0].rank(), 2);
        assert!(h[ 0].is_free());
    }

    #[test]
    fn kh_trefoil_mirror() {
        let l = Link::trefoil().mirror();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=3);

        assert_eq!(h[0].rank(), 2);
        assert!(h[0].is_free());

        assert!(h[1].is_zero());

        assert_eq!(h[2].rank(), 1);
        assert!(h[2].is_free());

        assert_eq!(h[3].rank(), 1);
        assert_eq!(h[3].tors(), &vec![2]);
    }

    #[test]
    fn kh_figure8() {
        let l = Link::figure8();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), -2..=2);

        assert_eq!(h[-2].rank(), 1);
        assert!(h[-2].is_free());

        assert_eq!(h[-1].rank(), 1);
        assert_eq!(h[-1].tors(), &vec![2]);

        assert_eq!(h[0].rank(), 2);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 1);
        assert!(h[1].is_free());

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].tors(), &vec![2]);
    }

    #[test]
    fn kh_empty_bigr() {
        let l = Link::empty();
        let h = KhHomologyBigraded::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);
        assert_eq!(h[(0,0)].rank(), 1);
        assert!(h[(0,0)].is_free());
    }

    #[test]
    fn kh_unknot_bigr() {
        let l = Link::unknot();
        let h = KhHomologyBigraded::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);
        assert_eq!(h.q_range(), -1..=1);

        assert_eq!(h[(0,-1)].rank(), 1);
        assert!(h[(0,-1)].is_free());
        assert_eq!(h[(0, 1)].rank(), 1);
        assert!(h[(0, 1)].is_free());
    }

    #[test]
    fn kh_unknot_bigr_red() {
        let l = Link::unknot();
        let h = KhHomologyBigraded::new(&l, &0, &0, true);

        assert_eq!(h.h_range(), 0..=0);
        assert_eq!(h.q_range(), 0..=0);

        assert_eq!(h[(0, 0)].rank(), 1);
        assert!(h[(0, 0)].is_free());
    }

    #[test]
    fn kh_trefoil_bigr() {
        let l = Link::trefoil();
        let h = KhHomologyBigraded::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), -3..=0);
        assert_eq!(h.q_range(), -9..=-1);

        assert_eq!(h[(-3,-9)].rank(), 1);
        assert!(h[(-3,-9)].is_free());
        assert_eq!(h[(-2,-7)].rank(), 0);
        assert_eq!(h[(-2,-7)].tors(), &vec![2]);
        assert_eq!(h[(-2,-5)].rank(), 1);
        assert!(h[(-2,-5)].is_free());
        assert_eq!(h[( 0,-3)].rank(), 1);
        assert!(h[( 0,-3)].is_free());
        assert_eq!(h[( 0,-1)].rank(), 1);
        assert!(h[( 0,-1)].is_free());
    }

    #[test]
    fn kh_trefoil_mirror_bigr() {
        let l = Link::trefoil().mirror();
        let h = KhHomologyBigraded::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=3);
        assert_eq!(h.q_range(), 1..=9);

        assert_eq!(h[(0, 1)].rank(), 1);
        assert!(h[(0, 1)].is_free());
        assert_eq!(h[(0, 3)].rank(), 1);
        assert!(h[(0, 3)].is_free());
        assert_eq!(h[(2, 5)].rank(), 1);
        assert!(h[(2, 5)].is_free());
        assert_eq!(h[(3, 7)].rank(), 0);
        assert_eq!(h[(3, 7)].tors(), &vec![2]);
        assert_eq!(h[(3, 9)].rank(), 1);
        assert!(h[(3, 9)].is_free());
    }

    #[test]
    fn kh_trefoil_bigr_red() {
        let l = Link::trefoil();
        let h = KhHomologyBigraded::new(&l, &0, &0, true);

        assert_eq!(h.h_range(), -3..=0);
        assert_eq!(h.q_range(), -8..=-2);

        assert_eq!(h[(-3,-8)].rank(), 1);
        assert!(h[(-3,-8)].is_free());
        assert_eq!(h[(-2,-6)].rank(), 1);
        assert!(h[(-2,-6)].is_free());
        assert_eq!(h[( 0,-2)].rank(), 1);
        assert!(h[( 0,-2)].is_free());
    }

    #[test]
    fn kh_figure8_bigr() {
        let l = Link::figure8();
        let h = KhHomologyBigraded::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), -2..=2);
        assert_eq!(h.q_range(), -5..=5);

        assert_eq!(h[(-2,-5)].rank(), 1);
        assert!(h[(-2,-5)].is_free());
        assert_eq!(h[(-1,-3)].rank(), 0);
        assert_eq!(h[(-1,-3)].tors(), &vec![2]);
        assert_eq!(h[(-1,-1)].rank(), 1);
        assert!(h[(-1,-1)].is_free());
        assert_eq!(h[( 0,-1)].rank(), 1);
        assert!(h[( 0,-1)].is_free());
        assert_eq!(h[( 0, 1)].rank(), 1);
        assert!(h[( 0, 1)].is_free());
        assert_eq!(h[( 1, 1)].rank(), 1);
        assert!(h[( 1, 1)].is_free());
        assert_eq!(h[( 2, 3)].rank(), 0);
        assert_eq!(h[( 2, 3)].tors(), &vec![2]);
        assert_eq!(h[( 2, 5)].rank(), 1);
        assert!(h[( 2, 5)].is_free());
    }

    #[test]
    fn kh_figure8_bigr_red() {
        let l = Link::figure8();
        let h = KhHomologyBigraded::new(&l, &0, &0, true);

        assert_eq!(h.h_range(), -2..=2);
        assert_eq!(h.q_range(), -4..=4);

        assert_eq!(h[(-2,-4)].rank(), 1);
        assert!(h[(-2,-4)].is_free());
        assert_eq!(h[(-1,-2)].rank(), 1);
        assert!(h[(-1,-2)].is_free());
        assert_eq!(h[( 0, 0)].rank(), 1);
        assert!(h[( 0, 0)].is_free());
        assert_eq!(h[( 1, 2)].rank(), 1);
        assert!(h[( 1, 2)].is_free());
        assert_eq!(h[( 2, 4)].rank(), 1);
        assert!(h[( 2, 4)].is_free());
    }

    #[test]
    fn into_bigr() {
        let l = Link::trefoil();
        let (h, t) = (0, 0);
        let kh = KhHomology::new(&l, &h, &t, false).into_bigraded();
 
        assert_eq!(kh.h_range(), -3..=0);
        assert_eq!(kh.q_range(), -9..=-1);
 
        assert_eq!(kh[(-3,-9)].rank(), 1);
        assert!(kh[(-3,-9)].is_free());
        assert_eq!(kh[(-2,-7)].rank(), 0);
        assert_eq!(kh[(-2,-7)].tors(), &vec![2]);
        assert_eq!(kh[(-2,-5)].rank(), 1);
        assert!(kh[(-2,-5)].is_free());
        assert_eq!(kh[( 0,-3)].rank(), 1);
        assert!(kh[( 0,-3)].is_free());
        assert_eq!(kh[( 0,-1)].rank(), 1);
        assert!(kh[( 0,-1)].is_free());
    }
 
    #[test]
    fn into_bigr_bn() {
        type R = FF2;
        type P = HPoly<'H', R>;

        let l = Link::trefoil();
        let (h, t) = (P::variable(), P::zero());
        let kh = KhHomology::new(&l, &h, &t, false).into_bigraded();

        assert_eq!(kh.h_range(), -2..=0);
        assert_eq!(kh.q_range(), -7..=-1);

        assert_eq!(kh[(-2,-7)].rank(), 0);
        assert_eq!(kh[(-2,-7)].tors(), &vec![h.clone()]);
        assert_eq!(kh[(-2,-5)].rank(), 0);
        assert_eq!(kh[(-2,-5)].tors(), &vec![h.clone()]);
        assert_eq!(kh[( 0,-3)].rank(), 1);
        assert!(kh[( 0,-3)].is_free());
        assert_eq!(kh[( 0,-1)].rank(), 1);
        assert!(kh[( 0,-1)].is_free());
    }
 }