use std::ops::{RangeInclusive, Index};
use std::sync::OnceLock;
use delegate::delegate;

use yui_homology::{ToSeqString, ToTableString, GrMod1, GrMod2, Summand};
use yui_core::{EucRing, EucRingOps, IteratorExt};
use yui_link::Link;

use crate::kh::{KhChainExt, KhState};
use crate::misc::make_gen_grid;

use super::{KhAlg, KhChain, KhComplex};

#[derive(Clone)]
pub struct KhHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    inner: GrMod1<KhState, R>,
    str: KhAlg<R>,
    deg_shift: (isize, isize),
    reduced: bool,
    canon_cycles: Vec<KhChain<R>>,
    gen_grid: OnceLock<GrMod2<KhState, R>>,
}

impl<R> KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        let c = KhComplex::new(l, h, t, reduced); 
        Self::from(&c)
    }
    
    pub fn new_no_simplify(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        let c = KhComplex::new_no_simplify(l, h, t, reduced); 
        Self::from(&c)
    }
    
    pub(crate) fn new_impl(inner: GrMod1<KhState, R>, str: KhAlg<R>, deg_shift: (isize, isize), reduced: bool, canon_cycles: Vec<KhChain<R>>) -> Self {
        Self { inner, str, deg_shift, reduced, canon_cycles, gen_grid: OnceLock::new() }
    }

    pub fn inner(&self) -> &GrMod1<KhState, R> { 
        &self.inner
    }

    delegate! {
        to self.inner {
            pub fn support(&self) -> impl Iterator<Item = &isize> + '_;
            pub fn is_supported(&self, i: isize) -> bool;
        }
    }

    pub fn str(&self) -> &KhAlg<R> { 
        &self.str
    }

    pub fn deg_shift(&self) -> (isize, isize) { 
        self.deg_shift
    }

    pub fn is_reduced(&self) -> bool { 
        self.reduced
    }

    pub fn h_range(&self) -> RangeInclusive<isize> {
        self.support().filter(|&&i|
            !self[i].is_zero()
        ).copied().range().unwrap_or(0..=-1)
    }

    pub fn q_range(&self) -> RangeInclusive<isize> {
        self.support().flat_map(|&i|
            self[i].generators().map(|z| z.q_deg())
        ).range().unwrap_or(0..=-1)
    }

    pub fn delta_range(&self) -> RangeInclusive<isize> {
        self.support().flat_map(|&i|
            self[i].generators().map(|z| 2 * z.h_deg() - z.q_deg())
        ).range().unwrap_or(0..=-1)
    }

    pub fn canon_cycles(&self) -> &Vec<KhChain<R>> { 
        &self.canon_cycles
    }

    pub fn truncated(&self, range: RangeInclusive<isize>) -> Self {
        Self::new_impl(
            self.inner.truncated(range), 
            self.str.clone(), 
            self.deg_shift, 
            self.reduced, 
            self.canon_cycles.clone()
        )
    }

    fn gen_grid(&self) -> &GrMod2<KhState, R> {
        self.gen_grid.get_or_init(|| make_gen_grid(self.inner()))
    }
}

impl<R> From<&KhComplex<R>> for KhHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    fn from(c: &KhComplex<R>) -> Self {
        KhHomology::new_impl(
            c.inner().reduced().homology(), 
            c.str().clone(), 
            c.deg_shift(), 
            c.is_reduced(),
            c.canon_cycles().clone()
        )
    }
}

impl<R> Index<isize> for KhHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = Summand<KhState, R>;

    delegate! {
        to self.inner {
            fn index(&self, index: isize) -> &Self::Output;
        }
    }
}

impl<R> Index<(isize, isize)> for KhHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = Summand<KhState, R>;

    delegate! {
        to self.gen_grid() {
            fn index(&self, index: (isize, isize)) -> &Self::Output;
        }
    }
}

impl<R> ToSeqString<isize> for KhHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    delegate! {
        to self.inner { 
            fn label(&self) -> String;
            fn indices(&self) -> Vec<isize>;
            fn entry_at(&self, i: &isize) -> String;
        }
    }
}

impl<R> ToTableString<isize> for KhHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    fn labels(&self) -> (String, String) { 
        ("i".to_string(), "j".to_string())
    }

    fn indices(&self) -> (Vec<isize>, Vec<isize>) { 
        (self.h_range().collect(), self.q_range().step_by(2).collect())
    }

    fn entry_at(&self, i: &isize, j: &isize) -> String {
        if self[(*i, *j)].is_zero() {
            ".".to_string()
        } else {
            self[(*i, *j)].to_string()
        }
    }
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use yui_core::poly::Poly;
    use yui_core::num::FF2;
    
    use yui_link::Link;
    use super::*;
    
    #[test]
    fn kh_empty() {
        let l = Link::empty();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);
        assert_eq!(h.q_range(), 0..=0);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot_old();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);
        assert_eq!(h.q_range(), -1..=1);
        
        assert_eq!(h[0].rank(), 2);
        assert!(h[0].is_free());
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::test_data("3_1");
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), -3..=0);
        assert_eq!(h.q_range(), -9..=-1);

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
        let l = Link::test_data("3_1").mirror();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=3);
        assert_eq!(h.q_range(), 1..=9);

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
        let l = Link::test_data("4_1");
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), -2..=2);
        assert_eq!(h.q_range(), -5..=5);

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
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);
        assert_eq!(h.q_range(), 0..=0);

        let h = h.gen_grid();
        assert_eq!(h[(0,0)].rank(), 1);
        assert!(h[(0,0)].is_free());
    }

    #[test]
    fn kh_unknot_bigr() {
        let l = Link::unknot_old();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);
        assert_eq!(h.q_range(), -1..=1);

        let h = h.gen_grid();
        assert_eq!(h[(0,-1)].rank(), 1);
        assert!(h[(0,-1)].is_free());
        assert_eq!(h[(0, 1)].rank(), 1);
        assert!(h[(0, 1)].is_free());
    }

    #[test]
    fn kh_unknot_bigr_red() {
        let l = Link::unknot_old();
        let h = KhHomology::new(&l, &0, &0, true);

        let h = h.gen_grid();
        assert_eq!(h[(0, 0)].rank(), 1);
        assert!(h[(0, 0)].is_free());
    }

    #[test]
    fn kh_trefoil_bigr() {
        let l = Link::test_data("3_1");
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), -3..=0);
        assert_eq!(h.q_range(), -9..=-1);

        let h = h.gen_grid();
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
        let l = Link::test_data("3_1").mirror();
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=3);
        assert_eq!(h.q_range(), 1..=9);

        let h = h.gen_grid();
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
        let l = Link::test_data("3_1");
        let h = KhHomology::new(&l, &0, &0, true);

        assert_eq!(h.h_range(), -3..=0);
        assert_eq!(h.q_range(), -8..=-2);

        let h = h.gen_grid();
        assert_eq!(h[(-3,-8)].rank(), 1);
        assert!(h[(-3,-8)].is_free());
        assert_eq!(h[(-2,-6)].rank(), 1);
        assert!(h[(-2,-6)].is_free());
        assert_eq!(h[( 0,-2)].rank(), 1);
        assert!(h[( 0,-2)].is_free());
    }

    #[test]
    fn kh_figure8_bigr() {
        let l = Link::test_data("4_1");
        let h = KhHomology::new(&l, &0, &0, false);

        assert_eq!(h.h_range(), -2..=2);
        assert_eq!(h.q_range(), -5..=5);

        let h = h.gen_grid();
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
        let l = Link::test_data("4_1");
        let h = KhHomology::new(&l, &0, &0, true);

        assert_eq!(h.h_range(), -2..=2);
        assert_eq!(h.q_range(), -4..=4);

        let h = h.gen_grid();
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
    fn bn_trefoil() {
        type R = FF2;
        type P = Poly<'H', R>;

        let l = Link::test_data("3_1");
        let (h, t) = (P::variable(), P::zero());
        let kh = KhHomology::new(&l, &h, &t, false);

        assert_eq!(kh.h_range(), -2..=0);
        assert_eq!(kh.q_range(), -7..=-1);

        let kh = kh.gen_grid();
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

 #[cfg(test)]
mod tests_v1 {
    use num_traits::Zero;
    use yui_core::poly::Poly;
    use yui_core::num::FF2;
    use yui_link::Link;
    use super::*;
    
    #[test]
    fn kh_empty() {
        let l = Link::empty();
        let h = KhHomology::new_no_simplify(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let h = KhHomology::new_no_simplify(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);
        
        assert_eq!(h[0].rank(), 2);
        assert!(h[0].is_free());
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::test_data("3_1");
        let h = KhHomology::new_no_simplify(&l, &0, &0, false);

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
        let l = Link::test_data("3_1").mirror();
        let h = KhHomology::new_no_simplify(&l, &0, &0, false);

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
        let l = Link::test_data("4_1");
        let h = KhHomology::new_no_simplify(&l, &0, &0, false);

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
        let h = KhHomology::new_no_simplify(&l, &0, &0, false);

        assert_eq!(h.h_range(), 0..=0);
        assert_eq!(h.q_range(), 0..=0);
        
        assert_eq!(h[(0,0)].rank(), 1);
        assert!(h[(0,0)].is_free());
    }

    #[test]
    fn kh_unknot_bigr() {
        let l = Link::unknot();
        let h = KhHomology::new_no_simplify(&l, &0, &0, false);

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
        let h = KhHomology::new_no_simplify(&l, &0, &0, true);

        assert_eq!(h.h_range(), 0..=0);
        assert_eq!(h.q_range(), 0..=0);

        assert_eq!(h[(0, 0)].rank(), 1);
        assert!(h[(0, 0)].is_free());
    }

    #[test]
    fn kh_trefoil_bigr() {
        let l = Link::test_data("3_1");
        let h = KhHomology::new_no_simplify(&l, &0, &0, false);

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
        let l = Link::test_data("3_1").mirror();
        let h = KhHomology::new_no_simplify(&l, &0, &0, false);

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
        let l = Link::test_data("3_1");
        let h = KhHomology::new_no_simplify(&l, &0, &0, true);

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
        let l = Link::test_data("4_1");
        let h = KhHomology::new_no_simplify(&l, &0, &0, false);

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
        let l = Link::test_data("4_1");
        let h = KhHomology::new_no_simplify(&l, &0, &0, true);

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
    fn trefoil_bn() {
        type R = FF2;
        type P = Poly<'H', R>;

        let l = Link::test_data("3_1");
        let (h, t) = (P::variable(), P::zero());
        let kh = KhHomology::new_no_simplify(&l, &h, &t, false);

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