use std::ops::{Index, RangeInclusive};
use delegate::delegate;
use yui_core::{EucRing, EucRingOps};
use yui_homology::{Grid2, GridTrait, Homology, Summand, SummandTrait};
use yui_link::InvLink;
use crate::kh::KhChainExt;
use crate::khi::{KhIComplex, KhIGen};
use crate::misc::{make_gen_grid, range_of};

use super::KhIChain;

#[derive(Clone)]
pub struct KhIHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    inner: Homology<KhIGen, R>,
    canon_cycles: Vec<KhIChain<R>>
}

impl<R> KhIHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(l: &InvLink, h: &R, t: &R, reduced: bool) -> Self {
        let c = KhIComplex::new(l, h, t, reduced);
        Self::from(&c)
    }

    pub fn new_no_simplify(l: &InvLink, h: &R, t: &R, reduced: bool) -> Self {
        let c = KhIComplex::new_no_simplify(l, h, t, reduced);
        Self::from(&c)
    }

    pub(crate) fn new_impl(inner: Homology<KhIGen, R>, canon_cycles: Vec<KhIChain<R>>) -> Self {
        Self { inner, canon_cycles }
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        range_of(self.support().filter(|&i| 
            !self[i].is_zero()
        ))
    }

    pub fn q_range(&self) -> RangeInclusive<isize> {
        range_of(self.support().flat_map(|i| 
            self[i].gens().map(|z| z.q_deg())
        ))
    }
    
    pub fn canon_cycles(&self) -> &[KhIChain<R>] { 
        &self.canon_cycles
    }

    pub fn inner(&self) -> &Homology<KhIGen, R> { 
        &self.inner
    }

    pub fn truncated(&self, range: RangeInclusive<isize>) -> Self {
        Self::new_impl(
            self.inner.truncated(range),
            self.canon_cycles.clone()
        )
    }

    pub fn gen_grid(self) -> Grid2<Summand<KhIGen, R>> { 
        make_gen_grid(self.inner())
    }
}

impl<R> From<&KhIComplex<R>> for KhIHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    fn from(c: &KhIComplex<R>) -> Self {
        KhIHomology::new_impl(
            c.inner().reduced().homology(), 
            c.canon_cycles().iter().cloned().collect()
        )
    }
}

impl<R> GridTrait<isize> for KhIHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Support = std::vec::IntoIter<isize>;
    type Item = Summand<KhIGen, R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Support;
            fn is_supported(&self, i: isize) -> bool;
            fn get(&self, i: isize) -> &Self::Item;
            fn get_default(&self) -> &Self::Item;
        }
    }
}

impl<R> Index<isize> for KhIHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = Summand<KhIGen, R>;

    delegate! { 
        to self.inner { 
            fn index(&self, index: isize) -> &Self::Output;
        }
    }
}

mod tests {
    #![allow(unused)]

    use itertools::Itertools;
    use yui_core::poly::HPoly;
    use yui_core::num::FF2;
    use num_traits::{Zero, One};
    use yui_homology::{ChainComplexTrait, DisplaySeq, DisplayTable, SummandTrait};
    use yui_link::Link;
    use super::*;

    #[test]
    fn khi() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        let (h, t) = (R::zero(), R::zero());
        let khi = KhIHomology::new(&l, &h, &t, false);

        assert_eq!(khi.h_range(), 0..=4);
        assert_eq!(khi[0].rank(), 2);
        assert_eq!(khi[1].rank(), 2);
        assert_eq!(khi[2].rank(), 2);
        assert_eq!(khi[3].rank(), 4);
        assert_eq!(khi[4].rank(), 2);
    }

    #[test]
    fn khi_fbn() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        let (h, t) = (R::one(), R::zero());
        let khi = KhIHomology::new(&l, &h, &t, false);

        assert_eq!(khi.h_range(), 0..=1);
        assert_eq!(khi[0].rank(), 2);
        assert_eq!(khi[1].rank(), 2);
        assert_eq!(khi[2].rank(), 0);
        assert_eq!(khi[3].rank(), 0);
        assert_eq!(khi[4].rank(), 0);
    }

    #[test]
    fn khi_bn() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        type P = HPoly<'H', R>;
        let (h, t) = (P::variable(), P::zero());

        let khi = KhIHomology::new(&l, &h, &t, false);

        assert_eq!(khi.h_range(), 0..=4);
        assert_eq!(khi[0].rank(), 2);
        assert_eq!(khi[1].rank(), 2);
        assert_eq!(khi[2].rank(), 0);
        assert_eq!(khi[3].rank(), 0);
        assert_eq!(khi[3].tors(), [P::variable(), P::variable()]);
        assert_eq!(khi[4].rank(), 0);
        assert_eq!(khi[4].tors(), [P::variable(), P::variable()]);
    }

    #[test]
    fn khi_red() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        let (h, t) = (R::zero(), R::zero());
        let khi = KhIHomology::new(&l, &h, &t, true);

        assert_eq!(khi.h_range(), 0..=4);
        assert_eq!(khi[0].rank(), 1);
        assert_eq!(khi[1].rank(), 1);
        assert_eq!(khi[2].rank(), 1);
        assert_eq!(khi[3].rank(), 2);
        assert_eq!(khi[4].rank(), 1);
    }

    #[test]
    fn khi_fbn_red() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        let (h, t) = (R::one(), R::zero());
        let khi = KhIHomology::new(&l, &h, &t, true);

        assert_eq!(khi.h_range(), 0..=1);
        assert_eq!(khi[0].rank(), 1);
        assert_eq!(khi[1].rank(), 1);
        assert_eq!(khi[2].rank(), 0);
        assert_eq!(khi[3].rank(), 0);
        assert_eq!(khi[4].rank(), 0);
    }

    #[test]
    fn khi_bn_red() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        type P = HPoly<'H', R>;
        let (h, t) = (P::variable(), P::zero());

        let khi = KhIHomology::new(&l, &h, &t, true);

        assert_eq!(khi.h_range(), 0..=4);
        assert_eq!(khi[0].rank(), 1);
        assert_eq!(khi[1].rank(), 1);
        assert_eq!(khi[2].rank(), 0);
        assert_eq!(khi[3].rank(), 0);
        assert_eq!(khi[3].tors(), [P::variable()]);
        assert_eq!(khi[4].rank(), 0);
        assert_eq!(khi[4].tors(), [P::variable()]);
    }

    #[test]
    fn khi_kh_bigr() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        let (h, t) = (R::zero(), R::zero());
        let khi = KhIHomology::new(&l, &h, &t, false).gen_grid();

        assert_eq!(khi[(0, 1)].rank(), 1);
        assert_eq!(khi[(0, 3)].rank(), 1);
        assert_eq!(khi[(1, 1)].rank(), 1);
        assert_eq!(khi[(1, 3)].rank(), 1);
        assert_eq!(khi[(2, 5)].rank(), 1);
        assert_eq!(khi[(2, 7)].rank(), 1);
        assert_eq!(khi[(3, 5)].rank(), 1);
        assert_eq!(khi[(3, 7)].rank(), 2);
        assert_eq!(khi[(3, 9)].rank(), 1);
        assert_eq!(khi[(4, 7)].rank(), 1);
        assert_eq!(khi[(4, 9)].rank(), 1);
    }

    #[test]
    fn khi_kh_red_bigr() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        let (h, t) = (R::zero(), R::zero());
        let khi = KhIHomology::new(&l, &h, &t, true).gen_grid();

        assert_eq!(khi[(0, 2)].rank(), 1);
        assert_eq!(khi[(1, 2)].rank(), 1);
        assert_eq!(khi[(2, 6)].rank(), 1);
        assert_eq!(khi[(3, 6)].rank(), 1);
        assert_eq!(khi[(3, 8)].rank(), 1);
        assert_eq!(khi[(4, 8)].rank(), 1);
    }
}

mod tests_v1 {
    #![allow(unused)]

    use itertools::Itertools;
    use yui_core::poly::HPoly;
    use yui_core::num::FF2;
    use num_traits::{Zero, One};
    use yui_homology::{ChainComplexTrait, DisplaySeq, DisplayTable, SummandTrait};
    use yui_link::Link;
    use super::*;

    #[test]
    fn khi() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        let (h, t) = (R::zero(), R::zero());
        let khi = KhIHomology::new_no_simplify(&l, &h, &t, false);

        assert_eq!(khi.h_range(), 0..=4);
        assert_eq!(khi[0].rank(), 2);
        assert_eq!(khi[1].rank(), 2);
        assert_eq!(khi[2].rank(), 2);
        assert_eq!(khi[3].rank(), 4);
        assert_eq!(khi[4].rank(), 2);
    }

    #[test]
    fn khi_fbn() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        let (h, t) = (R::one(), R::zero());
        let khi = KhIHomology::new_no_simplify(&l, &h, &t, false);

        assert_eq!(khi.h_range(), 0..=1);
        assert_eq!(khi[0].rank(), 2);
        assert_eq!(khi[1].rank(), 2);
        assert_eq!(khi[2].rank(), 0);
        assert_eq!(khi[3].rank(), 0);
        assert_eq!(khi[4].rank(), 0);
    }

    #[test]
    fn khi_bn() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        type P = HPoly<'H', R>;
        let (h, t) = (P::variable(), P::zero());

        let khi = KhIHomology::new_no_simplify(&l, &h, &t, false);

        assert_eq!(khi.h_range(), 0..=4);
        assert_eq!(khi[0].rank(), 2);
        assert_eq!(khi[1].rank(), 2);
        assert_eq!(khi[2].rank(), 0);
        assert_eq!(khi[3].rank(), 0);
        assert_eq!(khi[3].tors(), [P::variable(), P::variable()]);
        assert_eq!(khi[4].rank(), 0);
        assert_eq!(khi[4].tors(), [P::variable(), P::variable()]);
    }

    #[test]
    fn khi_red() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        let (h, t) = (R::zero(), R::zero());
        let khi = KhIHomology::new_no_simplify(&l, &h, &t, true);

        assert_eq!(khi.h_range(), 0..=4);
        assert_eq!(khi[0].rank(), 1);
        assert_eq!(khi[1].rank(), 1);
        assert_eq!(khi[2].rank(), 1);
        assert_eq!(khi[3].rank(), 2);
        assert_eq!(khi[4].rank(), 1);
    }

    #[test]
    fn khi_fbn_red() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        let (h, t) = (R::one(), R::zero());
        let khi = KhIHomology::new_no_simplify(&l, &h, &t, true);

        assert_eq!(khi.h_range(), 0..=1);
        assert_eq!(khi[0].rank(), 1);
        assert_eq!(khi[1].rank(), 1);
        assert_eq!(khi[2].rank(), 0);
        assert_eq!(khi[3].rank(), 0);
        assert_eq!(khi[4].rank(), 0);
    }

    #[test]
    fn khi_bn_red() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        type P = HPoly<'H', R>;
        let (h, t) = (P::variable(), P::zero());

        let khi = KhIHomology::new_no_simplify(&l, &h, &t, true);

        assert_eq!(khi.h_range(), 0..=4);
        assert_eq!(khi[0].rank(), 1);
        assert_eq!(khi[1].rank(), 1);
        assert_eq!(khi[2].rank(), 0);
        assert_eq!(khi[3].rank(), 0);
        assert_eq!(khi[3].tors(), [P::variable()]);
        assert_eq!(khi[4].rank(), 0);
        assert_eq!(khi[4].tors(), [P::variable()]);
    }

    #[test]
    fn khi_kh_bigr() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        let (h, t) = (R::zero(), R::zero());
        let khi = KhIHomology::new_no_simplify(&l, &h, &t, false).gen_grid();

        assert_eq!(khi[(0, 1)].rank(), 1);
        assert_eq!(khi[(0, 3)].rank(), 1);
        assert_eq!(khi[(1, 1)].rank(), 1);
        assert_eq!(khi[(1, 3)].rank(), 1);
        assert_eq!(khi[(2, 5)].rank(), 1);
        assert_eq!(khi[(2, 7)].rank(), 1);
        assert_eq!(khi[(3, 5)].rank(), 1);
        assert_eq!(khi[(3, 7)].rank(), 2);
        assert_eq!(khi[(3, 9)].rank(), 1);
        assert_eq!(khi[(4, 7)].rank(), 1);
        assert_eq!(khi[(4, 9)].rank(), 1);
    }

    #[test]
    fn khi_kh_red_bigr() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        type R = FF2;
        let (h, t) = (R::zero(), R::zero());
        let khi = KhIHomology::new_no_simplify(&l, &h, &t, true).gen_grid();

        assert_eq!(khi[(0, 2)].rank(), 1);
        assert_eq!(khi[(1, 2)].rank(), 1);
        assert_eq!(khi[(2, 6)].rank(), 1);
        assert_eq!(khi[(3, 6)].rank(), 1);
        assert_eq!(khi[(3, 8)].rank(), 1);
        assert_eq!(khi[(4, 8)].rank(), 1);
    }
}