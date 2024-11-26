use std::ops::{Index, RangeInclusive};
use cartesian::cartesian;
use delegate::delegate;
use yui::{EucRing, EucRingOps};
use yui_homology::{isize2, Grid2, GridTrait, Homology, Summand};
use yui_link::InvLink;
use crate::khi::{KhIComplex, KhIGen};
use crate::misc::{collect_gen_info, range_of};

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

    pub(crate) fn new_impl(inner: Homology<KhIGen, R>, canon_cycles: Vec<KhIChain<R>>) -> Self {
        Self { inner, canon_cycles }
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        range_of(self.support())
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

    pub fn into_bigraded(self) -> Grid2<Summand<KhIGen, R>> { 
        let table = collect_gen_info(self.inner());
        let h_range = range_of(table.keys().map(|i| i.0));
        let q_range = range_of(table.keys().map(|i| i.1)).step_by(2);
        let support = cartesian!(h_range, q_range.clone()).map(|(i, j)| 
            isize2(i, j)
        );

        Grid2::generate(support, move |idx| { 
            let i = idx.0;
            let Some(e) = table.get(&idx) else { 
                return Summand::zero()
            };
            
            let (rank, tors, indices) = e;
            let gens = self[i].raw_gens().clone(); 
            let trans = self[i].trans().sub(indices);
            Summand::new(gens, *rank, tors.clone(), trans)
        })
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
    use yui::poly::HPoly;
    use yui::FF2;
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

        assert_eq!(khi.h_range(), 0..=4);
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

        assert_eq!(khi.h_range(), 0..=4);
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

    // #[test]
    // fn khi_kh_bigr() { 
    //     let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

    //     type R = FF2;
    //     let (h, t) = (R::zero(), R::zero());
    //     let khi = KhIHomology::new(&l, &h, false).into_bigraded();

    //     assert_eq!(c[(0, 1)].rank(), 1);
    //     assert_eq!(c[(0, 3)].rank(), 1);
    //     assert_eq!(c[(1, 1)].rank(), 1);
    //     assert_eq!(c[(1, 3)].rank(), 1);
    //     assert_eq!(c[(2, 5)].rank(), 1);
    //     assert_eq!(c[(2, 7)].rank(), 1);
    //     assert_eq!(c[(3, 5)].rank(), 1);
    //     assert_eq!(c[(3, 7)].rank(), 2);
    //     assert_eq!(c[(3, 9)].rank(), 1);
    //     assert_eq!(c[(4, 7)].rank(), 1);
    //     assert_eq!(c[(4, 9)].rank(), 1);

    //     c.check_d_all();
    // }

    // #[test]
    // fn khi_kh_red_bigr() { 
    //     let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

    //     type R = FF2;
    //     let (h, t) = (R::zero(), R::zero());
    //     let khi = KhIHomology::new(&l, &h, true, false).into_bigraded();

    //     assert_eq!(c[(0, 2)].rank(), 1);
    //     assert_eq!(c[(1, 2)].rank(), 1);
    //     assert_eq!(c[(2, 6)].rank(), 1);
    //     assert_eq!(c[(3, 6)].rank(), 1);
    //     assert_eq!(c[(3, 8)].rank(), 1);
    //     assert_eq!(c[(4, 8)].rank(), 1);

    //     c.check_d_all();
    // }
}