use std::ops::{Index, RangeInclusive};
use std::sync::OnceLock;
use delegate::delegate;
use yui_core::abst::{EucRing, EucRingOps};
use yui_core::ext::IteratorExt;
use yui_homology::{ToSeqString, ToTableString, GrMod1, GrMod2, Summand};
use yui_link::InvLink;
use crate::kh::KhComplex;
use crate::khi::{KhIComplex, KhIGen, KhIGenExt};
use crate::tng::builder::SymBuildConfig;
use crate::util::Bigraded;

use super::KhIChain;
use yui_core::util::tex::TeX;
use yui_homology::tex::{ToTexSeq, ToTexTable};

#[derive(Clone)]
pub struct KhIHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    inner: GrMod1<KhIGen, R>,
    canon_cycles: Vec<KhIChain<R>>,
    deg_shift: (isize, isize),
    cache_bigr: OnceLock<GrMod2<KhIGen, R>>,
}

impl<R> KhIHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(l: &InvLink, h: &R, t: &R, reduced: bool) -> Self {
        let c = KhIComplex::new(l, h, t, reduced);
        Self::from(&c)
    }

    // builds one degree wider (for boundary maps), then truncates to `h_range`.
    pub fn new_partial(l: &InvLink, h: &R, t: &R, reduced: bool, h_range: Option<RangeInclusive<isize>>) -> Self {
        Self::new_with_config(l, h, t, reduced, SymBuildConfig { h_range, ..Default::default() })
    }

    // builds one degree wider (for boundary maps), then truncates to `config.h_range`;
    // the rest of `config` (e.g. `chunk_bound`) flows down to the complex build.
    pub fn new_with_config(l: &InvLink, h: &R, t: &R, reduced: bool, config: SymBuildConfig) -> Self {
        let Some(range) = config.h_range.clone() else {
            return Self::from(&KhIComplex::new_with_config(l, h, t, reduced, config));
        };
        let range = KhComplex::<R>::clamp_h_range(l.inner(), reduced, range); // resolve open ends before truncating
        let (a, b) = (*range.start(), *range.end());
        let cone_config = SymBuildConfig { h_range: Some((a - 1)..=(b + 1)), ..config };
        let c = KhIComplex::new_with_config(l, h, t, reduced, cone_config);
        Self::from_complex(&c, Some(a..=b))
    }

    fn from_complex(c: &KhIComplex<R>, range: Option<RangeInclusive<isize>>) -> Self {
        let reduced = c.inner().reduced();
        let homology = match &range {
            Some(r) => reduced.homology_in(r.clone()),
            None    => reduced.homology(),
        };
        // drop canon cycles whose h-degree falls outside the requested range (e.g. the Q-side at h+1).
        let canon_cycles = c.canon_cycles().iter()
            .filter(|z| range.as_ref().map_or(true, |r| r.contains(&c.h_deg_of_chain(z))))
            .cloned()
            .collect();
        KhIHomology::new_impl(homology, canon_cycles, c.deg_shift())
    }

    pub fn new_no_simplify(l: &InvLink, h: &R, t: &R, reduced: bool) -> Self {
        let c = KhIComplex::new_no_simplify(l, h, t, reduced);
        Self::from(&c)
    }

    pub(crate) fn new_impl(inner: GrMod1<KhIGen, R>, canon_cycles: Vec<KhIChain<R>>, deg_shift: (isize, isize)) -> Self {
        Self { inner, canon_cycles, deg_shift, cache_bigr: OnceLock::new() }
    }

    pub fn inner(&self) -> &GrMod1<KhIGen, R> {
        &self.inner
    }

    delegate! {
        to self.inner {
            pub fn support(&self) -> impl Iterator<Item = &isize> + '_;
            pub fn is_supported(&self, i: isize) -> bool;
        }
    }

    pub fn deg_shift(&self) -> (isize, isize) {
        self.deg_shift
    }

    pub fn h_deg_of(&self, x: &KhIGen) -> isize {
        self.deg_shift.0 + x.rel_h_deg()
    }

    pub fn q_deg_of(&self, x: &KhIGen) -> isize {
        self.deg_shift.1 + x.rel_q_deg()
    }

    pub fn h_deg_of_chain(&self, z: &KhIChain<R>) -> isize {
        z.keys().map(|x| self.h_deg_of(x)).min().unwrap_or(0)
    }

    pub fn q_deg_of_chain(&self, z: &KhIChain<R>) -> isize {
        z.keys().map(|x| self.q_deg_of(x)).min().unwrap_or(0)
    }

    pub fn h_range(&self) -> RangeInclusive<isize> {
        self.support().filter(|&&i|
            !self[i].is_zero()
        ).copied().range().unwrap_or(0..=-1)
    }

    pub fn q_range(&self) -> RangeInclusive<isize> {
        self.support().flat_map(|&i|
            self[i].generators().map(|z| self.q_deg_of_chain(&z))
        ).range().unwrap_or(0..=-1)
    }

    pub fn canon_cycles(&self) -> &[KhIChain<R>] {
        &self.canon_cycles
    }

    pub fn truncated(&self, range: RangeInclusive<isize>) -> Self {
        Self::new_impl(
            self.inner.truncated(range),
            self.canon_cycles.clone(),
            self.deg_shift,
        )
    }

    fn cached_bigraded(&self) -> &GrMod2<KhIGen, R> {
        self.cache_bigr.get_or_init(|| self.bigraded())
    }
}

impl<R> Bigraded<KhIGen, R> for KhIHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    fn base(&self) -> &GrMod1<KhIGen, R> { self.inner() }
    fn decomp_key(&self, z: &KhIChain<R>) -> isize { self.q_deg_of_chain(z) }
}

impl<R> From<&KhIComplex<R>> for KhIHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    fn from(c: &KhIComplex<R>) -> Self {
        KhIHomology::from_complex(c, None)
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

impl<R> Index<(isize, isize)> for KhIHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = Summand<KhIGen, R>;

    fn index(&self, index: (isize, isize)) -> &Self::Output {
        &self.cached_bigraded()[index]
    }
}

impl<R> ToSeqString<isize> for KhIHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    delegate! {
        to self.inner { 
            fn label(&self) -> String;
            fn indices(&self) -> Vec<isize>;
            fn entry_at(&self, i: &isize) -> String;
        }
    }
}

impl<R> ToTexSeq<isize> for KhIHomology<R>
where R: EucRing + TeX, for<'x> &'x R: EucRingOps<R> {
    fn tex_entry_at(&self, i: &isize) -> String {
        if self[*i].is_zero() {
            ".".to_string()
        } else {
            self[*i].tex_string()
        }
    }
}

impl<R> ToTableString<isize> for KhIHomology<R>
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

impl<R> ToTexTable<isize> for KhIHomology<R>
where R: EucRing + TeX, for<'x> &'x R: EucRingOps<R> {
    fn tex_entry_at(&self, i: &isize, j: &isize) -> String {
        if self[(*i, *j)].is_zero() {
            ".".to_string()
        } else {
            self[(*i, *j)].tex_string()
        }
    }
}

#[cfg(test)]
mod tests {
    use num_traits::{Zero, One};
    use yui_core::num::FF2;
    use yui_core::poly::Poly;
    use super::*;

    // the same assertions for both builders: `new` (simplified) and `new_no_simplify`.
    macro_rules! khi_homology_tests {
        ($build:expr) => {
            #[test]
            fn khi() {
                let l = InvLink::test_data("3_1");

                type R = FF2;
                let (h, t) = (R::zero(), R::zero());
                let khi = $build(&l, &h, &t, false);

                assert_eq!(khi.h_range(), 0..=4);
                assert_eq!(khi[0].rank(), 2);
                assert_eq!(khi[1].rank(), 2);
                assert_eq!(khi[2].rank(), 2);
                assert_eq!(khi[3].rank(), 4);
                assert_eq!(khi[4].rank(), 2);
            }

            #[test]
            fn khi_fbn() {
                let l = InvLink::test_data("3_1");

                type R = FF2;
                let (h, t) = (R::one(), R::zero());
                let khi = $build(&l, &h, &t, false);

                assert_eq!(khi.h_range(), 0..=1);
                assert_eq!(khi[0].rank(), 2);
                assert_eq!(khi[1].rank(), 2);
                assert_eq!(khi[2].rank(), 0);
                assert_eq!(khi[3].rank(), 0);
                assert_eq!(khi[4].rank(), 0);
            }

            #[test]
            fn khi_bn() {
                let l = InvLink::test_data("3_1");

                type R = FF2;
                type P = Poly<'H', R>;
                let (h, t) = (P::variable(), P::zero());
                let khi = $build(&l, &h, &t, false);

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
                let l = InvLink::test_data("3_1");

                type R = FF2;
                let (h, t) = (R::zero(), R::zero());
                let khi = $build(&l, &h, &t, true);

                assert_eq!(khi.h_range(), 0..=4);
                assert_eq!(khi[0].rank(), 1);
                assert_eq!(khi[1].rank(), 1);
                assert_eq!(khi[2].rank(), 1);
                assert_eq!(khi[3].rank(), 2);
                assert_eq!(khi[4].rank(), 1);
            }

            #[test]
            fn khi_fbn_red() {
                let l = InvLink::test_data("3_1");

                type R = FF2;
                let (h, t) = (R::one(), R::zero());
                let khi = $build(&l, &h, &t, true);

                assert_eq!(khi.h_range(), 0..=1);
                assert_eq!(khi[0].rank(), 1);
                assert_eq!(khi[1].rank(), 1);
                assert_eq!(khi[2].rank(), 0);
                assert_eq!(khi[3].rank(), 0);
                assert_eq!(khi[4].rank(), 0);
            }

            #[test]
            fn khi_bn_red() {
                let l = InvLink::test_data("3_1");

                type R = FF2;
                type P = Poly<'H', R>;
                let (h, t) = (P::variable(), P::zero());
                let khi = $build(&l, &h, &t, true);

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
                let l = InvLink::test_data("3_1");

                type R = FF2;
                let (h, t) = (R::zero(), R::zero());
                let khi = $build(&l, &h, &t, false);

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
                let l = InvLink::test_data("3_1");

                type R = FF2;
                let (h, t) = (R::zero(), R::zero());
                let khi = $build(&l, &h, &t, true);

                assert_eq!(khi[(0, 2)].rank(), 1);
                assert_eq!(khi[(1, 2)].rank(), 1);
                assert_eq!(khi[(2, 6)].rank(), 1);
                assert_eq!(khi[(3, 6)].rank(), 1);
                assert_eq!(khi[(3, 8)].rank(), 1);
                assert_eq!(khi[(4, 8)].rank(), 1);
            }
        };
    }

    // canon cycles outside the h-range are dropped: 3_1 has B-cycles at h0 and Q-cycles at h1,
    // so `..=0` keeps only the 2 B-cycles (was a crash when the Q-cycles indexed an unbuilt degree).
    #[test]
    fn canon_cycles_clipped_to_h_range() {
        let l = InvLink::test_data("3_1");
        type P = Poly<'H', FF2>;
        let (h, t) = (P::variable(), P::zero());

        assert_eq!(KhIHomology::new(&l, &h, &t, false).canon_cycles().len(), 4);
        let clipped = KhIHomology::new_partial(&l, &h, &t, false, Some(isize::MIN + 1 ..= 0));
        assert_eq!(clipped.canon_cycles().len(), 2);
    }

    #[test]
    fn khi_partial() {
        // partial KhI homology matches the full result across the whole window.
        let l = InvLink::test_data("3_1");

        type R = FF2;
        let (h, t) = (R::zero(), R::zero());
        let full = KhIHomology::new(&l, &h, &t, false);
        let part = KhIHomology::new_partial(&l, &h, &t, false, Some(1..=3));

        for i in 1..=3 {
            assert_eq!(part[i].rank(), full[i].rank(), "rank at {i}");
        }

        assert!(part[0].is_zero());
        assert!(part[4].is_zero());
    }

    mod v2 {
        use super::*;
        khi_homology_tests!(KhIHomology::new);
    }

    mod v1 {
        use super::*;
        khi_homology_tests!(KhIHomology::new_no_simplify);
    }
}
