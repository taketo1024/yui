use std::ops::{RangeInclusive, Index};
use std::sync::OnceLock;

use delegate::delegate;
use yui_core::lc::Lc;
use yui_core::{IteratorExt, Ring, RingOps, EucRing, EucRingOps};
use yui_link::Link;
use yui_homology::{ChainComplex1, ToSeqString, ToTableString, GrMod1, GrMod2, Summand};

use crate::kh::{KhGen, KhHomology};
use crate::tng::builder::BuildConfig;
use crate::util::Bigraded;

use super::KhAlg;
use yui_core::TeX;
use yui_homology::tex::{ToTexSeq, ToTexTable};

pub type KhChain<R> = Lc<KhGen, R>;
pub type KhComplexSummand<R> = Summand<KhGen, R>;

#[derive(Clone)]
pub struct KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    inner: ChainComplex1<KhGen, R>,
    alg: KhAlg<R>,
    deg_shift: (isize, isize),
    reduced: bool,
    canon_cycles: Vec<KhChain<R>>,
    cache_bigr: OnceLock<GrMod2<KhGen, R>>,
}

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        Self::new_partial(l, h, t, reduced, None)
    }

    // restricts the build to `h_range` (literal truncation; `None` = full).
    pub fn new_partial(l: &Link, h: &R, t: &R, reduced: bool, h_range: Option<RangeInclusive<isize>>) -> Self {
        Self::new_with_config(l, h, t, reduced, BuildConfig { h_range, ..Default::default() })
    }

    pub fn new_with_config(l: &Link, h: &R, t: &R, reduced: bool, config: BuildConfig) -> Self {
        use crate::tng::builder::TngComplexBuilder;

        assert!(!reduced || (!l.is_empty() && t.is_zero()));

        let config = BuildConfig {
            h_range: config.h_range.map(|r| Self::clamp_h_range(l, reduced, r)),
            ..config
        };
        let b = TngComplexBuilder::from_link(l, h, t, reduced).with_config(config).run();
        let canon_cycles = b.eval_elements();
        let inner = b.into_raw_complex(); // applies config.q_range on the no_full_deloop path

        KhComplex::from_raw_complex(l, h, t, reduced, inner, canon_cycles)
    }

    pub fn new_no_simplify(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        use super::cube::KhCube;

        assert!(!reduced || (!l.is_empty() && t.is_zero()));

        let base_pt = if reduced { l.base_pt() } else { None };
        let deg_shift = Self::deg_shift_for(l, reduced);
        let cube = KhCube::new(l, h, t, base_pt, deg_shift);
        let inner = cube.into_complex();

        let canon_cycles = if t.is_zero() && l.is_knot() {
            Self::make_canon_cycles(l, &R::zero(), h, reduced)
        } else {
            vec![]
        };

        KhComplex::from_raw_complex(l, h, t, reduced, inner, canon_cycles)
    }

    pub(crate) fn from_raw_complex(l: &Link, h: &R, t: &R, reduced: bool, inner: ChainComplex1<KhGen, R>, canon_cycles: Vec<KhChain<R>>) -> Self {
        let alg = KhAlg::new(h, t);
        let deg_shift = Self::deg_shift_for(l, reduced);

        KhComplex { inner, alg, deg_shift, reduced, canon_cycles, cache_bigr: OnceLock::new() }
    }

    pub fn deg_shift_for(l: &Link, reduced: bool) -> (isize, isize) {
        let (n_pos, n_neg) = l.n_signed_crossings();
        let (n_pos, n_neg) = (n_pos as isize, n_neg as isize);
        let h = -n_neg;
        let q = n_pos - 2 * n_neg;
        let e = if reduced { 1 } else { 0 };
        (h, q + e)
    }

    /// Clamp an h_range to the degree span, so open-ended ranges (`..=b` / `a..=`) don't overflow
    /// the build's `(a-1)..=(b+1)` widening. Span = `[deg_shift.0, deg_shift.0 + n_crossings + 1]`.
    pub fn clamp_h_range(l: &Link, reduced: bool, range: RangeInclusive<isize>) -> RangeInclusive<isize> {
        let lo = Self::deg_shift_for(l, reduced).0;
        let hi = lo + l.n_crossings() as isize + 1;
        (*range.start()).max(lo) ..= (*range.end()).min(hi)
    }

    pub fn inner(&self) -> &ChainComplex1<KhGen, R> {
        &self.inner
    }

    delegate! {
        to self.inner {
            pub fn support(&self) -> impl Iterator<Item = &isize> + '_;
            pub fn is_supported(&self, i: isize) -> bool;
            pub fn d_deg(&self) -> isize;
            pub fn d(&self, i: isize, z: &KhChain<R>) -> KhChain<R>;
            pub fn describe_d(&self) -> String;
            pub fn describe_d_at(&self, i: isize) -> String;
        }
    }

    pub fn alg(&self) -> &KhAlg<R> {
        &self.alg
    }

    pub fn deg_shift(&self) -> (isize, isize) { 
        self.deg_shift
    }

    pub fn is_reduced(&self) -> bool {
        self.reduced
    }

    pub fn h_deg_of(&self, x: &KhGen) -> isize {
        self.deg_shift.0 + x.rel_h_deg()
    }

    pub fn q_deg_of(&self, x: &KhGen) -> isize {
        self.deg_shift.1 + x.rel_q_deg()
    }

    pub fn h_deg_of_chain(&self, z: &KhChain<R>) -> isize {
        z.keys().map(|x| self.h_deg_of(x)).min().unwrap_or(0)
    }

    pub fn q_deg_of_chain(&self, z: &KhChain<R>) -> isize {
        z.keys().map(|x| self.q_deg_of(x)).min().unwrap_or(0)
    }

    pub fn h_range(&self) -> RangeInclusive<isize> {
        self.support().copied().range().unwrap_or(0..=-1)
    }

    pub fn q_range(&self) -> RangeInclusive<isize> {
        self.support().flat_map(|&i|
            self[i].raw_generators().iter().map(|x| self.q_deg_of(x))
        ).range().unwrap_or(0..=-1)
    }

    pub fn canon_cycles(&self) -> &Vec<KhChain<R>> { 
        &self.canon_cycles
    }

    fn cached_bigraded(&self) -> &GrMod2<KhGen, R> {
        self.cache_bigr.get_or_init(|| self.bigraded())
    }
}

impl<R> KhComplex<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology(&self) -> KhHomology<R> {
        self.into()
    }
}

impl<R> Bigraded<KhGen, R> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn base(&self) -> &GrMod1<KhGen, R> { self.inner.summands() }
    fn decomp_key(&self, z: &KhChain<R>) -> isize { self.q_deg_of_chain(z) }
}

impl<R> Index<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = KhComplexSummand<R>;

    delegate! {
        to self.inner {
            fn index(&self, index: isize) -> &Self::Output;
        }
    }
}

impl<R> Index<(isize, isize)> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = KhComplexSummand<R>;

    fn index(&self, index: (isize, isize)) -> &Self::Output {
        &self.cached_bigraded()[index]
    }
}


impl<R> ToSeqString<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    delegate! {
        to self.inner { 
            fn label(&self) -> String;
            fn indices(&self) -> Vec<isize>;
            fn entry_at(&self, i: &isize) -> String;
        }
    }
}

impl<R> ToTexSeq<isize> for KhComplex<R>
where R: Ring + TeX, for<'x> &'x R: RingOps<R> {
    fn tex_entry_at(&self, i: &isize) -> String {
        if self[*i].is_zero() {
            ".".to_string()
        } else {
            self[*i].tex_string()
        }
    }
}

impl<R> ToTableString<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
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

impl<R> ToTexTable<isize> for KhComplex<R>
where R: Ring + TeX, for<'x> &'x R: RingOps<R> {
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
        use yui_link::Link;

    use super::KhComplex;

    #[test]
    fn ckh_trefoil() {
        let l = Link::test_data("3_1").mirror();
        let c = KhComplex::new(&l, &0, &0, false);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c.q_range(), -9..=-1);

        assert_eq!(c[-3].rank(), 2);
        assert_eq!(c[-2].rank(), 2);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);

        c.inner().check_d_all();
    }

    #[test]
    fn ckh_trefoil_red() {
        let l = Link::test_data("3_1").mirror();
        let c = KhComplex::new(&l, &0, &0, true);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c.q_range(), -8..=-2);

        assert_eq!(c[-3].rank(), 1);
        assert_eq!(c[-2].rank(), 1);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 1);

        c.inner().check_d_all();
    }

    #[test]
    fn ckh_trefoil_bigr() {
        let l = Link::test_data("3_1").mirror();
        let c = KhComplex::new(&l, &0, &0, false);

        assert_eq!(c[(-3, -9)].rank(), 1);
        assert_eq!(c[(-3, -7)].rank(), 1);
        assert_eq!(c[(-2, -7)].rank(), 1);
        assert_eq!(c[(-2, -5)].rank(), 1);
        assert_eq!(c[(0, -3)].rank(), 1);
        assert_eq!(c[(0, -1)].rank(), 1);
    }

    #[test]
    fn ckh_trefoil_bigr_red() {
        let l = Link::test_data("3_1").mirror();
        let c = KhComplex::new(&l, &0, &0, true);

        assert_eq!(c[(-3, -8)].rank(), 1);
        assert_eq!(c[(-2, -6)].rank(), 1);
        assert_eq!(c[(0, -2)].rank(), 1);
    }
}

#[cfg(test)]
mod tests_v1 {
        use yui_link::Link;

    use super::KhComplex;

    #[test]
    fn ckh_trefoil() {
        let l = Link::test_data("3_1").mirror();
        let c = KhComplex::new_no_simplify(&l, &0, &0, false);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c[-3].rank(), 8);
        assert_eq!(c[-2].rank(), 12);
        assert_eq!(c[-1].rank(), 6);
        assert_eq!(c[ 0].rank(), 4);    

        c.inner().check_d_all();
    }

    #[test]
    fn ckh_trefoil_red() {
        let l = Link::test_data("3_1").mirror();
        let c = KhComplex::new_no_simplify(&l, &0, &0, true);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c[-3].rank(), 4);
        assert_eq!(c[-2].rank(), 6);
        assert_eq!(c[-1].rank(), 3);
        assert_eq!(c[ 0].rank(), 2);

        c.inner().check_d_all();
    }
}
