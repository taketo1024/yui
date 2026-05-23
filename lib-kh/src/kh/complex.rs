use std::ops::{RangeInclusive, Index};
use std::sync::OnceLock;

use delegate::delegate;
use yui_core::lc::Lc;
use yui_core::{IteratorExt, Ring, RingOps, EucRing, EucRingOps};
use yui_link::Link;
use yui_homology::{ChainComplex1, ToSeqString, ToTableString, GrMod2, Summand};

use crate::kh::{KhGen, KhHomology};
use crate::misc::decomp_by_q_deg;

use super::KhAlg;

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
    gen_grid: OnceLock<GrMod2<KhGen, R>>,
}

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        use crate::kh::internal::v2::builder::TngComplexBuilder;

        assert!(!reduced || (!l.is_empty() && t.is_zero()));

        TngComplexBuilder::build_kh_complex(l, h, t, reduced)
    }

    pub fn new_no_simplify(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        use crate::kh::internal::v1::cube::KhCube;

        assert!(!reduced || (!l.is_empty() && t.is_zero()));

        let base_pt = if reduced { l.base_pt() } else { None };
        let deg_shift = Self::deg_shift_for(l, reduced);

        let alg = KhAlg::new(h, t);
        let cube = KhCube::new(l, h, t, base_pt, deg_shift);
        let complex = cube.into_complex();

        let canon_cycles = if t.is_zero() && l.is_knot() {
            Self::make_canon_cycles(l, &R::zero(), h, reduced)
        } else {
            vec![]
        };

        KhComplex::new_impl(complex, alg, deg_shift, reduced, canon_cycles)
    }

    pub(crate) fn new_impl(inner: ChainComplex1<KhGen, R>, alg: KhAlg<R>, deg_shift: (isize, isize), reduced: bool, canon_cycles: Vec<KhChain<R>>) -> Self {
        KhComplex { inner, alg, deg_shift, reduced, canon_cycles, gen_grid: OnceLock::new() }
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

    pub fn inner(&self) -> &ChainComplex1<KhGen, R> {
        &self.inner
    }

    fn gen_grid(&self) -> &GrMod2<KhGen, R> {
        self.gen_grid.get_or_init(|| 
            decomp_by_q_deg(self.inner.summands(), |z| self.q_deg_of_chain(z))
        )
    }

    pub fn deg_shift_for(l: &Link, reduced: bool) -> (isize, isize) {
        let (n_pos, n_neg) = l.n_signed_crossings();
        let (n_pos, n_neg) = (n_pos as isize, n_neg as isize);
        let h = -n_neg;
        let q = n_pos - 2 * n_neg;
        let e = if reduced { 1 } else { 0 };
        (h, q + e)
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
}

impl<R> KhComplex<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology(&self) -> KhHomology<R> {
        self.into()
    }
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

    delegate! { 
        to self.gen_grid() {
            fn index(&self, index: (isize, isize)) -> &Self::Output;
        }
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

#[cfg(test)]
mod tests {
        use yui_link::Link;

    use super::KhComplex;

    #[test]
    fn ckh_trefoil() {
        let l = Link::test_data("3_1");
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
        let l = Link::test_data("3_1");
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
        let l = Link::test_data("3_1");
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
        let l = Link::test_data("3_1");
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
        let l = Link::test_data("3_1");
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
        let l = Link::test_data("3_1");
        let c = KhComplex::new_no_simplify(&l, &0, &0, true);

        assert_eq!(c.h_range(), -3..=0);
        assert_eq!(c[-3].rank(), 4);
        assert_eq!(c[-2].rank(), 6);
        assert_eq!(c[-1].rank(), 3);
        assert_eq!(c[ 0].rank(), 2);

        c.inner().check_d_all();
    }
}
