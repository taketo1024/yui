use std::ops::{Index, RangeInclusive};
use std::sync::Arc;

use delegate::delegate;
use itertools::Itertools;
use num_traits::Zero;
use yui_core::{Ring, RingOps};
use yui_core::lc::{LcKey, Lc};
use yui_matrix::sparse::{SpMat, SpVec};

use crate::algo::ChainReducer;
use crate::{ToSeqString, ToTableString, GenericChainComplexBase, GrMod, AddInd, isize2, isize3};
use super::Summand;

#[cfg(feature = "multithread")]
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

pub type ChainComplex <X, R> = ChainComplexBase<isize,  X, R>;
pub type ChainComplex2<X, R> = ChainComplexBase<isize2, X, R>;
pub type ChainComplex3<X, R> = ChainComplexBase<isize3, X, R>;

#[derive(Clone)]
pub struct ChainComplexBase<I, X, R>
where 
    I: AddInd,
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    summands: GrMod<I, X, R>,
    d_deg: I,
    d_map: Arc<dyn Fn(I, &Lc<X, R>) -> Lc<X, R> + Send + Sync>,
}

impl<I, X, R> ChainComplexBase<I, X, R>
where
    I: AddInd,
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    pub fn new<F>(summands: GrMod<I, X, R>, d_deg: I, d_map: F) -> Self
    where F: Fn(I, &Lc<X, R>) -> Lc<X, R> + Send + Sync + 'static {
        assert!(summands.iter().all(|(_, s)| s.is_free()));

        let d_map = Arc::new(d_map);
        Self { summands, d_deg, d_map }
    }

    pub fn zero() -> Self {
        Self::new(GrMod::default(), I::zero(), |_, _| Lc::zero())
    }

    pub fn summands(&self) -> &GrMod<I, X, R> {
        &self.summands
    }

    delegate! {
        to self.summands {
            pub fn support(&self) -> impl Iterator<Item = &I> + '_;
            pub fn is_supported(&self, i: I) -> bool;
        }
    }

    pub fn d_deg(&self) -> I {
        self.d_deg
    }

    pub(crate) fn raw_d(&self) -> Arc<dyn Fn(I, &Lc<X, R>) -> Lc<X, R> + Send + Sync> {
        self.d_map.clone()
    }

    pub fn d(&self, i: I, z: &Lc<X, R>) -> Lc<X, R> {
        (self.d_map)(i, z)
    }

    pub fn d_matrix(&self, i: I) -> SpMat<R> {
        let m = self[i + self.d_deg].rank();
        let n = self[i].rank();

        cfg_if::cfg_if! {
            if #[cfg(feature = "multithread")] {
                let cols = (0..n).into_par_iter().map(|j|
                    self.d_matrix_col(i, j)
                ).collect::<Vec<_>>();
                SpMat::from_col_vecs(m, cols)
            } else {
                let cols = (0..n).map(|j|
                    self.d_matrix_col(i, j)
                );
                SpMat::from_col_vecs(m, cols)
            }
        }
    }

    #[inline(never)] // for profilability
    fn d_matrix_col(&self, i: I, j: usize) -> SpVec<R> { 
        let z = self[i].generator(j);
        let w = self.d(i, &z);
        self[i + self.d_deg].vectorize(&w)
    }

    pub fn describe_d_at(&self, i: I) -> String {
        let c0 = &self[i];
        let c1 = &self[i + self.d_deg];
        let d = self.d_matrix(i).into_dense();
        format!("d[{i}]: {c0} -> {c1}\n{d}")
    }

    pub fn describe_d(&self) -> String {
        self.support().filter_map(|&i|
            if self[i].rank() > 0 && self[i + self.d_deg].rank() > 0 && !self.d_matrix(i).is_zero() {
                Some(self.describe_d_at(i))
            } else {
                None
            }
        ).join("")
    }

    pub fn as_generic(&self) -> GenericChainComplexBase<I, R> {
        GenericChainComplexBase::from_d_matrices(
            self.d_deg,
            self.support().map(|&i| (i, self.d_matrix(i)))
        )
    }

    pub fn reduced(&self) -> ChainComplexBase<I, X, R> { 
        let r = ChainReducer::reduce(self, true);

        let summands = GrMod::generate(
            self.summands.support().copied(),
            |i| {
                let c = &self[i];
                Summand::new(
                    c.raw_generators().clone(), 
                    r.rank(i).unwrap(), 
                    vec![], 
                    c.trans().merged(r.trans(i).unwrap())
                )
            } 
        );

        let d_deg = self.d_deg;
        let d_map = self.d_map.clone();
        Self { summands, d_deg, d_map }
    }

    pub fn reduced_generic(&self) -> GenericChainComplexBase<I, R> { 
        let r = ChainReducer::reduce(self, false);
        r.into_complex()
    }

    #[cfg(debug_assertions)]
    fn check_d_for(&self, i0: I, x: &X) {
        let i1 = i0 + self.d_deg();
        assert!(self.is_supported(i0), "Not supported: {i0}.");

        let dx = self.d(i0, &Lc::from(x.clone()));
        let ddx = self.d(i1, &dx);

        assert!(ddx.is_zero(), "d² is non-zero for {x} at {i0}.\n  dx: {dx}\n  ddx: {ddx}.");
    }

    #[cfg(debug_assertions)]
    pub fn check_d_at(&self, i0: I) {
        for x in self[i0].raw_generators().iter() {
            self.check_d_for(i0, x);
        }
    }

    #[cfg(debug_assertions)]
    pub fn check_d_all(&self) {
        for &i in self.support() {
            self.check_d_at(i);
        }
    }
}

impl<X, R> ChainComplex<X, R>
where 
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    pub fn truncated(&self, range: RangeInclusive<isize>) -> Self { 
        let d_deg = self.d_deg;
        let d_map = self.d_map.clone();
        let summands = self.summands.truncated(range.clone());

        Self::new(summands, d_deg, move |i, z| 
            if range.contains(&(i + d_deg)) { 
                d_map(i, z)
            } else { 
                Lc::zero()
            }
        )
    }
}

impl<I, X, R> Index<I> for ChainComplexBase<I, X, R>
where I: AddInd, X: LcKey, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = Summand<X, R>;
    fn index(&self, i: I) -> &Self::Output {
        &self.summands[i]
    }
}

impl<X, R> Index<(isize, isize)> for ChainComplex2<X, R>
where X: LcKey, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = Summand<X, R>;
    fn index(&self, i: (isize, isize)) -> &Self::Output {
        &self[isize2::from(i)]
    }
}

impl<X, R> Index<(isize, isize, isize)> for ChainComplex3<X, R>
where X: LcKey, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = Summand<X, R>;
    fn index(&self, i: (isize, isize, isize)) -> &Self::Output {
        &self[isize3::from(i)]
    }
}

impl<X, R> ToSeqString<isize> for ChainComplex<X, R>
where X: LcKey, R: Ring, for<'x> &'x R: RingOps<R> {
    delegate! {
        to self.summands { 
            fn label(&self) -> String;
            fn indices(&self) -> Vec<isize>;
            fn entry_at(&self, i: &isize) -> String;
        }
    }
}

impl<X, R> ToTableString<isize> for ChainComplex2<X, R>
where X: LcKey, R: Ring, for<'x> &'x R: RingOps<R> {
    delegate! {
        to self.summands { 
            fn labels(&self) -> (String, String);
            fn indices(&self) -> (Vec<isize>, Vec<isize>);
            fn entry_at(&self, i: &isize, j: &isize) -> String;
        }
    }
}