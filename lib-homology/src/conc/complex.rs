//! [`ChainComplex<I, X, R>`]: the central chain-complex type, together with
//! homology computation via SNF.

use std::collections::HashMap;
use std::ops::{Index, RangeInclusive};
use std::sync::Arc;

use delegate::delegate;
use itertools::Itertools;
use num_traits::Zero;
use yui_core::{EucRing, EucRingOps, Ring, RingOps};
use yui_core::lc::{LcKey, Lc};

use yui_matrix::MatTrait;
use yui_matrix::sparse::{SpMat, SpVec};

use crate::algo::{ChainReducer, HomologyCalc};
use crate::{ToSeqString, ToTableString, GenericChainComplex, GenericGrMod, GenericSummand, GrMod, AddInd, isize2, isize3};
use super::Summand;

#[cfg(feature = "multithread")]
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

pub type ChainComplex1<X, R> = ChainComplex<isize,  X, R>;
pub type ChainComplex2<X, R> = ChainComplex<isize2, X, R>;
pub type ChainComplex3<X, R> = ChainComplex<isize3, X, R>;

/// An `I`-graded chain complex: a [`GrMod`] of [`Summand`]s, a degree shift
/// `d_deg`, and a differential closure `Fn(I, &Lc<X, R>) -> Lc<X, R>`. Optionally
/// caches per-index `SpMat<R>`s so [`Self::d_matrix`] becomes a clone.
#[derive(Clone)]
pub struct ChainComplex<I, X, R>
where
    I: AddInd,
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    summands: GrMod<I, X, R>,
    d_deg: I,
    d_map: Arc<dyn Fn(I, &Lc<X, R>) -> Lc<X, R> + Send + Sync>,
    d_matrices: Arc<HashMap<I, SpMat<R>>>,
}

impl<I, X, R> ChainComplex<I, X, R>
where
    I: AddInd,
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    pub fn new<F>(summands: GrMod<I, X, R>, d_deg: I, d_map: F) -> Self
    where F: Fn(I, &Lc<X, R>) -> Lc<X, R> + Send + Sync + 'static {
        assert!(summands.iter().all(|(_, s)| s.is_free()));

        let d_map = Arc::new(d_map);
        let d_matrices = Arc::new(HashMap::new());
        Self { summands, d_deg, d_map, d_matrices }
    }

    /// Build from summands and explicit differential matrices (in the summands' generator order);
    /// the matrices are cached, and the differential closure is derived from them.
    pub fn new_with_d_matrices(summands: GrMod<I, X, R>, d_deg: I, matrices: impl IntoIterator<Item = (I, SpMat<R>)>) -> Self {
        let d_matrices: Arc<HashMap<I, SpMat<R>>> = Arc::new(matrices.into_iter().collect());

        let s = summands.clone();
        let ds = d_matrices.clone();
        let mut new = Self::new(summands, d_deg, move |i, z| {
            let Some(d) = ds.get(&i) else {
                return Lc::zero();
            };
            let v = s[i].vectorize(z);
            s[i + d_deg].devectorize(&(d * v))
        });

        new.d_matrices = d_matrices;

        #[cfg(debug_assertions)]
        new.check_d_matrices();

        new
    }

    pub(crate) fn with_d_matrices(mut self, matrices: impl IntoIterator<Item = (I, SpMat<R>)>) -> Self {
        let map: HashMap<I, SpMat<R>> = matrices.into_iter().collect();
        self.d_matrices = Arc::new(map);

        #[cfg(debug_assertions)]
        self.check_d_matrices();

        self
    }

    /// Each cached d-matrix's shape must match the summand ranks at its endpoints.
    /// Callable in release; construction only runs it under `debug_assertions`.
    pub fn check_d_matrices(&self) {
        for (&i, m) in self.d_matrices.iter() {
            let (n_rows, n_cols) = m.shape();
            assert_eq!(n_cols, self[i].rank(),
                "d_matrix at {i}: n_cols {n_cols} != rank(C[{i}]) {}", self[i].rank());
            let j = i + self.d_deg;
            assert_eq!(n_rows, self[j].rank(),
                "d_matrix at {i}: n_rows {n_rows} != rank(C[{j}]) {}", self[j].rank());
        }
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

    // Relabel generators by an injective `f` with inverse `g`; the differential is conjugated by them.
    pub fn map_keys<Y, F, G>(&self, f: F, g: G) -> ChainComplex<I, Y, R>
    where
        Y: LcKey,
        F: Fn(&X) -> Y + Send + Sync + 'static,
        G: Fn(&Y) -> X + Send + Sync + 'static,
    {
        let summands = self.summands.map_keys(&f);
        let d = self.d_map.clone();
        ChainComplex::new(summands, self.d_deg, move |i, z: &Lc<Y, R>| {
            let zx = z.clone().map_keys(|y| g(&y));
            d(i, &zx).map_keys(|x| f(&x))
        })
    }

    /// Returns the differential matrix `d_i: C_i → C_{i + d_deg}` in the
    /// stored basis. Hits the precomputed cache if populated; otherwise builds
    /// the matrix by applying `d_map` to each basis element.
    pub fn d_matrix(&self, i: I) -> SpMat<R> {
        if let Some(d) = self.d_matrices.get(&i) {
            return d.clone();
        }

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

    /// Forget the symbolic generators and return an equivalent
    /// [`GenericChainComplex`] built from the differential matrices.
    pub fn as_generic(&self) -> GenericChainComplex<I, R> {
        GenericChainComplex::from_d_matrices(
            self.d_deg,
            self.support().map(|&i| (i, self.d_matrix(i)))
        )
    }

    /// Reduce the complex via pivot cancellation (Schur complement), preserving
    /// symbolic generators. The returned complex has fewer raw generators but
    /// the same homology; basis-change tracking lets cycles pull back to the
    /// original generators.
    pub fn reduced(&self) -> ChainComplex<I, X, R> {
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

        let matrices = self.summands.support()
            .filter_map(|&i| r.matrix(i).map(|m| (i, m.clone())))
            .collect::<Vec<_>>();

        let d_deg = self.d_deg;
        let d_map = self.d_map.clone();
        Self { summands, d_deg, d_map, d_matrices: Arc::new(HashMap::new()) }
            .with_d_matrices(matrices)
    }

    /// Like [`Self::reduced`] but skips the basis-change tracking and discards
    /// symbolic generators, returning a [`GenericChainComplex`].
    pub fn reduced_generic(&self) -> GenericChainComplex<I, R> {
        let r = ChainReducer::reduce(self, false);
        r.into_generic_complex()
    }

    #[cfg(any(test, feature = "test-utils"))]
    fn check_d_for(&self, i0: I, x: &X) {
        let i1 = i0 + self.d_deg();
        assert!(self.is_supported(i0), "Not supported: {i0}.");

        let dx = self.d(i0, &Lc::from(x.clone()));
        let ddx = self.d(i1, &dx);

        assert!(ddx.is_zero(), "d² is non-zero for {x} at {i0}.\n  dx: {dx}\n  ddx: {ddx}.");
    }

    #[cfg(any(test, feature = "test-utils"))]
    pub fn check_d_at(&self, i0: I) {
        for x in self[i0].raw_generators().iter() {
            self.check_d_for(i0, x);
        }
    }

    #[cfg(any(test, feature = "test-utils"))]
    pub fn check_d_all(&self) {
        for &i in self.support() {
            self.check_d_at(i);
        }
    }
}

impl<I, X, R> ChainComplex<I, X, R>
where
    I: AddInd,
    X: LcKey,
    R: EucRing, for<'x> &'x R: EucRingOps<R>,
{
    pub fn homology_at(&self, i: I) -> Summand<X, R> {
        let c = &self[i];
        let d0 = self.d_matrix(i - self.d_deg());
        let d1 = self.d_matrix(i);
        let (rank, tors, trans) = HomologyCalc::calculate(d0, d1, true);
        let trans = trans.unwrap();

        Summand::new(
            c.raw_generators().clone(),
            rank,
            tors,
            c.trans().merged(&trans)
        )
    }

    pub fn homology(&self) -> GrMod<I, X, R> {
        self.homology_in(self.support().copied())
    }

    /// Homology at the given indices only, using the full differentials.
    pub fn homology_in(&self, support: impl IntoIterator<Item = I>) -> GrMod<I, X, R> {
        GrMod::generate_filtered(
            support.into_iter(),
            |i| {
                let hi = self.homology_at(i);
                (!hi.is_zero()).then_some(hi)
            }
        )
    }

    /// Compute the homology at index `i` without tracking the basis change.
    /// Faster than `homology_at` when only the rank/torsion are needed.
    pub fn generic_homology_at(&self, i: I) -> GenericSummand<I, R> {
        let d0 = self.d_matrix(i - self.d_deg());
        let d1 = self.d_matrix(i);
        let (rank, tors, _) = HomologyCalc::calculate(d0, d1, false);
        GenericSummand::generate(i, rank, tors, None)
    }

    pub fn generic_homology(&self) -> GenericGrMod<I, R> {
        self.generic_homology_in(self.support().copied())
    }

    /// Trans-free homology at the given indices only, using the full differentials.
    pub fn generic_homology_in(&self, support: impl IntoIterator<Item = I>) -> GenericGrMod<I, R> {
        GrMod::generate_filtered(
            support.into_iter(),
            |i| {
                let hi = self.generic_homology_at(i);
                (!hi.is_zero()).then_some(hi)
            }
        )
    }
}

impl<X, R> ChainComplex1<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    pub fn truncated(&self, range: RangeInclusive<isize>) -> Self {
        let d_deg = self.d_deg;
        let d_map = self.d_map.clone();
        let summands = self.summands.truncated(range.clone());

        // cached d-matrices with both endpoints inside the window stay valid.
        let matrices = self.d_matrices.iter().filter_map(|(&i, m)|
            (range.contains(&i) && range.contains(&(i + d_deg))).then(|| (i, m.clone()))
        ).collect_vec();

        Self::new(summands, d_deg, move |i, z|
            if range.contains(&(i + d_deg)) {
                d_map(i, z)
            } else {
                Lc::zero()
            }
        ).with_d_matrices(matrices)
    }
}

impl<I, X, R> Index<I> for ChainComplex<I, X, R>
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

impl<X, R> ToSeqString<isize> for ChainComplex1<X, R>
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

mod tex_impl {
    use super::*;
    use yui_core::TeX;
    use crate::utils::tex::TeXTable;

    impl<X, R> TeXTable<isize2> for ChainComplex2<X, R>
    where X: LcKey, R: Ring + TeX, for<'x> &'x R: RingOps<R> {
        delegate! {
            to self.summands {
                fn tex_table(&self, caption: &str, head: &str) -> String;
            }
        }
    }
}