use std::hash::BuildHasher;
use itertools::Itertools;
use indexmap::IndexSet;
use yui_core::{Ring, RingOps};
use yui_core::lc::{Lc, LcKey};
use yui_matrix::sparse::SpMat;

#[cfg(feature = "multithread")]
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

pub fn make_matrix<X, Y, R, F, S1, S2>(from: &IndexSet<X, S1>, to: &IndexSet<Y, S2>, f: F) -> SpMat<R>
where
    X: LcKey, Y: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(&X) -> Lc<Y, R> + Send + Sync,
    S1: BuildHasher + Sync, S2: BuildHasher + Sync,
{
    cfg_if::cfg_if! {
        if #[cfg(feature="multithread")] {
            make_matrix_m(from, to, f)
        } else {
            make_matrix_s(from, to, f)
        }
    }
}

#[cfg(not(feature = "multithread"))]
pub(crate) fn make_matrix_s<X, Y, R, F, S1, S2>(from: &IndexSet<X, S1>, to: &IndexSet<Y, S2>, f: F) -> SpMat<R>
where
    X: LcKey, Y: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(&X) -> Lc<Y, R>,
    S1: BuildHasher, S2: BuildHasher,
{
    let (m, n) = (to.len(), from.len());

    let entries = (0..n).flat_map(|j| {
        let x = &from[j];
        let fx = f(x);

        fx.iter().map(|(y, a)| {
            let i = to.get_index_of(y).unwrap();
            (i, j, a.clone())
        }).collect_vec()
    });

    SpMat::from_entries((m, n), entries)
}

#[cfg(feature = "multithread")]
pub(crate) fn make_matrix_m<X, Y, R, F, S1, S2>(from: &IndexSet<X, S1>, to: &IndexSet<Y, S2>, f: F) -> SpMat<R>
where
    X: LcKey, Y: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(&X) -> Lc<Y, R> + Send + Sync,
    S1: BuildHasher + Sync, S2: BuildHasher + Sync,
{
    let (m, n) = (to.len(), from.len());

    let entries = (0..n).into_par_iter().flat_map(|j| {
        let x = &from[j];
        let ys = f(x);
        ys.iter().map(|(y, a)| {
            let i = to.get_index_of(y).unwrap();
            (i, j, a.clone())
        }).collect_vec()
    }).collect::<Vec<_>>();

    SpMat::from_entries((m, n), entries)
}
