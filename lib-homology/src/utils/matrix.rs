use itertools::Itertools;
use yui_core::{Ring, RingOps, IndexList};
use yui_core::lc::{Lc, Gen};
use yui_matrix::sparse::SpMat;

#[cfg(feature = "multithread")]
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

pub fn make_matrix<X, Y, R, F>(from: &IndexList<X>, to: &IndexList<Y>, f: F) -> SpMat<R>
where 
    X: Gen, Y: Gen, 
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(&X) -> Lc<Y, R> + Send + Sync
{
    cfg_if::cfg_if! { 
        if #[cfg(feature="multithread")] { 
            make_matrix_m(from, to, f)
        } else { 
            make_matrix_s(from, to, f)
        }
    }
}

pub fn make_matrix_s<X, Y, R, F>(from: &IndexList<X>, to: &IndexList<Y>, f: F) -> SpMat<R>
where 
    X: Gen, Y: Gen, 
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(&X) -> Lc<Y, R> 
{
    let (m, n) = (to.len(), from.len());

    let entries = (0..n).flat_map(|j| {
        let x = &from[j];
        let fx = f(x);
        
        fx.iter().map(|(y, a)| { 
            let i = to.index_of(y).unwrap();
            (i, j, a.clone())
        }).collect_vec()
    });
    
    SpMat::from_entries((m, n), entries)
}

#[cfg(feature = "multithread")]
pub fn make_matrix_m<X, Y, R, F>(from: &IndexList<X>, to: &IndexList<Y>, f: F) -> SpMat<R>
where 
    X: Gen, Y: Gen, 
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(&X) -> Lc<Y, R> + Send + Sync
{
    let (m, n) = (to.len(), from.len());

    let entries = (0..n).into_par_iter().flat_map(|j| {
        let x = &from[j];
        let ys = f(x);
        ys.iter().map(|(y, a)| { 
            let i = to.index_of(y).unwrap();
            (i, j, a.clone())
        }).collect_vec()
    }).collect::<Vec<_>>();

    SpMat::from_entries((m, n), entries)
}