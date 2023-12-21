use itertools::Itertools;
use yui::{Ring, RingOps, IndexList};
use yui::lc::{Lc, Gen};
use yui_matrix::sparse::SpMat;

#[cfg(feature = "multithread")]
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

pub fn make_matrix<X, Y, R, F>(from: &IndexList<X>, to: &IndexList<Y>, f: F) -> SpMat<R>
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
pub fn make_matrix_async<X, Y, R, F>(from: &IndexList<X>, to: &IndexList<Y>, f: F) -> SpMat<R>
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