use std::hash::Hash;

use itertools::Itertools;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use yui_core::{Ring, RingOps, IndexList};
use yui_matrix::sparse::SpMat;

pub fn make_matrix<X, Y, R, F>(from: &IndexList<X>, to: &IndexList<Y>, f: F) -> SpMat<R>
where 
    X: Hash + Eq, Y: Hash + Eq, 
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(&X) -> Vec<(Y, R)> 
{
    let (m, n) = (to.len(), from.len());

    let entries = (0..n).flat_map(|j| {
        let x = &from[j];
        let fx = f(x);
        
        fx.iter().filter_map(|(y, a)| { 
            if a.is_zero() { 
                return None
            }
            let i = to.index_of(&y).unwrap();
            Some((i, j, a.clone()))
        }).collect_vec()
    });
    
    SpMat::from_entries((m, n), entries)
}

pub fn make_matrix_async<X, Y, R, F>(from: &IndexList<X>, to: &IndexList<Y>, f: F) -> SpMat<R>
where 
    X: Hash + Eq + Send + Sync, 
    Y: Hash + Eq + Send + Sync, 
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(&X) -> Vec<(Y, R)> + Send + Sync
{
    let (m, n) = (to.len(), from.len());

    let entries = (0..n).into_par_iter().flat_map(|j| {
        let x = &from[j];
        let ys = f(x);
        ys.iter().filter_map(|(y, a)| { 
            if a.is_zero() { 
                return None
            }
            let i = to.index_of(&y).unwrap();
            Some((i, j, a.clone()))
        }).collect_vec()
    }).collect::<Vec<_>>();

    SpMat::from_entries((m, n), entries)
}