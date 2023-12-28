use std::fmt::Display;

use ahash::AHashSet;
use itertools::Itertools;
use nalgebra::{Scalar, ClosedAdd};
use num_traits::Zero;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use sprs::PermOwned;
use yui::UnionFind;

cfg_if::cfg_if! { if #[cfg(feature = "multithread")] { 
    use std::sync::Mutex;
}}

use crate::MatTrait;

use super::SpMat;
use super::util::perm_for_indices;

pub fn dir_sum_indices<R>(a: &SpMat<R>) -> Option<(Vec<Vec<usize>>, Vec<Vec<usize>>)>
where R: Clone + Scalar + Zero + ClosedAdd + Send + Sync + Display {
    let (m, n) = a.shape();
    let cols = group_cols(a);
    let rows = cols.iter().map(|c| rows_in(a, c)).collect_vec();

    if (rows.len(), cols.len()) == (1, 1) && (rows[0].len(), cols[0].len()) == (m, n) { 
        None
    } else { 
        Some((rows, cols))
    }
}

pub fn dir_sum_decomp<R>(a: &SpMat<R>) -> Option<(Vec<SpMat<R>>, PermOwned, PermOwned)>
where R: Clone + Scalar + Zero + ClosedAdd + Send + Sync + Display {
    let Some((rows, cols)) = dir_sum_indices(a) else { 
        return None
    };

    let (m, n) = a.shape();
    let p = perm_for_indices(m, rows.iter().flatten());
    let q = perm_for_indices(n, cols.iter().flatten());
    let s = decomp_by(a, &rows, &cols, &p, &q);

    // println!("{} -> {} summands:\n{}-----", a.clone().into_dense(), s.len(), s.iter().map(|s| s.clone().into_dense()).join("  (+)\n"));

    Some((s, p, q))
}

fn group_cols<R>(a: &SpMat<R>) -> Vec<Vec<usize>>
where R: Zero + Send + Sync {
    // non-empty cols.
    let cols = (0..a.ncols()).filter(|&j|
        !a.inner().col(j).row_indices().is_empty()
    ).collect_vec();
    
    let l = cols.len();
    if l == 0 { // prevent panic caused by `l - 1`
        return vec![]
    }

    let u = UnionFind::new(l);
    let res: Vec<Vec<usize>>;

    cfg_if::cfg_if! { 
        if #[cfg(feature = "multithread")] { 
            let u = Mutex::new(u);

            (0 .. l - 1).into_par_iter().for_each(|i|
                (i + 1 .. l).into_par_iter().for_each(|j| {
                    if !u.lock().unwrap().is_same(i, j) && col_intersects(a, cols[i], cols[j]) { 
                        u.lock().unwrap().union(i, j)
                    }
                })
            );

            res = u.into_inner().unwrap().group();
        } else { 
            let mut u = u;

            (0 .. l - 1).for_each(|i|
                (i + 1 .. l).for_each(|j| {
                    println!("{}, {} -> {}", cols[i], cols[j], col_intersects(a, cols[i], cols[j]));
                    if !u.is_same(i, j) && col_intersects(a, cols[i], cols[j]) { 
                        u.union(i, j)
                    }
                })
            );

            res = u.group();
        }
    }

    res.into_iter().map(|list|
        list.into_iter().map(|i| cols[i]).collect()
    ).collect()
}

fn col_intersects<R>(a: &SpMat<R>, j1: usize, j2: usize) -> bool { 
    use std::cmp::Ordering::*;

    let col1 = a.inner().col(j1);
    let col2 = a.inner().col(j2);
    let mut itr1 = col1.row_indices().iter();
    let mut itr2 = col2.row_indices().iter();

    let Some(mut i1) = itr1.next() else { 
        return false
    };
    let Some(mut i2) = itr2.next() else { 
        return false
    };

    loop { 
        match usize::cmp(&i1, &i2) {
            Less => {
                if let Some(j1) = itr1.next() { 
                    i1 = j1
                } else { 
                    return false
                }
            },
            Equal => {
                return true
            },
            Greater => {
                if let Some(j2) = itr2.next() { 
                    i2 = j2
                } else { 
                    return false
                }
            },
        }
    }
}

fn rows_in<R>(a: &SpMat<R>, cols: &[usize]) -> Vec<usize> { 
    let mut set = AHashSet::new();

    cols.iter().for_each(|&j| { 
        set.extend(a.inner().col(j).row_indices());
    });

    set.into_iter().sorted().collect()
}

fn decomp_by<R>(a: &SpMat<R>, rows: &Vec<Vec<usize>>, cols: &Vec<Vec<usize>>, p: &PermOwned, q: &PermOwned) -> Vec<SpMat<R>>
where R: Clone + Zero + Scalar + ClosedAdd {
    assert_eq!(rows.len(), cols.len());

    let l = rows.len();
    let row_inds = rows.iter().map(|r| AHashSet::from_iter(r)).collect_vec();
    let row_offsets = rows.iter().fold(vec![0], |mut res, next| { 
        let offset = res.last().unwrap();
        res.push(offset + next.len());
        res
    });
    let col_offsets = cols.iter().fold(vec![0], |mut res, next| { 
        let offset = res.last().unwrap();
        res.push(offset + next.len());
        res
    });
    let mut entries = vec![vec![]; l];
    a.iter().for_each(|(i, j, a)| { 
        let Some(k) = row_inds.iter().position(|s| s.contains(&i)) else { 
            panic!();
        };
        let i = p.at(i) - row_offsets[k];
        let j = q.at(j) - col_offsets[k];
        let e = (i, j, a.clone());

        entries[k].push(e)
    });

    entries.into_iter().enumerate().map(|(k, e)| {
        let shape = (row_offsets[k + 1] - row_offsets[k], col_offsets[k + 1] - col_offsets[k]);
        SpMat::from_entries(shape, e)
    }).collect()
}

#[cfg(test)]
mod tests {
    use crate::sparse::SpMat;

    use super::dir_sum_decomp;
 
    #[test]
    fn test() { 
        let a = SpMat::from_dense_data((5, 7), [
            0, 0, 0, 1, 0, 0, 0, 
            0, 0, 3, 0, 0, 0, 1, 
            1, 0, 0, 2, 0, 0, 0,  
            0, 0, 0, 0, 1, 0, 0, 
            0, 0, 0, 0, 0, 0, 0, 
        ]);
        let decomp = dir_sum_decomp(&a);

        assert!(decomp.is_some());

        let (s, p, q) = decomp.unwrap();

        assert_eq!(s.len(), 3);
        assert_eq!(s, vec![
            SpMat::from_dense_data((2, 2), [0, 1, 1, 2]),
            SpMat::from_dense_data((1, 2), [3, 1]),
            SpMat::from_dense_data((1, 1), [1])
        ]);

        let b = a.permute(p.view(), q.view());

        assert_eq!(b, SpMat::from_dense_data((5, 7), [
            0, 1, 0, 0, 0, 0, 0,
            1, 2, 0, 0, 0, 0, 0,
            0, 0, 3, 1, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0,
        ]));
    }

    #[test]
    fn test2() { 
        let a = SpMat::from_dense_data((9, 10), [
             0,  0,  0,  0, -1,  0,  0,  0,  0, -1,
            -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,
             0,  0, -1,  0,  0,  0,  0,  0,  0,  0,
             0,  0,  0,  0,  0, -1, -1,  0,  0,  0,
             0,  0,  0, -1,  0,  0,  0,  0,  0,  0,
             1,  0,  0,  0,  0,  0,  0, -1, -1,  0,
             0,  0,  0,  1,  1,  1,  0,  0,  0,  0,
             0,  1,  0,  0,  0,  0,  1,  1,  0,  0,
             0,  0,  1,  0,  0,  0,  0,  0,  1,  1,
        ]);
        let decomp = dir_sum_decomp(&a);
        assert!(decomp.is_none());
    }
}