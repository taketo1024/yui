use itertools::Itertools;
use nalgebra_sparse::csc::CscCol;
use yui::UnionFind;

use crate::MatTrait;

use super::SpMat;


pub fn decomp<R>(a: &SpMat<R>) {
    let n = a.ncols();
    if n == 0 { 
        return
    }

    let mut u = UnionFind::new(n);

    (0 .. n-1).for_each(|i|
        (i+1 .. n).for_each(|j| {
            if !u.is_same(i, j) && intersects(a.inner().col(i), a.inner().col(j)) { 
                u.union(i, j)
            }
        })
    );

    let g = u.group();

    if g.len() > 1 { 
        println!("{:?} -> {}", a.shape(), g.iter().map(|l| l.len().to_string()).join(" + "));
    }
}

fn intersects<R>(col1: CscCol<R>, col2: CscCol<R>) -> bool { 
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
            std::cmp::Ordering::Less => {
                if let Some(j1) = itr1.next() { 
                    i1 = j1
                } else { 
                    return false
                }
            },
            std::cmp::Ordering::Equal => {
                return true
            },
            std::cmp::Ordering::Greater => {
                if let Some(j2) = itr2.next() { 
                    i2 = j2
                } else { 
                    return false
                }
            },
        }
    }
}

