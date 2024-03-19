use either::Either;
use log::debug;
use num_traits::Zero;
use yui::{Ring, RingOps};

use super::*;

cfg_if::cfg_if! {
    if #[cfg(feature = "multithread")] {
        use std::cell::RefCell;
        use std::sync::Arc;
        use thread_local::ThreadLocal;
        use rayon::prelude::*;
    }
}

const LOG_THRESHOLD: usize = 10_000;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum TriangularType { 
    Upper, Lower
}

impl TriangularType { 
    pub fn is_upper(&self) -> bool { 
        match self { 
            Self::Upper => true,
            Self::Lower => false
        }
    }

    pub fn tranpose(&self) -> Self { 
        match self { 
            Self::Upper => Self::Lower,
            Self::Lower => Self::Upper
        }
    }
}

pub fn inv_triangular<R>(t: TriangularType, a: &SpMat<R>) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let e = SpMat::id(a.nrows());
    solve_triangular(t, a, &e)
}

// solve ax = y.
pub fn solve_triangular<R>(t: TriangularType, a: &SpMat<R>, y: &SpMat<R>) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert_eq!(a.nrows(), y.nrows());
    debug_assert!(a.is_triang(t));

    cfg_if::cfg_if! { 
        if #[cfg(feature = "multithread")] { 
            solve_triangular_m(t, a, y)
        } else { 
            solve_triangular_s(t, a, y)
        }
    }
}

// solve xa = y.
pub fn solve_triangular_left<R>(t: TriangularType, a: &SpMat<R>, y: &SpMat<R>) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    solve_triangular(t.tranpose(), &a.transpose(), &y.transpose()).transpose()
}

pub fn solve_triangular_vec<R>(t: TriangularType, a: &SpMat<R>, b: &SpVec<R>) -> SpVec<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert_eq!(a.nrows(), b.dim());
    debug_assert!(a.is_triang(t));

    let diag = collect_diag(a);
    let mut b = b.to_dense();

    _solve_triangular(t, a, &diag, &mut b)
}

#[allow(unused)]
fn solve_triangular_s<R>(t: TriangularType, a: &SpMat<R>, y: &SpMat<R>) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    debug!("solve triangular, y: {:?}", y.shape());

    let (n, k) = (a.nrows(), y.ncols());
    let diag = collect_diag(a);
    let mut b = vec![R::zero(); n];

    let cols = (0..k).map(|j| { 
        copy_into(y.col_vec(j), &mut b);
        _solve_triangular(t, a, &diag, &mut b)
    });

    SpMat::from_col_vecs(n, cols)
}

#[cfg(feature = "multithread")]
fn solve_triangular_m<R>(t: TriangularType, a: &SpMat<R>, y: &SpMat<R>) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    use yui::util::sync::SyncCounter;

    debug!("solve triangular, y: {:?}", y.shape());

    let (n, k) = (a.nrows(), y.ncols());
    let diag = collect_diag(a);
    let tl_b = Arc::new(ThreadLocal::new());

    let report = should_report(y);
    let counter = SyncCounter::new();

    let cols = (0..k).into_par_iter().map(|j| { 
        let mut b = tl_b.get_or(|| 
            RefCell::new(vec![R::zero(); n])
        ).borrow_mut();

        copy_into(y.col_vec(j), &mut b);
        let col = _solve_triangular(t, a, &diag, &mut b);

        if report { 
            let c = counter.incr();
            if (c > 0 && c % LOG_THRESHOLD == 0) || c == k { 
                debug!("  solved {c}/{k}.");
            }
        }

        col
    }).collect::<Vec<_>>();

    SpMat::from_col_vecs(n, cols)
}

#[inline(never)] // for profilability
fn _solve_triangular<R>(t: TriangularType, a: &SpMat<R>, diag: &[&R], b: &mut [R]) -> SpVec<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let mut entries = vec![];

    let itr = diag.iter().enumerate();
    let itr = if t.is_upper() { 
        Either::Left(itr.rev())
    } else { 
        Either::Right(itr)
    };

    for (j, u) in itr { // u = a_jj
        if b[j].is_zero() { continue }

        let uinv = u.inv().unwrap();
        let x_j = &b[j] * &uinv; // non-zero

        for (i, a_ij) in a.col_vec(j).iter() {
            if a_ij.is_zero() { continue }
            b[i] -= a_ij * &x_j;
        }

        entries.push((j, x_j));
    }

    debug_assert!(b.iter().all(|b_i| 
        b_i.is_zero())
    );

    if t.is_upper() { 
        entries.reverse()
    };

    SpVec::from_sorted_entries(a.ncols(), entries)
}

fn collect_diag<'a, R>(a: &'a SpMat<R>) -> Vec<&'a R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    a.iter().filter_map(|(i, j, a)| 
        if i == j { Some(a) } else { None }
    ).collect()
}

fn copy_into<R>(vec: SpVec<R>, x: &mut [R])
where R: Clone + Zero { 
    vec.iter().for_each(|(i, r)| x[i] = r.clone())
}

#[allow(unused)]
fn should_report<R>(a: &SpMat<R>) -> bool { 
    usize::min(a.nrows(), a.ncols()) > LOG_THRESHOLD && log::max_level() >= log::LevelFilter::Trace
}

#[cfg(test)]
mod tests { 
    use super::*;
    use super::TriangularType::{Upper, Lower};

    #[test]
    fn solve_upper() { 
        let u = SpMat::from_dense_data((5, 5), vec![
            1, -2, 1,  3, 5,
            0, -1, 4,  2, 1,
            0,  0, 1,  0, 3,
            0,  0, 0, -1, 5,
            0,  0, 0,  0, 1
        ]);
        let x = SpVec::from(vec![1,2,3,4,5]);
        let b = SpVec::from(vec![37,23,18,21,5]);
        assert_eq!(solve_triangular_vec(Upper, &u, &b), x);
    }

    #[test]
    fn inv_upper() { 
        let u = SpMat::from_dense_data((5, 5), [
            1, -2, 1,  3, 5,
            0, -1, 4,  2, 1,
            0,  0, 1,  0, 3,
            0,  0, 0, -1, 5,
            0,  0, 0,  0, 1
        ]);
        let uinv = inv_triangular(Upper, &u);
        let e = &u * &uinv;
        assert!(e.is_id());
    }

    #[test]
    fn solve_lower() { 
        let l = SpMat::from_dense_data((5, 5), [
            1,  0, 0,  0, 0,
           -2, -1, 0,  0, 0,
            1,  4, 1,  0, 0,
            3,  2, 0, -1, 0,
            5,  1, 3,  5, 1
        ]);
        let x = SpVec::from(vec![1,2,3,4,5]);
        let b = SpVec::from(vec![1,-4,12,3,41]);
        assert_eq!(solve_triangular_vec(Lower, &l, &b), x);
    }

    #[test]
    fn inv_lower() { 
        let l = SpMat::from_dense_data((5, 5), [
            1,  0, 0,  0, 0,
           -2, -1, 0,  0, 0,
            1,  4, 1,  0, 0,
            3,  2, 0, -1, 0,
            5,  1, 3,  5, 1
        ]);
        let linv = inv_triangular(Lower, &l);
        let e = &l * &linv;
        assert!(e.is_id());
    }
}