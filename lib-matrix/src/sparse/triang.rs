use either::Either;
use log::*;
use yui_core::{Ring, RingOps};

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

    pub fn transpose(&self) -> Self {
        match self { 
            Self::Upper => Self::Lower,
            Self::Lower => Self::Upper
        }
    }

    fn str(&self) -> &'static str { 
        match self { 
            Self::Upper => "upper",
            Self::Lower => "lower"
        }
    }
}

pub fn inv_triangular<R>(t: TriangularType, a: &SpMat<R>) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let e = SpMat::id(a.n_rows());
    solve_triangular(t, a, &e)
}

// solve ax = y.
pub fn solve_triangular<R>(t: TriangularType, a: &SpMat<R>, y: &SpMat<R>) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let n = a.n_rows();
    let cols = solve_triangular_with(t, a, y, |_, x| x);
    SpMat::from_col_vecs(n, cols)
}

// Solve `ax = y` column by column, invoking `f(j, x_j)` on each solved
// column. Hoists the diagonal collection and RHS buffer out of the per-column
// loop, and (under `multithread`) reuses one buffer per thread.
pub(crate) fn solve_triangular_with<R, F, T>(
    t: TriangularType, a: &SpMat<R>, y: &SpMat<R>, f: F
) -> Vec<T>
where
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(usize, SpVec<R>) -> T + Sync,
    T: Send,
{
    assert_eq!(a.n_rows(), y.n_rows());
    debug_assert!(a.is_triang(t));

    cfg_if::cfg_if! {
        if #[cfg(feature = "multithread")] {
            solve_triangular_m(t, a, y, f)
        } else {
            solve_triangular_s(t, a, y, f)
        }
    }
}

// solve xa = y.
pub fn solve_triangular_left<R>(t: TriangularType, a: &SpMat<R>, y: &SpMat<R>) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    solve_triangular(t.transpose(), &a.transpose(), &y.transpose()).transpose()
}

pub fn solve_triangular_vec<R>(t: TriangularType, a: &SpMat<R>, b: &SpVec<R>) -> SpVec<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert_eq!(a.n_rows(), b.dim());
    debug_assert!(a.is_triang(t));

    debug!("solve {} triangular-vec", t.str());
    debug!("  a: {:?}", a.shape());

    let n = a.n_rows();
    let diag = collect_diag(t, a);
    let mut b_buf = vec![R::zero(); n];
    scatter_into(b.data(), &mut b_buf);

    _solve_triangular(t, a, &diag, &mut b_buf)
}

#[allow(unused)]
fn solve_triangular_s<R, F, T>(
    t: TriangularType, a: &SpMat<R>, y: &SpMat<R>, f: F
) -> Vec<T>
where
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(usize, SpVec<R>) -> T,
{
    debug!("solve {} triangular", t.str());
    debug!("  a: {:?}, y: {:?}", a.shape(), y.shape());

    let (n, k) = (a.n_rows(), y.n_cols());
    let diag = collect_diag(t, a);
    let mut b = vec![R::zero(); n];

    (0..k).map(|j| {
        scatter_into(y.col_data(j), &mut b);
        let x = _solve_triangular(t, a, &diag, &mut b);
        f(j, x)
    }).collect()
}

#[cfg(feature = "multithread")]
fn solve_triangular_m<R, F, T>(
    t: TriangularType, a: &SpMat<R>, y: &SpMat<R>, f: F
) -> Vec<T>
where
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(usize, SpVec<R>) -> T + Sync,
    T: Send,
{
    use yui_core::util::sync::SyncCounter;

    debug!("solve {} triangular (threads: {})", t.str(), rayon::current_num_threads());
    debug!("  a: {:?}, y: {:?}", a.shape(), y.shape());

    let (n, k) = (a.n_rows(), y.n_cols());
    let diag = collect_diag(t, a);
    let tl_b = Arc::new(ThreadLocal::new());

    let report = should_report(y);
    let counter = SyncCounter::new();

    (0..k).into_par_iter().map(|j| {
        let mut b = tl_b.get_or(||
            RefCell::new(vec![R::zero(); n])
        ).borrow_mut();

        scatter_into(y.col_data(j), &mut b);
        let x = _solve_triangular(t, a, &diag, &mut b);
        let result = f(j, x);

        if report {
            let c = counter.incr();
            if (c > 0 && c % LOG_THRESHOLD == 0) || c == k {
                trace!("  solved {c}/{k}.");
            }
        }

        result
    }).collect()
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

        let (idx, val) = a.col_data(j);
        for (&i, a_ij) in idx.iter().zip(val.iter()) {
            if a_ij.is_zero() { continue }
            b[i] -= a_ij * &x_j;
        }

        entries.push((j, x_j));
    }

    debug_assert!(b.iter().all(|b_i|
        b_i.is_zero())
    );

    let entries = if t.is_upper() {
        Either::Left(entries.into_iter().rev())
    } else {
        Either::Right(entries.into_iter())
    };

    SpVec::from_sorted_entries(a.n_cols(), entries)
}

fn collect_diag<'a, R>(t: TriangularType, a: &'a SpMat<R>) -> Vec<&'a R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (col_offsets, row_indices, values) = a.csc_data();
    (0..a.n_cols()).map(|j| {
        let p = if t.is_upper() {
            col_offsets[j + 1] - 1
        } else {
            col_offsets[j]
        };
        assert_eq!(row_indices[p], j, "broken input: missing diagonal at column {j}");
        &values[p]
    }).collect()
}

fn scatter_into<R: Clone>(data: (&[usize], &[R]), dst: &mut [R]) {
    let (idx, val) = data;
    for (&i, v) in idx.iter().zip(val.iter()) {
        dst[i] = v.clone();
    }
}

#[allow(unused)]
fn should_report<R>(a: &SpMat<R>) -> bool { 
    usize::min(a.n_rows(), a.n_cols()) > LOG_THRESHOLD && log::max_level() >= log::LevelFilter::Debug
}

#[cfg(test)]
mod tests { 
    use super::*;
    use super::TriangularType::{Upper, Lower};

    #[test]
    fn solve_upper() { 
        let u = SpMat::from_row_major((5, 5), vec![
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
        let u = SpMat::from_row_major((5, 5), [
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
        let l = SpMat::from_row_major((5, 5), [
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
        let l = SpMat::from_row_major((5, 5), [
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