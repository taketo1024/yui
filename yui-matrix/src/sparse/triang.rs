use std::cell::RefCell;
use std::sync::{Arc, Mutex};
use either::Either;
use log::info;
use rayon::prelude::*;
use thread_local::ThreadLocal;
use yui::{Ring, RingOps};
use super::*;

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
}

pub fn inv_triangular<R>(t: TriangularType, a: &SpMat<R>) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let e = SpMat::id(a.rows());
    solve_triangular(t, a, &e)
}

pub fn solve_triangular<R>(t: TriangularType, a: &SpMat<R>, y: &SpMat<R>) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert_eq!(a.rows(), y.rows());

    if t.is_upper() { 
        debug_assert!(is_upper_tri(a));
    } else { 
        // TODO
    }

    if crate::config::is_multithread_enabled() { 
        solve_triangular_m(t, a, y)
    } else { 
        solve_triangular_s(t, a, y)
    }
}

pub fn solve_triangular_vec<R>(t: TriangularType, a: &SpMat<R>, b: &SpVec<R>) -> SpVec<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert_eq!(a.rows(), b.dim());

    if t.is_upper() { 
        debug_assert!(is_upper_tri(a));
    } else { 
        // TODO
    }

    let n = a.rows();
    let diag = collect_diag(t, a);
    let mut b = b.to_dense().to_vec();

    let x = _solve_triangular(t, a, &diag, &mut b);

    SpVec::from_entries(n, x)
}

fn solve_triangular_s<R>(t: TriangularType, a: &SpMat<R>, y: &SpMat<R>) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    info!("solve triangular: a = {:?}, y = {:?}", a.shape(), y.shape());

    let (n, k) = (a.rows(), y.cols());
    let diag = collect_diag(t, a);
    let mut b = vec![R::zero(); n];

    let entries = (0..k).fold(vec![], |mut entries, j| { 
        copy_into(y.col_view(j).iter(), &mut b);

        let res = _solve_triangular(t, a, &diag, &mut b);
        entries.extend(res.into_iter().map(|(i, x)| (i, j, x)));
        entries
    });

    SpMat::from_entries((n, k), entries)
}

fn solve_triangular_m<R>(t: TriangularType, a: &SpMat<R>, y: &SpMat<R>) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let nth = std::thread::available_parallelism().map(|x| x.get()).unwrap_or(1);

    info!("solve triangular: a = {:?}, y = {:?} (multi-thread: {nth})", a.shape(), y.shape());

    let (n, k) = (a.rows(), y.cols());
    let diag = collect_diag(t, a);
    let tl_b = Arc::new(ThreadLocal::new());
    let entries = Mutex::new(vec![]);

    let report = log::max_level() >= log::LevelFilter::Info && k >= 10_000;
    let col_count = Mutex::new(0);
    
    (0..k).into_par_iter().for_each(|j| { 
        let mut b = tl_b.get_or(|| 
            RefCell::new(vec![R::zero(); n])
        ).borrow_mut();

        copy_into(y.col_view(j).iter(), &mut b);

        let res = _solve_triangular(t, a, &diag, &mut b);
        {
            let mut entries = entries.lock().unwrap();
            entries.extend(res.into_iter().map(|(i, x)| (i, j, x)));
        }

        if report { 
            let c = {
                let mut c = col_count.lock().unwrap();
                *c += 1;
                *c
            };
            if c > 0 && c % 10_000 == 0 { 
                info!("  solved {c}/{k}");
            }
        }
    });

    let entries = entries.into_inner().unwrap();
    SpMat::from_entries((n, k), entries)
}

#[inline(never)] // for profilability
fn _solve_triangular<R>(t: TriangularType, a: &SpMat<R>, diag: &[&R], b: &mut [R]) -> Vec<(usize, R)>
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

    entries
}

fn collect_diag<'a, R>(t: TriangularType, a: &'a SpMat<R>) -> Vec<&'a R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    let a = &a.cs_mat();
    let n = a.rows();
    let indptr = a.indptr();
    let data = a.data();

    if t.is_upper() { 
        (0..n).map( |i| {
            let p = indptr.index(i + 1);
            &data[p - 1]
        }).collect()
    } else { 
        (0..n).map( |i| {
            let p = indptr.index(i);
            &data[p]
        }).collect()
    }
}

fn is_upper_tri<R>(u: &SpMat<R>) -> bool
where R: Ring, for<'x> &'x R: RingOps<R> {
    u.rows() == u.cols() && 
    u.iter().all(|(i, j, a)| 
        i < j || (i == j && a.is_unit()) || (i > j && a.is_zero())
    )
}

fn copy_into<'a, Itr, R>(itr: Itr, x: &mut [R])
where Itr: Iterator<Item = (usize, &'a R)>, R: Clone + 'a { 
    itr.for_each(|(i, r)| x[i] = r.clone())
}

#[cfg(test)]
mod tests { 
    use super::*;
    use super::TriangularType::{Upper, Lower};

    #[test]
    fn solve_upper() { 
        let u = SpMat::from_vec((5, 5), vec![
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
        let u = SpMat::from_vec((5, 5), vec![
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
        let l = SpMat::from_vec((5, 5), vec![
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
        let l = SpMat::from_vec((5, 5), vec![
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