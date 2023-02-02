use std::cell::RefCell;
use std::mem::replace;
use std::ops::DerefMut;
use std::sync::{Arc, Mutex};
use either::Either;
use log::info;
use sprs::{CsMat, CsVec, TriMat};
use rayon::prelude::*;
use thread_local::ThreadLocal;
use crate::math::traits::{Ring, RingOps};
use crate::math::matrix::sparse::CsVecExt;

use super::sparse::CsMatExt;

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

pub fn inv_triangular<R>(t: TriangularType, a: &CsMat<R>) -> CsMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let e = CsMat::id(a.rows());
    solve_triangular(t, a, &e)
}

pub fn solve_triangular<R>(t: TriangularType, a: &CsMat<R>, y: &CsMat<R>) -> CsMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert!(a.is_csc());
    assert!(y.is_csc());
    assert_eq!(a.rows(), y.rows());

    if t.is_upper() { 
        debug_assert!(is_upper_tri(a));
    } else { 
        // TODO
    }

    const MULTI_THREAD: bool = true;
    if MULTI_THREAD { 
        solve_triangular_m(t, a, y)
    } else { 
        solve_triangular_s(t, a, y)
    }
}

fn solve_triangular_s<R>(t: TriangularType, a: &CsMat<R>, y: &CsMat<R>) -> CsMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    info!("solve triangular: a = {:?}, y = {:?}", a.shape(), y.shape());

    let (n, k) = (a.rows(), y.cols());
    let diag = diag(t, &a);

    let mut count = 0;
    let mut indptr = vec![0];
    let mut indices = vec![];
    let mut data = vec![];

    let mut x = vec![R::zero(); n];
    let mut b = vec![R::zero(); n];

    for j in 0..k { 
        for (i, r) in y.outer_view(j).unwrap().iter() { 
            b[i] = r.clone();
        }

        solve_triangular_into(t, &a, &diag, &mut b, &mut x);

        for i in 0..n { 
            if x[i].is_zero() { continue }
            let x_i = replace(&mut x[i], R::zero());

            count += 1;
            indices.push(i);
            data.push(x_i);
        }

        indptr.push(count);
    }

    CsMat::new_csc((n, k), indptr, indices, data)
}

fn solve_triangular_m<R>(t: TriangularType, a: &CsMat<R>, y: &CsMat<R>) -> CsMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let nth = std::thread::available_parallelism().map(|x| x.get()).unwrap_or(1);

    info!("solve triangular: a = {:?}, y = {:?} (multi-thread: {nth})", a.shape(), y.shape());

    let report = log::max_level() >= log::LevelFilter::Info;
    let count = Mutex::new(0);
    
    let (n, k) = (a.rows(), y.cols());
    let diag = diag(t, &a);

    let trip = Mutex::new( TriMat::new((n, k)) );
    let tls1 = Arc::new(ThreadLocal::new());
    let tls2 = Arc::new(ThreadLocal::new());

    (0..k).into_par_iter().for_each(|j| { 
        let mut x_st = tls1.get_or(|| {
            let x = vec![R::zero(); n];
            RefCell::new(x)
        }).borrow_mut();

        let mut b_st = tls2.get_or(|| {
            let b = vec![R::zero(); n];
            RefCell::new(b)
        }).borrow_mut();

        let x = x_st.deref_mut();
        let b = b_st.deref_mut();

        for (i, r) in y.outer_view(j).unwrap().iter() { 
            b[i] = r.clone();
        }

        solve_triangular_into(t, &a, &diag, b, x);

        { 
            let mut trip = trip.lock().unwrap();
            for i in 0..n { 
                if x[i].is_zero() { continue }
                let x_i = replace(&mut x[i], R::zero());
                trip.add_triplet(i, j, x_i);
            }
        }

        if report { 
            let c = {
                let mut count = count.lock().unwrap();
                *count += 1;
                *count
            };

            if c > 0 && c % 10_000 == 0 { 
                info!("  solved {c}/{k}");
            }
        }
    });

    trip.into_inner().unwrap().to_csc()
}

pub fn solve_triangular_vec<R>(t: TriangularType, a: &CsMat<R>, b: &CsVec<R>) -> CsVec<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert!(a.is_csc());
    assert_eq!(a.rows(), b.dim());

    if t.is_upper() { 
        debug_assert!(is_upper_tri(a));
    } else { 
        // TODO
    }

    let n = a.rows();
    let diag = diag(t, &a);

    let mut x = vec![R::zero(); n];
    let mut b = b.to_dense().to_vec();

    solve_triangular_into(t, a, &diag, &mut b, &mut x);

    CsVec::from_vec(x)
}

fn solve_triangular_into<R>(t: TriangularType, a: &CsMat<R>, diag: &Vec<&R>, b: &mut Vec<R>, x: &mut Vec<R>)
where R: Ring, for<'x> &'x R: RingOps<R> {
    debug_assert!(x.iter().all(|x_i| 
        x_i.is_zero())
    ); 

    let n = a.rows();
    let range = if t.is_upper() { 
        Either::Left((0..n).rev())
    } else { 
        Either::Right(0..n)
    };

    for j in range {
        if b[j].is_zero() { continue }

        let u_jj_inv = diag[j].inv().unwrap();
        let x_j = &b[j] * &u_jj_inv; // non-zero

        for (i, u_ij) in a.outer_view(j).unwrap().iter() {
            if u_ij.is_zero() { continue }
            b[i] -= u_ij * &x_j;
        }

        x[j] = x_j;
    }

    debug_assert!(b.iter().all(|b_i| 
        b_i.is_zero())
    );
}

fn diag<'a, R>(t: TriangularType, a: &'a CsMat<R>) -> Vec<&'a R> {
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

fn is_upper_tri<R>(u: &CsMat<R>) -> bool
where R: Ring, for<'x> &'x R: RingOps<R> {
    u.rows() == u.cols() && 
    u.iter().all(|(a, (i, j))| 
        i < j || (i == j && a.is_unit()) || (i > j && a.is_zero())
    )
}

#[cfg(test)]
mod tests { 
    use sprs::CsVec;

    use crate::math::matrix::sparse::{CsMatExt, CsVecExt};
    use super::*;
    use super::TriangularType::{Upper, Lower};

    #[test]
    fn solve_upper() { 
        let u = CsMat::csc_from_vec((5, 5), vec![
            1, -2, 1,  3, 5,
            0, -1, 4,  2, 1,
            0,  0, 1,  0, 3,
            0,  0, 0, -1, 5,
            0,  0, 0,  0, 1
        ]);
        let x = CsVec::from_vec(vec![1,2,3,4,5]);
        let b = CsVec::from_vec(vec![37,23,18,21,5]);
        assert_eq!(solve_triangular_vec(Upper, &u, &b), x);
    }

    #[test]
    fn inv_upper() { 
        let u = CsMat::csc_from_vec((5, 5), vec![
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
        let l = CsMat::csc_from_vec((5, 5), vec![
            1,  0, 0,  0, 0,
           -2, -1, 0,  0, 0,
            1,  4, 1,  0, 0,
            3,  2, 0, -1, 0,
            5,  1, 3,  5, 1
        ]);
        let x = CsVec::from_vec(vec![1,2,3,4,5]);
        let b = CsVec::from_vec(vec![1,-4,12,3,41]);
        assert_eq!(solve_triangular_vec(Lower, &l, &b), x);
    }

    #[test]
    fn inv_lower() { 
        let l = CsMat::csc_from_vec((5, 5), vec![
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