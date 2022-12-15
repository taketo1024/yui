use std::mem::replace;

use itertools::Itertools;
use sprs::{CsMat, TriMat, CsVec};

use crate::math::{traits::{Ring, RingOps}, matrix::sparse::CsVecExt};

use super::CsMatElem;

pub fn inv_upper_tri<R>(u: &CsMat<R>) -> CsMat<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> {
    assert!(u.is_csc());
    debug_assert!(is_upper_tri(u));

    let n = u.rows();
    let zero = R::zero();
    let diag = diag(&u, &zero);

    let mut trip = TriMat::new(u.shape());  // TODO: directly construct CsMat.
    let mut x = vec![R::zero(); n];
    let mut e = vec![R::zero(); n];

    for j in 0..n { 
        e[j] = R::one();

        solve_upper_tri_into(&u, &diag, &mut e, &mut x);

        for i in 0..n { 
            if x[i].is_zero() { continue }

            let x_i = replace(&mut x[i], R::zero());
            trip.add_triplet(i, j, x_i);
        }
    }

    trip.to_csc()
}

pub fn solve_upper_tri<R>(u: &CsMat<R>, b: &CsVec<R>) -> CsVec<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> {
    assert!(u.is_csc());
    debug_assert!(is_upper_tri(u));

    let n = u.rows();
    let zero = R::zero();
    let diag = diag(&u, &zero);

    let mut x = vec![R::zero(); n];
    let mut b = b.to_dense().to_vec();

    solve_upper_tri_into(u, &diag, &mut b, &mut x);

    CsVec::from_vec(x)
}

fn solve_upper_tri_into<R>(u: &CsMat<R>, diag: &Vec<&R>, b: &mut Vec<R>, x: &mut Vec<R>)
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> {
    debug_assert!(x.iter().all(|x_i| 
        x_i.is_zero())
    ); 

    let n = u.rows();
    for j in (0..n).rev() {
        if b[j].is_zero() { continue }

        let Some(u_jj_inv) = diag[j].inv() else { panic!() };
        let x_j = &b[j] * &u_jj_inv; // non-zero

        for (i, u_ij) in u.outer_view(j).unwrap().iter() {
            if u_ij.is_zero() { continue }
            b[i] -= u_ij * &x_j;
        }

        x[j] = x_j;
    }

    debug_assert!(b.iter().all(|b_i| 
        b_i.is_zero())
    );
}

fn diag<'a, R>(u: &'a CsMat<R>, zero: &'a R) -> Vec<&'a R> {
    (0..u.rows()).map( |i| 
        u.get(i, i).unwrap_or(&zero)
    ).collect_vec()
}

fn is_upper_tri<R>(u: &CsMat<R>) -> bool
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> {
    u.rows() == u.cols() && 
    u.iter().all(|(a, (i, j))| 
        i < j || (i == j && a.is_unit()) || (i > j && a.is_zero())
    )
}

#[cfg(test)]
mod tests { 
    use num_traits::Zero;
    use sprs::CsVec;

    use crate::math::{matrix::sparse::{CsMatExt, CsVecExt}, traits::RingMethods};

    use super::*;

    #[test]
    fn solve() { 
        let u = CsMat::csc_from_vec((5, 5), vec![
            1, -2, 1,  3, 5,
            0, -1, 4,  2, 1,
            0,  0, 1,  0, 3,
            0,  0, 0, -1, 5,
            0,  0, 0,  0, 1
        ]);
        let x = CsVec::from_vec(vec![1,2,3,4,5]);
        let b = CsVec::from_vec(vec![37,23,18,21,5]);
        assert_eq!(solve_upper_tri(&u, &b), x);
    }

    #[test]
    fn inv() { 
        let u = CsMat::csc_from_vec((5, 5), vec![
            1, -2, 1,  3, 5,
            0, -1, 4,  2, 1,
            0,  0, 1,  0, 3,
            0,  0, 0, -1, 5,
            0,  0, 0,  0, 1
        ]);
        let uinv = inv_upper_tri(&u);
        let e = &u * &uinv;

        assert!(e.iter().all(|(a, (i, j))|
            (i == j && a.is_unit()) || (i != j && a.is_zero())
        ))
    }
}