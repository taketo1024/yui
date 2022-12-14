use std::mem::replace;
use itertools::Itertools;
use sprs::{CsMat, CsVec};
use crate::math::traits::{Ring, RingOps};
use crate::math::matrix::sparse::CsVecExt;

pub fn inv_upper_tri<R>(u: &CsMat<R>) -> CsMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert!(u.is_csc());
    debug_assert!(is_upper_tri(u));

    let shape = u.shape();
    let n = u.rows();
    let diag = diag(&u);

    let mut count = 0;
    let mut indptr = vec![0];
    let mut indices = vec![];
    let mut data = vec![];

    let mut x = vec![R::zero(); n];
    let mut e = vec![R::zero(); n];

    for j in 0..n { 
        e[j] = R::one();

        solve_upper_tri_into(&u, &diag, &mut e, &mut x);

        for i in 0..=j { 
            if x[i].is_zero() { continue }
            let x_i = replace(&mut x[i], R::zero());

            count += 1;
            indices.push(i);
            data.push(x_i);
        }

        indptr.push(count);
    }

    CsMat::new_csc(shape, indptr, indices, data)
}

pub fn solve_upper_tri<R>(u: &CsMat<R>, b: &CsVec<R>) -> CsVec<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert!(u.is_csc());
    debug_assert!(is_upper_tri(u));

    let n = u.rows();
    let diag = diag(&u);

    let mut x = vec![R::zero(); n];
    let mut b = b.to_dense().to_vec();

    solve_upper_tri_into(u, &diag, &mut b, &mut x);

    CsVec::from_vec(x)
}

fn solve_upper_tri_into<R>(u: &CsMat<R>, diag: &Vec<&R>, b: &mut Vec<R>, x: &mut Vec<R>)
where R: Ring, for<'x> &'x R: RingOps<R> {
    debug_assert!(x.iter().all(|x_i| 
        x_i.is_zero())
    ); 

    let n = u.rows();
    for j in (0..n).rev() {
        if b[j].is_zero() { continue }

        let u_jj_inv = diag[j].inv().unwrap();
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

fn diag<'a, R>(u: &'a CsMat<R>) -> Vec<&'a R> {
    let n = u.rows();
    let indptr = u.indptr();
    let data = u.data();

    (0..n).map( |i| {
        let p = indptr.index(i + 1);
        &data[p - 1]
    }).collect_vec()
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