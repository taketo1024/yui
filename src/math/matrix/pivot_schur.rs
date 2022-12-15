use sprs::CsMat;

use crate::math::traits::{Ring, RingOps};
use super::CsMatElem;
use super::sparse::CsMatExt;
use super::triang::inv_upper_tri;
use super::pivot::{find_pivots, perms_by_pivots};

pub fn pivot_schur<R>(a: CsMat<R>) -> CsMat<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> {
    pivot_schur_into(a.clone())
}

pub fn pivot_schur_into<R>(a: CsMat<R>) -> CsMat<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> {
    pivot_schur_into_step(a, 1)
}

fn pivot_schur_into_step<R>(a: CsMat<R>, step: usize) -> CsMat<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> {
    let pivs = find_pivots(&a);
    let r = pivs.len();

    if r == 0 { 
        return a
    }

    let (p, q) = perms_by_pivots(&a, &pivs);
    let a = a.permute(p.view(), q.view());
    let s = schur_compl(a, r);

    if s.is_zero() { 
        s
    } else {
        pivot_schur_into_step(s, step + 1)
    }
}

fn schur_compl<R>(a: CsMat<R>, r: usize) -> CsMat<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> {
    assert!(a.is_csc());

    let (m, n) = a.shape();

    // a = [u b]  ~>  s = d - c u⁻¹ b.
    //     [c d]
    //
    // [u⁻¹     ] [u b] [id -u⁻¹b] = [id  ]
    // [-cu⁻¹ id] [c d] [      id]   [   s]

    let (a1, a2) = (a.slice_outer(0..r), a.slice_outer(r..n));

    let u = a1.submatrix(0..r, 0..r);
    let b = a2.submatrix(0..r, 0..n-r);
    let c = a1.submatrix(r..m, 0..r);
    let d = a2.submatrix(r..m, 0..n-r);

    let uinv = inv_upper_tri(&u);
    let s = schur_compl_block(&uinv, &b, &c, &d);

    s
}

fn schur_compl_block<R>(ainv: &CsMat<R>, b: &CsMat<R>, c: &CsMat<R>, d: &CsMat<R>) -> CsMat<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> {
    assert_eq!(ainv.rows(), b.rows());
    assert_eq!(c.rows(), d.rows());
    assert_eq!(ainv.cols(), c.cols());
    assert_eq!(b.cols(), d.cols());
    d - &( &(c * ainv) * b)
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn _schur_compl_block() {
        let ainv = CsMat::csc_from_vec((2, 2), vec![2,-5,-1,3]);
        let b = CsMat::csc_from_vec((2, 3), vec![1,2,3,4,5,6]);
        let c = CsMat::csc_from_vec((3, 2), vec![1,0,1,2,3,7]);
        let d = CsMat::csc_from_vec((3, 3), vec![9,8,2,5,3,1,3,4,7]);

        let s = schur_compl_block(&ainv, &b, &c, &d);
        assert_eq!(s, CsMat::csc_from_vec((3,3), vec![27,29,26,1,-2,-5,-20,-24,-26]));
    }

    #[test]
    fn _schur_compl() {
        let a = CsMat::csc_from_vec((5, 6), vec![
            1, 2, 3, 4, 5, 6,
            0,-1, 2, 2, 3, 2,
            0, 0, 1, 4, 5,-3,
            1, 2, 0,-3, 2, 1,
            3, 2, 3, 0, 2, 8
        ]);
        let s = schur_compl(a, 3);
        assert_eq!(s, CsMat::csc_from_vec((3,3), vec![5,36,12,45,-14,-60]));
    }
}