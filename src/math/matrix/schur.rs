use sprs::CsMat;

use crate::math::traits::{Ring, RingOps};
use super::CsMatElem;
use super::sparse::CsMatExt;
use super::triang::inv_upper_tri;

// a = [u b]  ~>  s = d - c u⁻¹ b.
//     [c d]
//
// [u⁻¹     ] [u b] [id -u⁻¹b] = [id  ]
// [-cu⁻¹ id] [c d] [      id]   [   s]

pub fn schur_partial_upper_triang<R>(a: CsMat<R>, r: usize) -> CsMat<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> {
    assert!(a.is_csc());

    let [u, b, c, d] = a.divide4(r, r);
    let uinv = inv_upper_tri(&u);
    let s = schur(&uinv, &b, &c, &d);

    s
}

fn schur<R>(ainv: &CsMat<R>, b: &CsMat<R>, c: &CsMat<R>, d: &CsMat<R>) -> CsMat<R>
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

        let s = schur(&ainv, &b, &c, &d);
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
        let s = schur_partial_upper_triang(a, 3);
        assert_eq!(s, CsMat::csc_from_vec((2,3), vec![5,12,-14,36,45,-60]));
    }
}