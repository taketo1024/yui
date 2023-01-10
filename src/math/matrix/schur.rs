use sprs::{CsMat, CsVec};

use crate::math::traits::{Ring, RingOps};
use super::sparse::CsMatExt;
use super::triang::inv_upper_tri;

pub fn schur_partial_upper_triang<R>(a: CsMat<R>, r: usize) -> Schur<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert!(a.is_csc());

    let [u, b, c, d] = a.divide4(r, r);
    let uinv = inv_upper_tri(&u);
    let s = Schur::new(uinv, b, c, d);

    s
}

// A = [a b]  ~>  s = d - c a⁻¹ b.
//     [c d]
//
// [a⁻¹     ] [a b] [id -a⁻¹b] = [id  ]
// [-ca⁻¹ id] [c d] [      id]   [   s]

pub struct Schur<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    ainv: CsMat<R>,
    b: CsMat<R>,
    c: CsMat<R>,
    d: CsMat<R>
}

impl<R> Schur<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(ainv: CsMat<R>, b: CsMat<R>, c: CsMat<R>, d: CsMat<R>) -> Self {
        assert_eq!(ainv.rows(), b.rows());
        assert_eq!(c.rows(), d.rows());
        assert_eq!(ainv.cols(), c.cols());
        assert_eq!(b.cols(), d.cols());
        Self { ainv, b, c, d }
    }

    pub fn complement(&self) -> CsMat<R> {
        let (ainv, b, c, d) = self.blocks();
        d - &( &(c * ainv) * b)
    }

    fn blocks(&self) -> (&CsMat<R>, &CsMat<R>, &CsMat<R>, &CsMat<R>) {
        (&self.ainv, &self.b, &self.c, &self.d)
    }
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

        let s = Schur::new(ainv, b, c, d);
        assert_eq!(s.complement(), CsMat::csc_from_vec((3,3), vec![27,29,26,1,-2,-5,-20,-24,-26]));
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
        assert_eq!(s.complement(), CsMat::csc_from_vec((2,3), vec![5,12,-14,36,45,-60]));
    }
}