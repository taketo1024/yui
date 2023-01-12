use sprs::{CsMat, CsVec};

use crate::math::matrix::sparse::CsVecExt;
use crate::math::traits::{Ring, RingOps};
use super::sparse::CsMatExt;
use super::triang::inv_upper_tri;

pub fn schur_partial_upper_triang<R>(a: CsMat<R>, r: usize) -> Schur<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert!(a.is_csc());

    let [u, b, c, d] = a.divide4(r, r);
    let uinv = inv_upper_tri(&u);
    let s = Schur::new(u, uinv, b, c, d);

    s
}

// A = [a b]  ~>  s = d - c a⁻¹ b.
//     [c d]
//
// [a⁻¹     ] [a b] [id -a⁻¹b] = [id  ]
// [-ca⁻¹ id] [c d] [      id]   [   s]

pub struct Schur<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    a: CsMat<R>,
    ainv: CsMat<R>,
    b: CsMat<R>,
    c: CsMat<R>,
    d: CsMat<R>
}

impl<R> Schur<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(a: CsMat<R>, ainv: CsMat<R>, b: CsMat<R>, c: CsMat<R>, d: CsMat<R>) -> Self {
        let (r, m, n) = (a.rows(), c.rows(), b.cols());
        assert_eq!(a.shape(), ainv.shape());
        assert_eq!(a.cols(), r);
        assert_eq!(b.rows(), r);
        assert_eq!(c.cols(), r);
        assert_eq!(d.rows(), m);
        assert_eq!(d.cols(), n);

        Self { a, ainv, b, c, d }
    }

    pub fn rows(&self) -> usize {
        self.a.rows() + self.c.rows()
    }

    pub fn cols(&self) -> usize {
        self.a.cols() + self.b.cols()
    }

    pub fn complement(&self) -> CsMat<R> {
        let (ainv, b, c, d) = (&self.ainv, &self.b, &self.c, &self.d);
        d - &( &(c * ainv) * b)
    }

    // returns: [-ca⁻¹ id][v1] = v2 - ca⁻¹v1
    //                    [v2] 
    pub fn trans_vec(&self, v: CsVec<R>) -> CsVec<R> { 
        assert_eq!(v.dim(), self.rows());
        
        let (ainv, c) = (&self.ainv, &self.c);
        let r = self.ainv.rows();
        let (v1, v2) = v.divide2(r);

        &v2 - &(c * &(ainv * &v1))
    }

    // p = [a⁻¹     ]
    //     [-ca⁻¹ id]
    pub fn p(&self) -> CsMat<R> { 
        let (r, m) = (self.ainv.rows(), self.rows());
        let (ainv, c) = (&self.ainv, &self.c);

        let zero = |shape| CsMat::<R>::zero(shape);
        let id = |n| CsMat::<R>::id(n);

        CsMat::combine4([
            ainv, 
            &zero((r, m - r)), 
            &( &zero((m - r, r)) - &(c * ainv) ), // MEMO: unary - not supported.
            &id(m - r)
        ])
    }

    // p⁻¹ = [a    ]
    //       [c  id]
    pub fn pinv(&self) -> CsMat<R> { 
        let (r, m) = (self.a.rows(), self.rows());
        let (a, c) = (&self.a, &self.c);

        let zero = |shape| CsMat::<R>::zero(shape);
        let id = |n| CsMat::<R>::id(n);

        CsMat::combine4([
            a, 
            &zero((r, m - r)), 
            c,
            &id(m - r)
        ])
    }

    // q = [id -a⁻¹b]
    //     [      id]
    pub fn q(&self) -> CsMat<R> { 
        let (r, n) = (self.a.rows(), self.cols());
        let (ainv, b) = (&self.ainv, &self.b);

        let zero = |shape| CsMat::<R>::zero(shape);
        let id = |n| CsMat::<R>::id(n);

        CsMat::combine4([
            &id(r), 
            &( &zero((r, n - r)) - &(ainv * b) ), // MEMO: unary - not supported., 
            &zero((n - r, r)),
            &id(n - r)
        ])
    }

    // q⁻¹ = [id a⁻¹b]
    //       [     id]
    pub fn qinv(&self) -> CsMat<R> { 
        let (r, n) = (self.a.rows(), self.cols());
        let (ainv, b) = (&self.ainv, &self.b);

        let zero = |shape| CsMat::<R>::zero(shape);
        let id = |n| CsMat::<R>::id(n);

        CsMat::combine4([
            &id(r), 
            &(ainv * b),
            &zero((n - r, r)),
            &id(n - r)
        ])
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn schur_compl_block() {
        let x = CsMat::csc_from_vec((5,5), vec![
            3,5,1,2,3,
            1,2,4,5,6,
            1,0,9,8,2,
            1,2,5,3,1,
            3,7,3,4,7
        ]);
        let [a,b,c,d] = x.divide4(2,2);
        let ainv = CsMat::csc_from_vec((2, 2), vec![2,-5,-1,3]);

        let s = Schur::new(a, ainv, b, c, d);
        assert_eq!(
            s.complement().to_dense(), 
            CsMat::csc_from_vec((3,3), vec![27,29,26,1,-2,-5,-20,-24,-26]).to_dense()
        );
    }

    #[test]
    fn trans_mats() {
        let x = CsMat::csc_from_vec((5,5), vec![
            3,5,1,2,3,
            1,2,4,5,6,
            1,0,9,8,2,
            1,2,5,3,1,
            3,7,3,4,7
        ]);
        let [a,b,c,d] = x.divide4(2,2);
        let ainv = CsMat::csc_from_vec((2, 2), vec![2,-5,-1,3]);

        let s = Schur::new(a, ainv, b, c, d);

        let id = |n| CsMat::<i32>::id(n);
        let (p, pinv) = (s.p(), s.pinv());
        let (q, qinv) = (s.q(), s.qinv());

        assert_eq!((&p * &pinv).to_dense(), id(5).to_dense());
        assert_eq!((&q * &qinv).to_dense(), id(5).to_dense());

        let y = &(&p * &x) * &q;
        let [e,f,g,h] = y.divide4(2,2);

        assert_eq!(e.to_dense(), id(2).to_dense());
        assert!(f.is_zero());
        assert!(g.is_zero());
        assert_eq!(h.to_dense(), s.complement().to_dense());
    }

    #[test]
    fn schur_compl() {
        let a = CsMat::csc_from_vec((5, 6), vec![
            1, 2, 3, 4, 5, 6,
            0,-1, 2, 2, 3, 2,
            0, 0, 1, 4, 5,-3,
            1, 2, 0,-3, 2, 1,
            3, 2, 3, 0, 2, 8
        ]);
        let s = schur_partial_upper_triang(a.clone(), 3);

        assert_eq!(s.complement(), CsMat::csc_from_vec((2,3), vec![5,12,-14,36,45,-60]));

        let id = |n| CsMat::<i32>::id(n);
        let (p, pinv) = (s.p(), s.pinv());
        let (q, qinv) = (s.q(), s.qinv());

        assert_eq!((&p * &pinv).to_dense(), id(5).to_dense());
        assert_eq!((&q * &qinv).to_dense(), id(6).to_dense());

        let y = &(&p * &a) * &q;
        let [e,f,g,h] = y.divide4(3,3);

        assert_eq!(e.to_dense(), id(3).to_dense());
        assert!(f.is_zero());
        assert!(g.is_zero());
        assert_eq!(h.to_dense(), s.complement().to_dense());
    }

    #[test]
    fn trans_vec() {
        let a = CsMat::csc_from_vec((5, 6), vec![
            1, 2, 3, 4, 5, 6,
            0,-1, 2, 2, 3, 2,
            0, 0, 1, 4, 5,-3,
            1, 2, 0,-3, 2, 1,
            3, 2, 3, 0, 2, 8
        ]);

        let s = schur_partial_upper_triang(a, 3);
        let v = CsVec::new(5, (0..5).collect(), (0..5).collect());
        let w = s.trans_vec(v);
        assert_eq!(w, CsVec::new(2, vec![0,1], vec![9, 28]));
    }
}