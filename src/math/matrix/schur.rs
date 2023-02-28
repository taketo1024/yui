use crate::math::matrix::sp_vec::SpVec;
use crate::math::matrix::triang::TriangularType;
use crate::math::traits::{Ring, RingOps};
use super::sparse::*;
use super::triang::{solve_triangular_vec, solve_triangular_with};

// A = [a b]  ~>  s = d - c a⁻¹ b.
//     [c d]
//
// [a⁻¹     ] [a b] [id -a⁻¹b] = [id  ]
// [-ca⁻¹ id] [c d] [      id]   [   s]

pub struct SchurLT<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    r: usize,
    compl: SpMat<R>,
    pinv: SpMat<R>
}

impl<'a, R> SchurLT<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_partial_lower(a: &SpMatView<R>, r: usize) -> Self {
        assert!(r <= a.rows());
        assert!(r <= a.cols());

        let pinv = Self::compute_pinv(a, r);
        let compl = Self::compute_compl(a, r, &pinv);

        Self { r, compl, pinv }
    }

    // p⁻¹ = [a    ] : lower triangular
    //       [c  id]
    fn compute_pinv(a: &SpMatView<R>, r: usize) -> SpMat<R> {
        let m = a.rows();
        SpMat::generate((m, m), |set| { 
            // [a; c]
            for (i, j, x) in a.submatrix_cols(0..r).iter() {
                set(i, j, x.clone());
            }
            // [0; id]
            for i in r..m { 
                set(i, i, R::one());
            }
        })
    }

    // s = d - c a⁻¹ b is given by the x2 of 
    //
    //   [a    ][x1] = [b]
    //   [c  id][x2]   [d].
    //
    fn compute_compl(a: &SpMatView<R>, r: usize, pinv: &SpMat<R>) -> SpMat<R> { 
        let (m, n) = a.shape();
        let bd = a.submatrix_cols(r..n).to_owned();

        SpMat::generate((m - r, n - r), |set| { 
            solve_triangular_with(TriangularType::Lower, pinv, &bd, |i, j, x|
                if i >= r { 
                    set(i - r, j, x)
                }
            );
        })
    }

    pub fn complement(&self) -> &SpMat<R> {
        &self.compl
    }

    pub fn complement_into(self) -> SpMat<R> {
        self.compl
    }

    // Given 
    //
    //   v = [v1], 
    //       [v2]
    //
    // returns: 
    //
    //   [-ca⁻¹ id][v1] = v2 - ca⁻¹v1
    //             [v2] 
    //
    // which is given by the x2 of
    //
    //   [a    ][x1] = [v1]
    //   [c  id][x2]   [v2].
    //
    pub fn trans_vec(&self, v: SpVec<R>) -> SpVec<R> { 
        let r = self.r;
        let m = r + self.compl.rows();

        assert_eq!(m, v.dim());

        let pinv = &self.pinv;
        let x = solve_triangular_vec(TriangularType::Lower, pinv, &v);
        let x2 = x.subvec(r..m).to_owned();

        x2
    }
}

/*
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
*/