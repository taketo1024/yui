use once_cell::sync::OnceCell;
use sprs::CsVec;

use crate::math::matrix::sparse::CsVecExt;
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
    ac: SpMat<R>, 
    bd: SpMat<R>, 
    pinv: OnceCell<SpMat<R>>
}

impl<R> SchurLT<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_partial_lower(a: SpMatView<R>, r: usize) -> Self {
        assert!(r <= a.rows());
        assert!(r <= a.cols());

        debug_assert!(a.iter().all(|(i, j, a)| { 
            (j < r && (i >= j || a.is_zero())) || (j >= r) 
        }));

        // TODO improve performance!
        
        let ac = a.submatrix_cols(0 .. r).to_owned();
        let bd = a.submatrix_cols(r .. a.cols()).to_owned();
        let pinv = OnceCell::new();
        Self { ac, bd, pinv }
    }

    // p⁻¹ = [a    ] : lower triangular
    //       [c  id]
    pub fn pinv(&self) -> &SpMat<R> { 
        self.pinv.get_or_init(|| { 
            let zero = |n, m| SpMat::<R>::zero((n, m));
            let id = |n| SpMat::<R>::id(n);
    
            let (n, r) = self.ac.shape();
            let right = zero(r, n - r).stack(&id(n - r));
            let pinv = self.ac.concat(&right);

            pinv
        })
    }

    // s = d - c a⁻¹ b is given by the x2 of 
    //
    //   [a    ][x1] = [b]
    //   [c  id][x2]   [d].
    //
    pub fn complement(&self) -> SpMat<R> {
        let (m, r) = self.ac.shape();
        let k = self.bd.cols();
        let pinv = self.pinv();

        SpMat::generate((m - r, k), |init| { 
            solve_triangular_with(TriangularType::Lower, pinv, &self.bd, |i, j, x|
                if i >= r { 
                    init(i - r, j, x)
                }
            );
        })
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
    pub fn trans_vec(&self, v: CsVec<R>) -> CsVec<R> { 
        let (m, r) = self.ac.shape();
        let pinv = self.pinv();
        let x = solve_triangular_vec(TriangularType::Lower, pinv, &v);
        x.subvec(r..m)
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