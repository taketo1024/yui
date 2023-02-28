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
    pub fn from_partial_lower(a: SpMatView<R>, r: usize) -> Self {
        assert!(r <= a.rows());
        assert!(r <= a.cols());

        let pinv = Self::compute_pinv(&a, r);
        let compl = Self::compute_compl(&a, r, &pinv);

        Self { r, compl, pinv }
    }

    // p⁻¹ = [a    ] : lower triangular
    //       [c  id]
    fn compute_pinv(a: &SpMatView<R>, r: usize) -> SpMat<R> {
        let m = a.rows();
        SpMat::generate((m, m), |set| { 
            // [a; c]
            for (i, j, x) in a.submat_cols(0..r).iter() {
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
        let bd = a.submat_cols(r..n).to_owned();

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

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn schur_compl() {
        let a = SpMat::from_vec((6, 5), vec![
            1, 0, 0, 1, 3,
            2,-1, 0, 2, 2,
            3, 2, 1, 0, 3,
            4, 2, 4,-3, 0,
            5, 3, 5, 2, 2,
            6, 2,-3, 1, 8
        ]);
        let s = SchurLT::from_partial_lower(a.view(), 3);
        assert_eq!(s.complement(), &SpMat::from_vec((3,2), vec![5,36,12,45,-14,-60]));
    }

    #[test]
    fn trans_vec() {
        let a = SpMat::from_vec((6, 5), vec![
            1, 0, 0, 1, 3,
            2,-1, 0, 2, 2,
            3, 2, 1, 0, 3,
            4, 2, 4,-3, 0,
            5, 3, 5, 2, 2,
            6, 2,-3, 1, 8
        ]);
        let s = SchurLT::from_partial_lower(a.view(), 3);

        let v = SpVec::from(vec![1,2,0,-3,2,1]);
        let w = s.trans_vec(v);
        assert_eq!(w, SpVec::from(vec![5,12,-14]));

        let v = SpVec::from(vec![3,2,3,0,2,8]);
        let w = s.trans_vec(v);
        assert_eq!(w, SpVec::from(vec![36,45,-60]));
    }
}