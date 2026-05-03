use std::ops::AddAssign;

use log::debug;
use nalgebra::Scalar;
use num_traits::{One, Zero};
use yui_core::{Ring, RingOps};
use super::*;
use super::triang::{TriangularType, solve_triangular_left, solve_triangular_with};

//                [a  b]
//                [c  d]
//            X ----------> Y
//  [1 -a⁻¹b] ^             | [1      ]
//  [     1 ] |   [a   ]    | [-ca⁻¹ 1]
//            |   [   s]    V
//            X ----------> Y
//       [0]  ^             | 
//       [1]  |             | [0  1]
//            |      s      V
//            X'----------> Y'
//
// s = d - c a⁻¹ b

pub struct Schur<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    s: SpMat<R>,
    t_src: Option<Trans<R>>,
    t_tgt: Option<Trans<R>>,
}

impl<R> Schur<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_partial_triangular(t: TriangularType, abcd: &SpMat<R>, r: usize, with_trans: bool) -> Self {
        assert!(r <= abcd.nrows());
        assert!(r <= abcd.ncols());

        let (m, n) = abcd.shape();
        let [a, b, c, d] = abcd.divide4((r, r));
        let (m_d, n_b) = (m - r, n - r);

        debug!("schur: a{:?}, b{:?}, c{:?}, d{:?}", a.shape(), b.shape(), c.shape(), d.shape());

        // Fuse the triangular solve `a · x_j = b_j` with the Schur update `s_j = d_j - c · x_j`. 
        // When `with_trans` is false, `x_j` is dropped immediately, so `a⁻¹b` is never materialized.
        let (s, t_src) = if with_trans {
            let pairs = solve_triangular_with(t, &a, &b, |j, x_j| {
                let s_j = d.col_vec(j) - &c * &x_j;
                (x_j, s_j)
            });
            let (x_cols, s_cols): (Vec<_>, Vec<_>) = pairs.into_iter().unzip();
            let x = SpMat::from_col_vecs(r, x_cols);   // x = a⁻¹b
            let s = SpMat::from_col_vecs(m_d, s_cols);

            let t_src_f = proj_mat(n, n_b);
            let t_src_b = (-x).stack(&id_mat(n_b)); // [-x, 1]^T

            (s, Some(Trans::new(t_src_f, t_src_b)))
        } else {
            let s_cols = solve_triangular_with(t, &a, &b, |j, x| {
                d.col_vec(j) - &c * x
            });
            (SpMat::from_col_vecs(m_d, s_cols), None)
        };

        debug!("schur: {:?}", s.shape());

        let t_tgt = with_trans.then(|| {
            let mut t_tgt_f = -solve_triangular_left(t, &a, &c); // (-x)a = c
            t_tgt_f.extend_cols(id_mat(m_d)); // [-ca⁻¹, 1]
            let t_tgt_b = incl_mat(m, m_d);   // [0, 1]^T
            Trans::new(t_tgt_f, t_tgt_b)
        });

        Self { s, t_src, t_tgt }
    }

    pub fn complement(&self) -> &SpMat<R> {
        &self.s
    }

    pub fn trans_src(&self) -> Option<&Trans<R>> { 
        self.t_src.as_ref()
    }

    pub fn trans_tgt(&self) -> Option<&Trans<R>> { 
        self.t_tgt.as_ref()
    }

    pub fn disassemble(self) -> (SpMat<R>, Option<Trans<R>>, Option<Trans<R>>) {
        (self.s, self.t_src, self.t_tgt)
    }
}

fn id_mat<R: Scalar + One>(n: usize) -> SpMat<R> { 
    SpMat::<R>::id(n)
}

fn incl_mat<R: Scalar + One + Zero + AddAssign>(n: usize, k: usize) -> SpMat<R> {
    SpMat::from_entries((n, k), (0..k).map(|i| (n - k + i, i, R::one()))) // [0, 1]^T
}

fn proj_mat<R: Scalar + One + Zero + AddAssign>(n: usize, k: usize) -> SpMat<R> {
    SpMat::from_entries((k, n), (0..k).map(|i| (i, n - k + i, R::one()))) // [0, 1]
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn schur_lower() {
        let a = SpMat::from_dense_data((6, 5), [
            1, 0, 0, 1, 3,
            2,-1, 0, 2, 2,
            3, 2, 1, 0, 3,
            4, 2, 4,-3, 0,
            5, 3, 5, 2, 2,
            6, 2,-3, 1, 8
        ]);
        let sch = Schur::from_partial_triangular(TriangularType::Lower, &a, 3, false);
        let s = sch.complement();

        assert_eq!(s, &SpMat::from_dense_data((3,2), [
             5,  36, 
             12, 45,
            -14,-60
        ]));
        assert!(sch.trans_src().is_none());
        assert!(sch.trans_tgt().is_none());
    }

    #[test]
    fn schur_lower_with_trans() {
        let a = SpMat::from_dense_data((6, 5), [
            1, 0, 0, 1, 3,
            2,-1, 0, 2, 2,
            3, 2, 1, 0, 3,
            4, 2, 4,-3, 0,
            5, 3, 5, 2, 2,
            6, 2,-3, 1, 8
        ]);
        let sch = Schur::from_partial_triangular(TriangularType::Lower, &a, 3, true);
        let s = sch.complement();

        assert_eq!(s, &SpMat::from_dense_data((3,2), [
             5,  36, 
             12, 45,
            -14,-60
        ]));
        assert!(sch.trans_src().is_some());
        assert!(sch.trans_tgt().is_some());

        let t_in  = sch.trans_src().unwrap().backward_mat();
        let t_out = sch.trans_tgt().unwrap().forward_mat();

        assert_eq!(t_in, SpMat::from_dense_data((5,2), [
            -1, -3,
             0, -4,
             3, 14,
             1,  0,
             0,  1
        ]));
        
        assert_eq!(t_out, SpMat::from_dense_data((3,6), [
             20, -6, -4, 1, 0, 0,
             24, -7, -5, 0, 1, 0,
            -31,  8,  3, 0, 0, 1
        ]));

        assert_eq!(&(t_out * &a * t_in), s);
    }

    #[test]
    fn schur_upper() {
        let a = SpMat::from_dense_data((5, 6), [
            1, 2, 3, 4, 5, 6,
            0, -1, 2, 2, 3, 2,
            0, 0, 1, 4, 5, -3,
            1, 2, 0, -3, 2, 1,
            3, 2, 3, 0, 2, 8,
        ]);
        let sch = Schur::from_partial_triangular(TriangularType::Upper, &a, 3, false);
        let s = sch.complement();

        assert_eq!(s, &SpMat::from_dense_data((2, 3), [
            5, 12,-14,
            36,45,-60
        ]));
        assert!(sch.trans_src().is_none());
        assert!(sch.trans_tgt().is_none());
    }

    #[test]
    fn schur_upper_with_trans() {
        let a = SpMat::from_dense_data((5, 6), [
            1, 2, 3, 4, 5, 6,
            0, -1, 2, 2, 3, 2,
            0, 0, 1, 4, 5, -3,
            1, 2, 0, -3, 2, 1,
            3, 2, 3, 0, 2, 8,
        ]);
        let sch = Schur::from_partial_triangular(TriangularType::Upper, &a, 3, true);
        let s = sch.complement();

        assert_eq!(s, &SpMat::from_dense_data((2, 3), [
            5, 12,-14,
            36,45,-60
        ]));
        assert!(sch.trans_src().is_some());
        assert!(sch.trans_tgt().is_some());

        let t_in  = sch.trans_src().unwrap().backward_mat();
        let t_out = sch.trans_tgt().unwrap().forward_mat();

        assert_eq!(t_in,  SpMat::from_dense_data((6,3), [
            20, 24, -31,
            -6, -7,   8,
            -4, -5,   3,
             1,  0,   0,
             0,  1,   0,
             0,  0,   1
        ]));

        assert_eq!(t_out, SpMat::from_dense_data((2, 5), [
            -1,  0,  3, 1, 0,
            -3, -4, 14, 0, 1
        ]));

        assert_eq!(&(t_out * &a * t_in), s);
    }
}