use std::ops::AddAssign;

use log::debug;
use nalgebra::Scalar;
use num_traits::{One, Zero};
use sprs::PermOwned;
use yui_core::{Ring, RingOps};
use crate::sparse::pivot::{PivotType, split_by_pqr};

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
    pub fn from_pivots(
        a: &SpMat<R>,
        t: PivotType,
        p: &PermOwned,
        q: &PermOwned,
        r: usize,
        with_trans_src: bool,
        with_trans_tgt: bool,
    ) -> Self {
        let (m, n) = a.shape();
        assert!(r <= m);
        assert!(r <= n);

        let t = if t == PivotType::Rows { TriangularType::Upper } else { TriangularType::Lower };
        let blocks = split_by_pqr(a, p, q, r);
        Self::from_blocks(t, blocks, with_trans_src, with_trans_tgt)
    }

    pub(crate) fn from_partial_triangular(
        t: TriangularType,
        a: SpMat<R>,
        r: usize,
        with_trans_src: bool,
        with_trans_tgt: bool,
    ) -> Self {
        let (m, n) = a.shape();
        assert!(r <= m);
        assert!(r <= n);

        let blocks = a.divide_into_blocks((r, r));
        Self::from_blocks(t, blocks, with_trans_src, with_trans_tgt)
    }

    pub(crate) fn from_blocks(
        t: TriangularType,
        blocks: [SpMat<R>; 4],
        with_trans_src: bool,
        with_trans_tgt: bool,
    ) -> Self {
        let [a, b, c, d] = blocks;
        assert!(a.is_square());

        let r = a.nrows();
        let (m_d, n_b) = (d.nrows(), b.ncols());

        debug!("compute schur: a{:?}, r: {r}", (m_d + r, n_b + r));

        // Compute `s` via one of three fused paths, picking whichever matches the
        // requested transforms — never materializing a `(m_d × n_b)` matmul:
        //   - with_trans_src:    right-solve fusion → `s` and `a⁻¹b` together.
        //   - with_trans_tgt only: left-solve fusion (transposed view) → `s` and `c·a⁻¹` together.
        //   - neither:           right-solve streaming, `s` only.
        let pairs = solve_triangular_with(t, &a, &b, |j, x_j| {
            let s_j = d.col_vec(j) - &c * &x_j;
            let x_j = (with_trans_src).then_some(x_j);
            (s_j, x_j)
        });
        
        let (s_cols, x_cols): (Vec<_>, Vec<_>) = pairs.into_iter().unzip();
        let s = SpMat::from_col_vecs(m_d, s_cols);

        let t_src = with_trans_src.then(|| {
            let x = SpMat::from_col_vecs(r, x_cols.into_iter().map(|x| x.unwrap())); // x = a⁻¹b
            let f = proj_mat(r + n_b, n_b);
            let b = (-x).stack(&id_mat(n_b)); // [-a⁻¹b, 1]^T
            Trans::new(f, b)
        });

        let t_tgt = with_trans_tgt.then(|| {
            let mut f = -solve_triangular_left(t, &a, &c); // f = -ca⁻¹
            f.extend_cols(id_mat(m_d)); // [-c·a⁻¹, 1]
            let b = incl_mat(r + m_d, m_d); // [0, 1]^T
            Trans::new(f, b)
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

    /// Returns `a⁻¹b` if it was retained (i.e. `with_trans_src=true`).
    pub fn ainvb(&self) -> Option<SpMat<R>> {
        self.t_src.as_ref().map(|trans| {
            // backward_mat = [-a⁻¹b ; I], shape (r + n_b, n_b); top r rows are -a⁻¹b.
            let bb = trans.backward_mat();
            let r = bb.nrows() - bb.ncols();
            -bb.submat_rows(0..r)
        })
    }

    /// Returns `c·a⁻¹` if it was retained (i.e. `with_trans_tgt=true`).
    pub fn ca_inv(&self) -> Option<SpMat<R>> {
        self.t_tgt.as_ref().map(|trans| {
            // forward_mat = [-c·a⁻¹, I], shape (m_d, r + m_d); first r cols are -c·a⁻¹.
            let f = trans.forward_mat();
            let r = f.ncols() - f.nrows();
            -f.submat_cols(0..r)
        })
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
        let sch = Schur::from_partial_triangular(TriangularType::Lower, a, 3, false, false);
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
        let sch = Schur::from_partial_triangular(TriangularType::Lower, a.clone(), 3, true, true);
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
        let sch = Schur::from_partial_triangular(TriangularType::Upper, a, 3, false, false);
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
        let sch = Schur::from_partial_triangular(TriangularType::Upper, a.clone(), 3, true, true);
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