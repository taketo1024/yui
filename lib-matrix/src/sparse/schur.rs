use std::ops::AddAssign;

use log::debug;
use nalgebra::Scalar;
use num_traits::{One, Zero};
use sprs::PermOwned;
use yui_core::{Ring, RingOps};
use crate::Perm;
use crate::sparse::pivot::PivotType;

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
    col_mult: Option<SpMat<R>>, // a⁻¹·b — column-elimination multiplier
    row_mult: Option<SpMat<R>>, // c·a⁻¹ — row-elimination multiplier
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
        let (p_new, q_new) = (Perm::new(p.vec()), Perm::new(q.vec()));
        let [a0, a1, a2, a3] = a.permute_and_split(&p_new, &q_new, r);
        Self::from_blocks(t, [&a0, &a1, &a2, &a3], with_trans_src, with_trans_tgt)
    }

    pub(crate) fn from_blocks(
        t: TriangularType,
        blocks: [&SpMat<R>; 4],
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
        let pairs = solve_triangular_with(t, a, b, |j, x_j| {
            let s_j = d.col_vec(j) - c * &x_j;
            let x_j = (with_trans_src).then_some(x_j);
            (s_j, x_j)
        });

        let (s_cols, x_cols): (Vec<_>, Vec<_>) = pairs.into_iter().unzip();
        let s = SpMat::from_col_vecs(m_d, s_cols);

        let col_mult = with_trans_src.then(|| {
            SpMat::from_col_vecs(r, x_cols.into_iter().map(|x| x.unwrap())) // a⁻¹b
        });

        let row_mult = with_trans_tgt.then(|| {
            solve_triangular_left(t, a, c) // c·a⁻¹
        });

        Self { s, col_mult, row_mult }
    }

    pub fn complement(&self) -> &SpMat<R> {
        &self.s
    }

    /// Returns `a⁻¹·b` if it was retained (i.e. `with_trans_src=true`).
    pub fn col_mult(&self) -> Option<&SpMat<R>> {
        self.col_mult.as_ref()
    }

    /// Returns `c·a⁻¹` if it was retained (i.e. `with_trans_tgt=true`).
    pub fn row_mult(&self) -> Option<&SpMat<R>> {
        self.row_mult.as_ref()
    }

    pub fn trans_src(&self) -> Option<Trans<R>> {
        self.col_mult.as_ref().map(|x| {
            let (r, n_b) = (x.nrows(), x.ncols());
            let f = proj_mat(r + n_b, n_b);
            let b = SpMat::stack(-x, id_mat(n_b)); // [-a⁻¹b ; I]
            Trans::new(f, b)
        })
    }

    pub fn trans_tgt(&self) -> Option<Trans<R>> {
        self.row_mult.as_ref().map(|y| {
            let (m_d, r) = (y.nrows(), y.ncols());
            let f = SpMat::concat(-y, id_mat(m_d)); // [-c·a⁻¹, I]
            let b = incl_mat(r + m_d, m_d);
            Trans::new(f, b)
        })
    }

    pub fn disassemble(self) -> (SpMat<R>, Option<SpMat<R>>, Option<SpMat<R>>) {
        (self.s, self.col_mult, self.row_mult)
    }

    pub fn into_s(self) -> SpMat<R> {
        self.s
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
        let sch = Schur::from_pivots(&a, PivotType::Cols, &PermOwned::identity(6), &PermOwned::identity(5), 3, false, false);
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
        let sch = Schur::from_pivots(&a, PivotType::Cols, &PermOwned::identity(6), &PermOwned::identity(5), 3, true, true);
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
        let sch = Schur::from_pivots(&a, PivotType::Rows, &PermOwned::identity(5), &PermOwned::identity(6), 3, false, false);
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
        let sch = Schur::from_pivots(&a, PivotType::Rows, &PermOwned::identity(5), &PermOwned::identity(6), 3, true, true);
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