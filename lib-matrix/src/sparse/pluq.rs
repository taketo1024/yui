// Sparse PLUQ decomposition & linear solver.
// Implemented with the help of Claude Code.

use log::debug;
use yui_core::{Ring, RingOps, Field, FieldOps};

use crate::{MatTrait, Perm};
use crate::dense::Mat;
use crate::dense::pluq::pluq as dense_pluq;
use super::SpMat;
use super::SpVec;
use super::pivot::{PivotFinderConfig, PivotType, find_pivots};
use super::schur::Schur;
use super::triang::{TriangularType, solve_triangular_vec};

/// Result of a sparse PLUQ decomposition.
///
/// Satisfies `p * A * q = l * u + s` where `s` is the
/// `(m - rank) × (n - rank)` Schur complement (bottom-right block).
pub struct SpPluq<R> {
    pub p: Perm,
    pub q: Perm,
    pub l: SpMat<R>,
    pub u: SpMat<R>,
    pub s: SpMat<R>,
}

impl<R> SpPluq<R> {
    /// Constructs a `PartialPluq` after asserting the shapes are mutually
    /// consistent: `l.n_cols() == u.n_rows() = r`, `l.n_rows() == p.dim() = m`,
    /// `u.n_cols() == q.dim() = n`, and `s.shape() == (m - r, n - r)`.
    pub fn new(p: Perm, q: Perm, l: SpMat<R>, u: SpMat<R>, s: SpMat<R>) -> Self {
        let r = l.n_cols();
        let m = l.n_rows();
        let n = u.n_cols();
        assert_eq!(r, u.n_rows(), "l.n_cols() must match u.n_rows()");
        assert_eq!(m, p.len(), "l.n_rows() must match p.len()");
        assert_eq!(n, q.len(), "u.n_cols() must match q.len()");
        assert_eq!(s.shape(), (m - r, n - r), "s shape must be (m - r, n - r)");
        Self { p, q, l, u, s }
    }

    pub fn rank(&self) -> usize { self.l.n_cols() }

    pub fn take_l(&mut self) -> SpMat<R> {
        std::mem::take(&mut self.l)
    }

    pub fn take_u(&mut self) -> SpMat<R> {
        std::mem::take(&mut self.u)
    }

    pub fn take_s(&mut self) -> SpMat<R> {
        std::mem::take(&mut self.s)
    }
}

/// Converts `a` into the trivial PLUQ whose Schur complement is `a` itself
impl<R> From<SpMat<R>> for SpPluq<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(a: SpMat<R>) -> Self {
        let (m, n) = a.shape();
        Self::new(
            Perm::id(m),
            Perm::id(n),
            SpMat::zero((m, 0)),
            SpMat::zero((0, n)),
            a,
        )
    }
}

/// Computes a partial PLUQ decomposition of `a` under the given pivot-finder
/// configuration.
///
/// Splits the permuted matrix into four blocks `[[a0|a1],[a2|a3]]` at row/col
/// `r`, then asks Schur to fuse the triangular solve with the Schur update.
///
/// Rows: top half `[a0|a1]` is `u` (upper triangular on the left). Schur produces
///   `l1 = a2·a0⁻¹` and `s = a3 - l1·a1`; final `l = [I_r; l1]`.
/// Cols: left half `[a0;a2]` is `l` (lower triangular on top). Schur produces
///   `u1 = a0⁻¹·a1` and `s = a3 - a2·u1`; final `u = [I_r | u1]`.
pub fn pre_pluq<R>(a: &SpMat<R>, config: PivotFinderConfig) -> SpPluq<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    debug!("compute sparse pre-pluq: {:?}", a.shape());

    let (m, n) = a.shape();
    let piv_type = config.piv_type;
    let (p, q, r) = find_pivots(a, config);

    if r == 0 {
        return SpPluq::new(Perm::id(m), Perm::id(n), SpMat::zero((m, 0)), SpMat::zero((0, n)), a.clone());
    }

    let [a0, a1, a2, a3] = a.permute_and_split(&p, &q, r);

    let (l, u, s) = match piv_type {
        PivotType::Rows => {
            let sch = Schur::from_blocks(TriangularType::Upper, [&a0, &a1, &a2, &a3], false, true);
            let (s, _, row_mult) = sch.disassemble();
            let l1 = row_mult.unwrap();
            let u = SpMat::h_stack(a0, a1);          // u = [a0 | a1]
            let l = SpMat::v_stack(SpMat::id(r), l1); // l = [I_r ; l1]
            (l, u, s)
        },
        PivotType::Cols => {
            let sch = Schur::from_blocks(TriangularType::Lower, [&a0, &a1, &a2, &a3], true, false);
            let (s, col_mult, _) = sch.disassemble();
            let u1 = col_mult.unwrap();
            let l = SpMat::v_stack(a0, a2);            // l = [a0 ; a2]
            let u = SpMat::h_stack(SpMat::id(r), u1); // u = [I_r | u1]
            (l, u, s)
        }
    };

    SpPluq::new(p, q, l, u, s)
}

/// Computes a full PLUQ decomposition of `a`.
///
/// Iterates `pre_pluq` on the current Schur complement while it keeps finding
/// sparse pivots; falls back to `dense_pluq_in` once the Schur complement is
/// non-zero but admits no further sparse pivots.
pub fn pluq<R>(a: &SpMat<R>, config: PivotFinderConfig) -> SpPluq<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    debug!("compute sparse pluq: {:?}", a.shape());

    let piv_type = config.piv_type;
    let mut pp = SpPluq::from(a.clone());

    while !pp.s.is_zero() {
        let pp_next = pre_pluq(&pp.s, config);
        if pp_next.rank() == 0 { break; }

        merge_pluq(&mut pp, pp_next);
    }

    if pp.s.is_zero() { return pp; }

    let pp_dense = dense_pluq_in(&pp.s, piv_type);
    merge_pluq(&mut pp, pp_dense);
    pp
}

fn dense_pluq_in<R>(s: &SpMat<R>, piv_type: PivotType) -> SpPluq<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let transpose = piv_type == PivotType::Rows;
    let (ms, ns) = s.shape();
    let (row_idx, col_idx, mat) = extract_dense(s, transpose);
    let (m0, n0) = (row_idx.len(), col_idx.len());

    let raw = dense_pluq(&mat);
    let dp = if transpose { raw.transpose() } else { raw };

    let p2 = extend_perm(ms, &row_idx, dp.p);
    let q2 = extend_perm(ns, &col_idx, dp.q);

    let mut l2 = SpMat::from(dp.l);
    l2.extend_by_zero(ms - m0, 0);

    let mut u2 = SpMat::from(dp.u);
    u2.extend_by_zero(0, ns - n0);

    let mut s2 = SpMat::from(dp.s);
    s2.extend_by_zero(ms - m0, ns - n0);

    SpPluq::new(p2, q2, l2, u2, s2)
}

// Extracts the compact dense submatrix of `s` using only its non-zero rows/cols.
// If `transpose` is true, returns S0^T (n0 × m0); otherwise returns S0 (m0 × n0).
// Also returns the sorted non-zero row/col indices of `s`.
fn extract_dense<R>(s: &SpMat<R>, transpose: bool) -> (Vec<usize>, Vec<usize>, Mat<R>)
where R: Ring, for<'x> &'x R: RingOps<R> {
    use std::collections::BTreeSet;

    let row_idx: Vec<usize> = s.iter_nz().map(|(i, _, _)| i).collect::<BTreeSet<_>>().into_iter().collect();
    let col_idx: Vec<usize> = s.iter_nz().map(|(_, j, _)| j).collect::<BTreeSet<_>>().into_iter().collect();
    let (m0, n0) = (row_idx.len(), col_idx.len());

    let row_perm = Perm::forward_indices(s.n_rows(), row_idx.iter().copied());
    let col_perm = Perm::forward_indices(s.n_cols(), col_idx.iter().copied());

    let shape = if transpose { (n0, m0) } else { (m0, n0) };
    let mut mat = Mat::zero(shape);

    for (i, j, v) in s.iter_nz() {
        let (ri, cj) = (row_perm.at(i), col_perm.at(j));
        if transpose {
            mat[(cj, ri)] = v.clone(); // S0^T[cj, ri] = S0[ri, cj]
        } else {
            mat[(ri, cj)] = v.clone();
        }
    }

    (row_idx, col_idx, mat)
}

// Merges two partial PLUQ decompositions. `pp1` has rank `r1` and shape (m, n);
// `pp2` is a partial PLUQ of `pp1.s` with rank `r2` and shape (m - r1, n - r1).
// Returns a partial PLUQ of the same matrix as `pp1` with rank `r1 + r2` and
// schur complement `pp2.s`.
fn merge_pluq<R>(pp1: &mut SpPluq<R>, pp2: SpPluq<R>)
where R: Ring, for<'x> &'x R: RingOps<R> {
    debug!("merge pluq: {} + {}", pp1.rank(), pp2.rank());

    let (m, n) = (pp1.l.n_rows(), pp1.u.n_cols());
    let r1 = pp1.rank();
    let r2 = pp2.rank();

    assert_eq!(pp2.l.n_rows(), m - r1);
    assert_eq!(pp2.u.n_cols(), n - r1);

    // MEMO: Even if r2 == 0, there could be non-trivial permutations
    // when R is not a field.

    pp1.l = {
        let [l0, l1] = pp1.take_l().v_split(r1);
        let l1 = l1.permute_rows(&pp2.p);
        let zero_tr = SpMat::zero((r1, r2));
        SpMat::block_combine([l0, zero_tr, l1, pp2.l])
    };

    pp1.u = {
        let [u0, u1] = pp1.take_u().h_split(r1);
        let u1 = u1.permute_cols(&pp2.q);
        let zero_bl = SpMat::zero((r2, r1));
        SpMat::block_combine([u0, u1, zero_bl, pp2.u])
    };

    pp1.s = pp2.s;
    pp1.p = merge_perm(&pp1.p, pp2.p);
    pp1.q = merge_perm(&pp1.q, pp2.q);
}

/// Solves `a * x = y` over a field using sparse PLUQ.
///
/// Returns `Some(x)` if a solution exists, `None` otherwise.
pub fn solve_pluq<R>(a: &SpMat<R>, y: &SpVec<R>) -> Option<SpVec<R>>
where R: Field, for<'x> &'x R: FieldOps<R> {
    debug!("solve pluq, a: {:?}", a.shape());

    assert_eq!(y.dim(), a.n_rows());

    if y.is_zero() {
        return Some(SpVec::zero(a.n_cols())); // y = 0 ⇒ x = 0 solves it — skip the factorization.
    }

    let pp = pluq(a, PivotFinderConfig {
        piv_type: PivotType::Rows,
        ..Default::default()
    });

    let y_dense = y.clone().into_dense();
    let yp = pp.p.apply_to(y_dense);
    let xq = solve_lu(&pp.l, &pp.u, &yp)?;
    let x = pp.q.apply_inv_to(xq);

    Some(SpVec::from(x))
}

// Solves `L * U * x = y` and returns `x` of length `n = u.n_cols()` with
// entries beyond `r = l.n_cols()` set to zero (free variables = 0).
//
// Requires the top r × r block of L to be unit lower triangular and the top
// r × r block of U to be invertible upper triangular.
//
// Returns `None` when `solve_l(l, y, true)` detects an inconsistent residual.
// When `l` is square (`l.n_rows() == r`) the residual is empty and the call
// always succeeds.
fn solve_lu<R>(l: &SpMat<R>, u: &SpMat<R>, y: &[R]) -> Option<Vec<R>>
where R: Field, for<'x> &'x R: FieldOps<R> {
    assert_eq!(l.n_cols(), u.n_rows());
    assert_eq!(y.len(), l.n_rows());

    let z = solve_l(l, y, true)?;
    let x = solve_u(u, &z);

    Some(x)
}

// Solves `l[0..r, 0..r] * z = y[0..r]` by forward substitution, where
// `r = l.n_cols()`. The top r × r block of L must be lower triangular with
// non-zero diagonal.
//
// If `check_consistency` is true and `r < y.len()`, also verifies the residual
// `y[r..] - l[r.., :] * z` is zero, returning `None` when it isn't. When `r ==
// y.len()` the residual is trivially empty so the check is skipped.
fn solve_l<R>(l: &SpMat<R>, y: &[R], check_consistency: bool) -> Option<Vec<R>>
where R: Field, for<'x> &'x R: FieldOps<R> {
    assert_eq!(l.n_rows(), y.len());
    let r = l.n_cols();

    let x = if r == y.len() {
        let y = SpVec::from(y.to_vec());
        solve_triangular_vec(TriangularType::Lower, l, &y).into_dense()
    } else {
        let l0 = l.submat(0..r, 0..r);
        let y0 = SpVec::from(y[..r].to_vec());
        let x = solve_triangular_vec(TriangularType::Lower, &l0, &y0).into_dense();

        if check_consistency && !is_consistent(l, y, &x) { 
            return None;
        }
        x
    };

    Some(x)
}

// check l * x == y
fn is_consistent<R>(l: &SpMat<R>, y: &[R], x: &[R]) -> bool
where R: Ring, for<'x> &'x R: RingOps<R> {
    is_consistent_upto(l, y, x, y.len())
}

fn is_consistent_upto<R>(l: &SpMat<R>, y: &[R], x: &[R], k: usize) -> bool
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert_eq!(l.n_rows(), y.len());
    assert_eq!(l.n_cols(), x.len());
    assert!(x.len() <= k && k <= y.len());

    let r = x.len();
    let mut res = y[r..k].to_vec();

    for (i, j, v) in l.iter_nz() {
        if r <= i && i < k {
            res[i - r] -= v * &x[j];
        }
    }

    res.iter().all(|v| v.is_zero())
}

// Solves `u[0..r, 0..r] * x[..r] = y` by back-substitution, where
// `r = u.n_rows()`, and returns `x` of length `n = u.n_cols()` with entries
// beyond `r` set to zero. The top r × r block of U must be upper triangular
// with non-zero diagonal.
fn solve_u<R>(u: &SpMat<R>, y: &[R]) -> Vec<R>
where R: Field, for<'x> &'x R: FieldOps<R> {
    let (r, n) = u.shape();
    assert_eq!(y.len(), r);
    assert!(n >= r);

    let mut x = if n == r {
        solve_triangular_vec(TriangularType::Upper, u, &SpVec::from(y.to_vec())).into_dense()
    } else {
        let u0 = u.submat(0..r, 0..r);
        solve_triangular_vec(TriangularType::Upper, &u0, &SpVec::from(y.to_vec())).into_dense()
    };

    x.resize(n, R::zero());
    x
}

/// Solves `a * x = y` over a field using an incremental sparse PLUQ.
///
/// Behaves like [`solve_pluq`] but designed for huge matrices: caps the initial
/// sparse pre-PLUQ at `max_piv` pivots, then incrementally processes the
/// remaining Schur complement `chunk` rows at a time. Returns `None` (without
/// completing the full PLUQ) as soon as a chunk reveals inconsistency.
pub fn solve_pluq_incr<R>(a: &SpMat<R>, y: &SpVec<R>, max_piv: usize, chunk: usize) -> Option<SpVec<R>>
where R: Field, for<'x> &'x R: FieldOps<R> {
    debug!("solve pluq (incremental), a: {:?}", a.shape());

    assert_eq!(y.dim(), a.n_rows());

    if y.is_zero() {
        return Some(SpVec::zero(a.n_cols())); // y = 0 ⇒ x = 0 solves it — skip the factorization.
    }

    let mut pp = pre_pluq(a, PivotFinderConfig {
        piv_type: PivotType::Rows,
        max_pivots: max_piv,
        ..Default::default()
    });
    let y_dense = y.clone().into_dense();
    let mut yp = pp.p.apply_to(y_dense);

    let mut step = 1;
    let total_step = (a.n_rows() - pp.rank()) / chunk + 1;

    while pp.s.n_rows() > 0 {
        debug!("(step {}/{})", step, total_step);
        debug!("  current rank: {}", pp.rank());

        let r_old = pp.rank();
        let (pp_next, r_next, c) = chunk_pluq(pp.take_s(), chunk);
        let p_next = pp_next.p.clone();
        
        merge_pluq(&mut pp, pp_next);

        // Apply the chunk's row perm to the tail of yp so it stays in sync with pp.l.
        let yp_tail = p_next.apply_to(yp[r_old..].to_vec());
        yp[r_old..].clone_from_slice(&yp_tail);

        // The top `k` rows of pp.s are zero rows (chunk's PLUQ leftover);
        // they demand `yp[r_new..r_new+k] == L[r_new..r_new+k, :] * z` for
        // consistency, regardless of future chunks.
        let k = c - r_next;
        let z = solve_l(&pp.l, &yp, false).unwrap();

        if !is_consistent_upto(&pp.l, &yp, &z, z.len() + k) {
            debug!("found inconsistency at step {}/{}.", step, total_step);
            return None;
        }

        trim_zero_rows(&mut pp, &mut yp, k);
        step += 1;
    }

    debug!("pluq complete, rank: {}", pp.rank());
    debug!("solve pluq..");

    let xq = solve_lu(&pp.l, &pp.u, &yp)?;
    let x = pp.q.apply_inv_to(xq);

    Some(SpVec::from(x))
}

// Takes the top `min(chunk_size, s.n_rows())` rows of `s`, runs `pluq` on them,
// and lifts the result to act on all of `s` via `extend_chunk_to_full`.
// Returns `(pp_chunk_full, r_chunk, c)`.
fn chunk_pluq<R>(s: SpMat<R>, chunk_size: usize) -> (SpPluq<R>, usize, usize)
where R: Ring, for<'x> &'x R: RingOps<R> {
    let c = chunk_size.min(s.n_rows());
    let [s_chunk, s_rest] = s.v_split(c);
    let pp_chunk = pluq(&s_chunk, PivotFinderConfig {
        piv_type: PivotType::Rows,
        ..Default::default()
    });
    let r_chunk = pp_chunk.rank();
    let pp_full = extend_chunk_to_full(pp_chunk, s_rest);
    (pp_full, r_chunk, c)
}

// Lifts a PLUQ of the top `c` rows of some matrix `s` (= `pp_chunk`) to a
// partial PLUQ acting on all rows of `s`, by absorbing the untouched bottom
// rows `s_rest` (shape (m_s - c, n_s)) into L (below the new pivots) and into
// the new schur complement.
//
// Sparse analog of `dense_pluq_in`, applied to a single chunk.
fn extend_chunk_to_full<R>(pp_chunk: SpPluq<R>, s_rest: SpMat<R>) -> SpPluq<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (c, n_s) = (pp_chunk.l.n_rows(), pp_chunk.u.n_cols());
    let r_chunk = pp_chunk.rank();
    let m_rest = s_rest.n_rows();
    let m_s = c + m_rest;

    assert_eq!(s_rest.n_cols(), n_s);
    assert_eq!(pp_chunk.s.shape(), (c - r_chunk, n_s - r_chunk));

    let s_rest_q = s_rest.permute_cols(&pp_chunk.q);
    let [s_rest_left, s_rest_right] = s_rest_q.h_split(r_chunk);
    let [u_top, u_right] = pp_chunk.u.clone().h_split(r_chunk);

    // Same Schur shape as pre_pluq's Rows branch: u_top (upper triangular) plays
    // the role of `a`, with `c = s_rest_left`, `b = u_right`, `d = s_rest_right`.
    let sch = Schur::from_blocks(
        TriangularType::Upper,
        [&u_top, &u_right, &s_rest_left, &s_rest_right],
        false, true
    );
    let (s_ext, _, row_mult) = sch.disassemble();
    let l_ext = row_mult.unwrap();

    let chunk_idx: Vec<usize> = (0..c).collect();
    let p = extend_perm(m_s, &chunk_idx, pp_chunk.p);
    let q = pp_chunk.q;
    let l = SpMat::v_stack(pp_chunk.l, l_ext);
    let u = pp_chunk.u;
    let s = SpMat::v_stack(pp_chunk.s, s_ext);

    SpPluq::new(p, q, l, u, s)
}

// Drops the top `k` zero rows of `pp.s` from `pp.l`, `pp.s`, `pp.p`, and `yp`.
//
// Caller is responsible for verifying via `is_consistent_upto` that the
// `k` rows being removed have zero residual; otherwise the resulting system
// would silently lose constraints.
fn trim_zero_rows<R>(pp: &mut SpPluq<R>, yp: &mut Vec<R>, k: usize)
where R: Ring, for<'x> &'x R: RingOps<R> {
    if k == 0 { return; }

    let r = pp.rank();
    let m = pp.l.n_rows();
    assert!(r + k <= m);
    assert_eq!(yp.len(), m);
    assert_eq!(pp.p.len(), m);

    // Drop rows [r..r+k] from pp.l: keep [0..r] and [r+k..m], shifted down.
    pp.l = pp.l.extract((m - k, r), |i, j| {
        if i < r {
            Some((i, j))
        } else if i < r + k {
            None
        } else {
            Some((i - k, j))
        }
    });

    // Drop the top k rows of pp.s.
    pp.s = pp.s.submat_rows(k..pp.s.n_rows());

    // Drop yp entries [r..r+k] to stay in sync with pp.l.
    yp.drain(r..r + k);

    // Drop pp.p entries that map to [r..r+k]; shift later positions down by k.
    let new_p_at: Vec<usize> = (0..m).filter_map(|i| {
        let pos = pp.p.at(i);
        if pos < r {
            Some(pos)
        } else if pos < r + k {
            None
        } else {
            Some(pos - k)
        }
    }).collect();
    pp.p = Perm::new(new_p_at);
}

// Composes perm1 with perm2: the first `r` positions stay, the rest are
// shifted by `r` and remapped by perm2 (where `r = perm1.len() - perm2.len()`).
fn merge_perm(perm1: &Perm, perm2: Perm) -> Perm {
    assert!(perm1.len() >= perm2.len());
    let r = perm1.len() - perm2.len();
    perm2.shift(r) * perm1
}

// Lifts a compact permutation (acting on compact_idx elements of [0..n]) to the
// full index space.  compact_idx[k] maps to compact_perm.at(k) (within [0..mr]);
// all other indices map to consecutive positions starting at mr (in sorted order).
fn extend_perm(n: usize, compact_idx: &[usize], compact_perm: Perm) -> Perm {
    let c = compact_idx.len();

    assert!(n >= c);
    assert_eq!(compact_perm.len(), c);

    let front = Perm::forward_indices(n, compact_idx.iter().copied());
    compact_perm.extend(n - c) * front
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::One;

    fn cfg(piv_type: PivotType) -> PivotFinderConfig {
        PivotFinderConfig { piv_type, ..Default::default() }
    }

    fn sample() -> SpMat<i32> {
        SpMat::from_row_major((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0,
        ])
    }

    // ---- pre_pluq ----

    fn check_pre_pluq_rows(a: &SpMat<i32>) {
        let pp = pre_pluq(a, cfg(PivotType::Rows));
        let (m, n) = a.shape();
        let r = pp.rank();

        assert_eq!(pp.l.shape(), (m, r));
        assert_eq!(pp.u.shape(), (r, n));
        assert_eq!(pp.s.shape(), (m - r, n - r));

        let paq = a.permute(&pp.p, &pp.q);
        let rem_full = SpMat::from_entries((m, n),
            pp.s.iter_nz().map(|(i, j, v)| (i + r, j + r, v.clone()))
        );
        assert_eq!(paq, &pp.l * &pp.u + &rem_full);

        let b = pp.u.clone().into_dense();
        for k in 0..r {
            assert!(b[(k, k)].is_one(), "u[{k},{k}] should be a pivot (=1)");
        }
        for j in 0..r {
            for i in j + 1..r {
                assert_eq!(b[(i, j)], 0, "u[{i},{j}] should be zero (below diagonal)");
            }
        }
    }

    fn check_pre_pluq_cols(a: &SpMat<i32>) {
        let pp = pre_pluq(a, cfg(PivotType::Cols));
        let (m, n) = a.shape();
        let r = pp.rank();

        assert_eq!(pp.l.shape(), (m, r));
        assert_eq!(pp.u.shape(), (r, n));
        assert_eq!(pp.s.shape(), (m - r, n - r));

        let paq = a.permute(&pp.p, &pp.q);
        let rem_full = SpMat::from_entries((m, n),
            pp.s.iter_nz().map(|(i, j, v)| (i + r, j + r, v.clone()))
        );
        assert_eq!(paq, &pp.l * &pp.u + &rem_full);

        let b = pp.l.clone().into_dense();
        for k in 0..r {
            assert!(b[(k, k)].is_one(), "l[{k},{k}] should be a pivot (=1)");
        }
        for i in 0..r {
            for j in i + 1..r {
                assert_eq!(b[(i, j)], 0, "l[{i},{j}] should be zero (above diagonal)");
            }
        }
    }

    #[test]
    fn test_pre_pluq_rows() { check_pre_pluq_rows(&sample()); }

    #[test]
    fn test_pre_pluq_cols() { check_pre_pluq_cols(&sample()); }

    #[test]
    fn test_pre_pluq_zero() {
        let a = SpMat::<i32>::zero((4, 5));
        let pp = pre_pluq(&a, cfg(PivotType::Rows));
        assert_eq!(pp.rank(), 0);
        assert_eq!(pp.l.shape(), (4, 0));
        assert_eq!(pp.u.shape(), (0, 5));
        assert_eq!(pp.s.shape(), (4, 5)); // (m-r, n-r) = (4, 5) when r=0
        assert_eq!(pp.s, a.permute(&pp.p, &pp.q));
    }

    #[test]
    fn test_pre_pluq_square_full_rank() {
        let a = SpMat::from_row_major((3, 3), [1, 0, 0, 0, 1, 0, 0, 0, 1]);
        let pp = pre_pluq(&a, cfg(PivotType::Rows));
        assert_eq!(pp.rank(), 3);
        assert_eq!(pp.s.shape(), (0, 0)); // full rank: Schur complement is empty
    }

    #[test]
    fn test_pre_pluq_rank_rows() {
        assert_eq!(pre_pluq(&sample(), cfg(PivotType::Rows)).rank(), 5);
    }

    #[test]
    fn test_pre_pluq_rank_cols() {
        assert_eq!(pre_pluq(&sample(), cfg(PivotType::Cols)).rank(), 6);
    }

    #[test]
    fn test_pre_pluq_rand_rows() { check_pre_pluq_rows(&SpMat::<i32>::rand((40, 60), 0.1)); }

    #[test]
    fn test_pre_pluq_rand_cols() { check_pre_pluq_cols(&SpMat::<i32>::rand((40, 60), 0.1)); }

    // ---- pluq ----

    fn check_pluq_rows(a: &SpMat<i32>) {
        let pp = pluq(a, cfg(PivotType::Rows));
        let (m, n) = a.shape();
        let r = pp.rank();

        assert_eq!(pp.l.shape(), (m, r));
        assert_eq!(pp.u.shape(), (r, n));
        assert_eq!(pp.s.shape(), (m - r, n - r));

        let paq = a.permute(&pp.p, &pp.q);
        let rem = SpMat::from_entries((m, n),
            pp.s.iter_nz().map(|(i, j, v)| (i + r, j + r, v.clone()))
        );
        assert_eq!(paq, &pp.l * &pp.u + &rem, "p*A*q != l*u + rest");

        // L is unit lower: 1s on diagonal
        let lb = pp.l.clone().into_dense();
        for k in 0..r {
            assert!(lb[(k, k)].is_one(), "l[{k},{k}] should be 1");
            for j in k + 1..r { assert_eq!(lb[(k, j)], 0, "l[{k},{j}] above diag"); }
        }
    }

    fn check_pluq_cols(a: &SpMat<i32>) {
        let pp = pluq(a, cfg(PivotType::Cols));
        let (m, n) = a.shape();
        let r = pp.rank();

        assert_eq!(pp.l.shape(), (m, r));
        assert_eq!(pp.u.shape(), (r, n));
        assert_eq!(pp.s.shape(), (m - r, n - r));

        let paq = a.permute(&pp.p, &pp.q);
        let rem = SpMat::from_entries((m, n),
            pp.s.iter_nz().map(|(i, j, v)| (i + r, j + r, v.clone()))
        );
        assert_eq!(paq, &pp.l * &pp.u + &rem, "p*A*q != l*u + rest");

        // U is unit upper: 1s on diagonal
        let ub = pp.u.clone().into_dense();
        for k in 0..r {
            assert!(ub[(k, k)].is_one(), "u[{k},{k}] should be 1");
            for i in k + 1..r { assert_eq!(ub[(i, k)], 0, "u[{i},{k}] below diag"); }
        }
    }

    // ---- extract_dense ----

    #[test]
    fn test_extract_dense_no_transpose() {
        // S has a zero row (row 1) and a zero col (col 1).
        // Non-zero entries: (0,0)=1, (0,2)=2, (2,0)=3, (2,2)=4.
        // row_idx=[0,2], col_idx=[0,2].
        // S0 (2×2) = [[1,2],[3,4]].
        let s = SpMat::from_row_major((3, 3), [1i32, 0, 2, 0, 0, 0, 3, 0, 4]);
        let (row_idx, col_idx, mat) = extract_dense(&s, false);
        assert_eq!(row_idx, vec![0usize, 2]);
        assert_eq!(col_idx, vec![0usize, 2]);
        assert_eq!(mat, crate::dense::Mat::from_row_major((2, 2), [1i32, 2, 3, 4]));
    }

    #[test]
    fn test_extract_dense_transpose() {
        // Same S, but with transpose=true.  S0^T (2×2) = [[1,3],[2,4]].
        let s = SpMat::from_row_major((3, 3), [1i32, 0, 2, 0, 0, 0, 3, 0, 4]);
        let (row_idx, col_idx, mat) = extract_dense(&s, true);
        assert_eq!(row_idx, vec![0usize, 2]);
        assert_eq!(col_idx, vec![0usize, 2]);
        assert_eq!(mat, crate::dense::Mat::from_row_major((2, 2), [1i32, 3, 2, 4]));
    }

    // ---- dense_pluq_in ----

    fn check_dense_pluq_in(s: &SpMat<i32>, piv_type: PivotType) {
        let pp = dense_pluq_in(s, piv_type);
        let (ms, ns) = s.shape();
        let r = pp.rank();

        assert_eq!(pp.l.shape(), (ms, r));
        assert_eq!(pp.u.shape(), (r, ns));
        assert_eq!(pp.s.shape(), (ms - r, ns - r));

        let psq = s.permute(&pp.p, &pp.q);
        let rem = SpMat::from_entries((ms, ns),
            pp.s.iter_nz().map(|(i, j, v)| (i + r, j + r, v.clone()))
        );
        assert_eq!(psq, &pp.l * &pp.u + &rem, "p*s*q != l*u + rest");
    }

    #[test]
    fn test_dense_pluq_in_cols_with_zero_row_and_col() {
        // S has a zero row (row 1) and a zero col (col 1).
        let s = SpMat::from_row_major((3, 3), [1i32, 0, 2, 0, 0, 0, 3, 0, 4]);
        check_dense_pluq_in(&s, PivotType::Cols);
    }

    #[test]
    fn test_dense_pluq_in_rows_with_zero_row_and_col() {
        let s = SpMat::from_row_major((3, 3), [1i32, 0, 2, 0, 0, 0, 3, 0, 4]);
        check_dense_pluq_in(&s, PivotType::Rows);
    }

    #[test]
    fn test_dense_pluq_in_all_zero() {
        let s = SpMat::<i32>::zero((4, 5));
        check_dense_pluq_in(&s, PivotType::Cols);
        check_dense_pluq_in(&s, PivotType::Rows);
    }

    #[test]
    fn test_dense_pluq_in_no_zero_rows_or_cols() {
        // No zero rows/cols: compact_dense gives the full matrix.
        let s = SpMat::from_row_major((3, 3), [1i32,2,3,4,5,6,7,8,9]);
        check_dense_pluq_in(&s, PivotType::Cols);
        check_dense_pluq_in(&s, PivotType::Rows);
    }

    #[test]
    fn test_pluq_rows() { check_pluq_rows(&sample()); }

    #[test]
    fn test_pluq_cols() { check_pluq_cols(&sample()); }

    #[test]
    fn test_pluq_zero() {
        let a = SpMat::<i32>::zero((4, 5));
        check_pluq_rows(&a);
        check_pluq_cols(&a);
    }

    #[test]
    fn test_pluq_rand_rows() { check_pluq_rows(&SpMat::<i32>::rand((40, 60), 0.1)); }

    #[test]
    fn test_pluq_rand_cols() { check_pluq_cols(&SpMat::<i32>::rand((40, 60), 0.1)); }

    // ---- solve_l ----

    use yui_core::num::Ratio;
    type R = Ratio<i64>;
    fn r(n: i64) -> R { R::from(n) }
    
    fn sp_mat(shape: (usize, usize), data: impl IntoIterator<Item = R>) -> SpMat<R> {
        SpMat::from_row_major(shape, data)
    }

    fn sp_vec(data: impl IntoIterator<Item = R>) -> SpVec<R> {
        SpVec::from(data.into_iter().collect::<Vec<_>>())
    }

    #[test]
    fn test_solve_l_square() {
        // l = [[2, 0], [3, 4]], y = [4, 11]
        // l[0..2,0..2]*x = [4,11] → x = [2, 5/4]
        let l = sp_mat((2, 2), [r(2), r(0), r(3), r(4)]);
        let y = vec![r(4), r(11)];
        let x = solve_l(&l, &y, true);
        assert_eq!(x, Some(vec![r(2), r(5)/r(4)]));
    }

    #[test]
    fn test_solve_l_rectangular_consistent() {
        // l = [[1,0],[2,1],[3,4]] (3×2 lower triangular with extra row), y = [1,2,3].
        // Forward sub on top 2×2: z = [1, 2 - 2*1] = [1, 0].
        // Residual at row 2: 3 - 3*1 - 4*0 = 0 → consistent.
        let l = sp_mat((3, 2), [r(1), r(0), r(2), r(1), r(3), r(4)]);
        let y = vec![r(1), r(2), r(3)];
        assert_eq!(solve_l(&l, &y, true), Some(vec![r(1), r(0)]));
    }

    #[test]
    fn test_solve_l_rectangular_inconsistent() {
        // Same l as above but y = [1,2,4]. Residual at row 2: 4 - 3 - 0 = 1 ≠ 0 → None.
        let l = sp_mat((3, 2), [r(1), r(0), r(2), r(1), r(3), r(4)]);
        let y = vec![r(1), r(2), r(4)];
        assert_eq!(solve_l(&l, &y, true), None);
    }

    #[test]
    fn test_solve_l_no_check() {
        // Same inconsistent input as above; with check=false the residual is ignored
        // and Some(z) is still returned (z is the forward-sub solution on the top).
        let l = sp_mat((3, 2), [r(1), r(0), r(2), r(1), r(3), r(4)]);
        let y = vec![r(1), r(2), r(4)];
        assert_eq!(solve_l(&l, &y, false), Some(vec![r(1), r(0)]));
    }

    #[test]
    fn test_solve_l_zero_cols_consistent() {
        // r = 0, y = 0 → returns Some(empty) (third branch with empty submat).
        let l: SpMat<R> = SpMat::zero((3, 0));
        assert_eq!(solve_l(&l, &[r(0); 3], true), Some(vec![]));
    }

    #[test]
    fn test_solve_l_zero_cols_inconsistent() {
        // r = 0 with non-zero y → residual = y ≠ 0 → None when check=true.
        let l: SpMat<R> = SpMat::zero((3, 0));
        assert_eq!(solve_l(&l, &[r(1), r(0), r(0)], true), None);
    }

    // ---- is_consistent / is_consistent_upto ----

    #[test]
    fn test_is_consistent_full() {
        // l = [[1,0],[2,1],[3,4]], x = [1, 0]:
        //   y = [1, 2, 3]: residual = [0, 0] → consistent.
        //   y = [1, 2, 4]: residual at row 2 = 1 ≠ 0 → inconsistent.
        let l = sp_mat((3, 2), [r(1), r(0), r(2), r(1), r(3), r(4)]);
        let x = vec![r(1), r(0)];
        assert!( is_consistent(&l, &[r(1), r(2), r(3)], &x));
        assert!(!is_consistent(&l, &[r(1), r(2), r(4)], &x));
    }

    #[test]
    fn test_is_consistent_upto_partial() {
        // l = [[1,0],[2,1],[3,4],[5,6]], x = [1, 0], y = [1, 2, 3, 99]:
        //   row-2 residual = 0; row-3 residual = 99 - 5 = 94.
        //   k=2: checks rows in [2..2] (none) → trivially true.
        //   k=3: checks row 2 only → true.
        //   k=4: checks rows 2 and 3 → false (row 3 fails).
        let l = sp_mat((4, 2), [r(1), r(0), r(2), r(1), r(3), r(4), r(5), r(6)]);
        let y = vec![r(1), r(2), r(3), r(99)];
        let x = vec![r(1), r(0)];
        assert!( is_consistent_upto(&l, &y, &x, 2));
        assert!( is_consistent_upto(&l, &y, &x, 3));
        assert!(!is_consistent_upto(&l, &y, &x, 4));
    }

    // ---- solve_u ----

    #[test]
    fn test_solve_u_square() {
        // u = [[1, 2], [0, 3]], y = [4, 6].
        // Back sub: x[1] = 6/3 = 2; x[0] = (4 - 2*2)/1 = 0.
        let u = sp_mat((2, 2), [r(1), r(2), r(0), r(3)]);
        let y = vec![r(4), r(6)];
        assert_eq!(solve_u(&u, &y), vec![r(0), r(2)]);
    }

    #[test]
    fn test_solve_u_rectangular() {
        // u = [[1, 2, 5, 6], [0, 3, 7, 8]], y = [4, 6].
        // Top 2×2 same as above → x[..2] = [0, 2]; trailing entries are zeros.
        let u = sp_mat((2, 4), [r(1), r(2), r(5), r(6), r(0), r(3), r(7), r(8)]);
        let y = vec![r(4), r(6)];
        assert_eq!(solve_u(&u, &y), vec![r(0), r(2), r(0), r(0)]);
    }

    #[test]
    fn test_solve_u_empty() {
        // r = 0, n = 0 → empty input/output.
        let u: SpMat<R> = SpMat::zero((0, 0));
        assert_eq!(solve_u(&u, &[]), Vec::<R>::new());
    }

    // ---- solve_pluq integration tests ----

    fn solve_check(a: &SpMat<R>, y: &SpVec<R>) -> SpVec<R> {
        let x = solve_pluq(a, y).expect("expected a solution");
        assert_eq!(&(a * &x), y, "A*x != y");
        x
    }

    #[test]
    fn test_solve_square() {
        let a = sp_mat((2, 2), [r(1), r(2), r(3), r(4)]);
        solve_check(&a, &sp_vec([r(5), r(6)]));
    }

    #[test]
    fn test_solve_overdetermined_consistent() {
        let a = sp_mat((3, 2), [r(1), r(0), r(0), r(1), r(1), r(1)]);
        solve_check(&a, &sp_vec([r(2), r(3), r(5)]));
    }

    #[test]
    fn test_solve_overdetermined_inconsistent() {
        let a = sp_mat((3, 2), [r(1), r(0), r(0), r(1), r(1), r(1)]);
        assert!(solve_pluq(&a, &sp_vec([r(1), r(1), r(0)])).is_none());
    }

    #[test]
    fn test_solve_underdetermined() {
        let a = sp_mat((2, 3), [r(1), r(0), r(2), r(0), r(1), r(3)]);
        solve_check(&a, &sp_vec([r(4), r(5)]));
    }

    #[test]
    fn test_solve_zero_rhs() {
        let a = sp_mat((2, 2), [r(1), r(2), r(3), r(4)]);
        let x = solve_check(&a, &sp_vec([r(0), r(0)]));
        assert_eq!(x, sp_vec([r(0), r(0)]));
    }

    #[test]
    fn test_solve_no_solution() {
        let a = sp_mat((2, 2), [r(1), r(2), r(2), r(4)]);
        assert!(solve_pluq(&a, &sp_vec([r(1), r(0)])).is_none());
    }

    #[test]
    fn test_solve_identity() {
        let a: SpMat<R> = SpMat::from_entries((4, 4), (0..4).map(|k| (k, k, r(1))));
        let y = sp_vec([r(1), r(2), r(3), r(4)]);
        let x = solve_check(&a, &y);
        assert_eq!(x, y);
    }

    // ---- extend_chunk_to_full ----

    fn check_extend_chunk(s: &SpMat<i32>, c: usize) {
        let (m, n) = s.shape();
        assert!(c <= m);
        let [s_top, s_rest] = s.clone().v_split(c);
        let pp_chunk = pluq(&s_top, cfg(PivotType::Rows));
        let pp = extend_chunk_to_full(pp_chunk, s_rest);
        let r = pp.rank();

        assert_eq!(pp.l.shape(), (m, r));
        assert_eq!(pp.u.shape(), (r, n));
        assert_eq!(pp.s.shape(), (m - r, n - r));

        let psq = s.permute(&pp.p, &pp.q);
        let rem = SpMat::from_entries((m, n),
            pp.s.iter_nz().map(|(i, j, v)| (i + r, j + r, v.clone()))
        );
        assert_eq!(psq, &pp.l * &pp.u + &rem, "p*s*q != l*u + rest (c = {c})");
    }

    #[test]
    fn test_extend_chunk_top() { check_extend_chunk(&sample(), 3); }

    #[test]
    fn test_extend_chunk_full() { check_extend_chunk(&sample(), 6); }

    #[test]
    fn test_extend_chunk_empty() { check_extend_chunk(&sample(), 0); }

    #[test]
    fn test_extend_chunk_rand() { check_extend_chunk(&SpMat::<i32>::rand((40, 60), 0.1), 17); }

    // ---- chunk_pluq ----

    #[test]
    fn test_chunk_pluq() {
        let s = sample();
        let (m, n) = s.shape();
        let (pp, r_chunk, c) = chunk_pluq(s.clone(), 3);
        let r = pp.rank();

        assert_eq!(c, 3);
        assert_eq!(r, r_chunk);
        assert_eq!(pp.l.shape(), (m, r));
        assert_eq!(pp.u.shape(), (r, n));
        assert_eq!(pp.s.shape(), (m - r, n - r));

        let psq = s.permute(&pp.p, &pp.q);
        let rem = SpMat::from_entries((m, n),
            pp.s.iter_nz().map(|(i, j, v)| (i + r, j + r, v.clone()))
        );
        assert_eq!(psq, &pp.l * &pp.u + &rem);
    }

    #[test]
    fn test_chunk_pluq_oversize() {
        let s = sample();
        let m = s.n_rows();
        let (_, _, c) = chunk_pluq(s, 100);
        assert_eq!(c, m);
    }

    // ---- solve_pluq_incr integration tests ----

    fn solve_incr_check(a: &SpMat<R>, y: &SpVec<R>, max_piv: usize, chunk: usize) -> SpVec<R> {
        let x = solve_pluq_incr(a, y, max_piv, chunk).expect("expected a solution");
        assert_eq!(&(a * &x), y, "A*x != y (max_piv={max_piv}, chunk={chunk})");
        x
    }

    #[test]
    fn test_solve_incr_square() {
        let a = sp_mat((2, 2), [r(1), r(2), r(3), r(4)]);
        // exercise different (max_piv, chunk) combinations
        for (mp, ch) in [(0, 1), (0, 2), (1, 1), (usize::MAX, 1), (usize::MAX, 100)] {
            solve_incr_check(&a, &sp_vec([r(5), r(6)]), mp, ch);
        }
    }

    #[test]
    fn test_solve_incr_overdetermined_consistent() {
        let a = sp_mat((3, 2), [r(1), r(0), r(0), r(1), r(1), r(1)]);
        solve_incr_check(&a, &sp_vec([r(2), r(3), r(5)]), 0, 2);
        solve_incr_check(&a, &sp_vec([r(2), r(3), r(5)]), 1, 1);
    }

    #[test]
    fn test_solve_incr_overdetermined_inconsistent() {
        let a = sp_mat((3, 2), [r(1), r(0), r(0), r(1), r(1), r(1)]);
        for (mp, ch) in [(0, 1), (0, 3), (usize::MAX, 1)] {
            assert!(solve_pluq_incr(&a, &sp_vec([r(1), r(1), r(0)]), mp, ch).is_none());
        }
    }

    #[test]
    fn test_solve_incr_underdetermined() {
        let a = sp_mat((2, 3), [r(1), r(0), r(2), r(0), r(1), r(3)]);
        solve_incr_check(&a, &sp_vec([r(4), r(5)]), 0, 1);
        solve_incr_check(&a, &sp_vec([r(4), r(5)]), 1, 1);
    }

    #[test]
    fn test_solve_incr_zero_rhs() {
        let a = sp_mat((2, 2), [r(1), r(2), r(3), r(4)]);
        let x = solve_incr_check(&a, &sp_vec([r(0), r(0)]), 0, 1);
        assert_eq!(x, sp_vec([r(0), r(0)]));
    }

    #[test]
    fn test_solve_incr_no_solution() {
        let a = sp_mat((2, 2), [r(1), r(2), r(2), r(4)]);
        for (mp, ch) in [(0, 1), (0, 2), (usize::MAX, 1)] {
            assert!(solve_pluq_incr(&a, &sp_vec([r(1), r(0)]), mp, ch).is_none());
        }
    }

    #[test]
    fn test_solve_incr_identity() {
        let a: SpMat<R> = SpMat::from_entries((4, 4), (0..4).map(|k| (k, k, r(1))));
        let y = sp_vec([r(1), r(2), r(3), r(4)]);
        let x = solve_incr_check(&a, &y, 0, 2);
        assert_eq!(x, y);
    }

    #[test]
    fn test_solve_incr_zero_matrix_zero_rhs() {
        let a = SpMat::<R>::zero((3, 4));
        // Ax = 0 with A=0 has any x as a solution; expect all-zero free vars.
        let x = solve_pluq_incr(&a, &sp_vec([r(0); 3]), 0, 1).expect("zero rhs is consistent");
        assert_eq!(x, sp_vec([r(0); 4]));
    }

    #[test]
    fn test_solve_incr_zero_matrix_nonzero_rhs() {
        let a = SpMat::<R>::zero((3, 4));
        assert!(solve_pluq_incr(&a, &sp_vec([r(1), r(0), r(0)]), 0, 1).is_none());
    }

    #[test]
    fn test_solve_incr_matches_solve_pluq() {
        // Random sparse system of moderate size; the two solvers should agree.
        let a: SpMat<R> = sp_mat((6, 9), [
            r(1), r(0), r(0), r(0), r(0), r(1), r(0), r(0), r(1),
            r(0), r(1), r(1), r(1), r(0), r(1), r(0), r(1), r(0),
            r(0), r(0), r(1), r(1), r(0), r(0), r(0), r(1), r(1),
            r(0), r(1), r(0), r(0), r(1), r(0), r(0), r(0), r(0),
            r(0), r(0), r(1), r(0), r(0), r(0), r(0), r(0), r(0),
            r(0), r(1), r(0), r(0), r(0), r(1), r(0), r(1), r(0),
        ]);
        let y = sp_vec([r(1), r(2), r(3), r(0), r(1), r(0)]);

        // Pick y that's reachable: y = a * (1, 1, ..., 1) is guaranteed consistent.
        let y_consistent = &a * &sp_vec(vec![r(1); 9]);

        for (mp, ch) in [(0, 1), (0, 3), (2, 2), (usize::MAX, 2)] {
            // Inputs may be inconsistent for this `y`; check both behave the same.
            let x_full = solve_pluq(&a, &y);
            let x_incr = solve_pluq_incr(&a, &y, mp, ch);
            assert_eq!(x_full.is_some(), x_incr.is_some(), "mp={mp}, ch={ch}");

            // Always-consistent y: both must succeed and produce a solution.
            solve_incr_check(&a, &y_consistent, mp, ch);
        }
    }

    // ---- merge_perm ----

    #[test]
    fn test_merge_perm() {
        // perm1 (size 5) = [2, 0, 3, 1, 4]; r = 2; perm2 (size 3) = [1, 2, 0].
        // For each i in 0..5, let j = perm1.at(i):
        //   i=0: j=2 ≥ r → r + perm2.at(0) = 2 + 1 = 3
        //   i=1: j=0 < r → 0
        //   i=2: j=3 ≥ r → r + perm2.at(1) = 2 + 2 = 4
        //   i=3: j=1 < r → 1
        //   i=4: j=4 ≥ r → r + perm2.at(2) = 2 + 0 = 2
        let perm1 = Perm::from_indices([2, 0, 3, 1, 4]);
        let perm2 = Perm::from_indices([1, 2, 0]);
        let p = merge_perm(&perm1, perm2);
        for (i, expected) in [3, 0, 4, 1, 2].iter().enumerate() {
            assert_eq!(p.at(i), *expected, "mismatch at i={i}");
        }
    }

    // ---- extend_perm ----

    #[test]
    fn test_extend_perm() {
        // compact_idx = [1, 3] in full space of size 5.
        // compact_perm swaps the two: at(0)=1, at(1)=0.
        // Expected:
        //   i=1 (compact_idx[0]) -> compact_perm.at(0) = 1
        //   i=3 (compact_idx[1]) -> compact_perm.at(1) = 0
        //   rest = [0,2,4] -> positions [2,3,4]
        //     i=0 -> 2,  i=2 -> 3,  i=4 -> 4
        let cp = Perm::from_indices([1, 0]);
        let idx = vec![1usize, 3];
        let p = extend_perm(5, &idx, cp);
        assert_eq!(p.at(0), 2);
        assert_eq!(p.at(1), 1);
        assert_eq!(p.at(2), 3);
        assert_eq!(p.at(3), 0);
        assert_eq!(p.at(4), 4);
    }

    #[test]
    fn test_extend_perm_identity() {
        // compact_idx = [0, 2, 5] with identity compact_perm.
        // extend_perm should equal Perm::forward_indices(7, [0,2,5]).
        let cp = Perm::id(3);
        let idx = vec![0usize, 2, 5];
        let p = extend_perm(7, &idx, cp);
        let expected = Perm::forward_indices(7, idx.iter().copied());
        for i in 0..7 {
            assert_eq!(p.at(i), expected.at(i), "mismatch at i={i}");
        }
    }

}
