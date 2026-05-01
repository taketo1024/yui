// Sparse PLUQ decomposition & linear solver.
// Implemented with the help of Claude Code.

use log::debug;
use sprs::PermOwned;
use sprs::PermView;
use yui_core::{Ring, RingOps, Field, FieldOps};

use crate::MatTrait;
use crate::dense::Mat;
use crate::dense::pluq::pluq as dense_pluq;
use super::SpMat;
use super::SpVec;
use super::pivot::{PivotFinderConfig, PivotType, find_pivots, perms_by_pivots};
use super::triang::{TriangularType, solve_triangular, solve_triangular_left, solve_triangular_vec};
use super::util::perm_for_indices;

/// Result of a sparse PLUQ decomposition.
///
/// Satisfies `p * A * q = l * u + s` where `s` is the
/// `(m - rank) × (n - rank)` Schur complement (bottom-right block).
pub struct SpPluq<R> {
    pub p: PermOwned,
    pub q: PermOwned,
    pub l: SpMat<R>,
    pub u: SpMat<R>,
    pub s: SpMat<R>,
}

impl<R> SpPluq<R> {
    /// Constructs a `PartialPluq` after asserting the shapes are mutually
    /// consistent: `l.ncols() == u.nrows() = r`, `l.nrows() == p.dim() = m`,
    /// `u.ncols() == q.dim() = n`, and `s.shape() == (m - r, n - r)`.
    pub fn new(p: PermOwned, q: PermOwned, l: SpMat<R>, u: SpMat<R>, s: SpMat<R>) -> Self {
        let r = l.ncols();
        let m = l.nrows();
        let n = u.ncols();
        assert_eq!(r, u.nrows(), "l.ncols() must match u.nrows()");
        assert_eq!(m, p.dim(), "l.nrows() must match p.dim()");
        assert_eq!(n, q.dim(), "u.ncols() must match q.dim()");
        assert_eq!(s.shape(), (m - r, n - r), "s shape must be (m - r, n - r)");
        Self { p, q, l, u, s }
    }

    pub fn rank(&self) -> usize { self.l.ncols() }
}

/// Computes a partial PLUQ decomposition of `a` under the given pivot-finder
/// configuration.
pub fn pre_pluq<R>(a: &SpMat<R>, config: PivotFinderConfig) -> SpPluq<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    debug!("compute sparse pluq: {:?}", a.shape());

    let piv_type = config.piv_type;
    let pivots = find_pivots(a, config);
    let (p, q) = perms_by_pivots(a, &pivots);
    let r = pivots.len();

    let paq = split_by_pqr(a, &p, &q, r);
    let (l, u, s) = build_lus(piv_type, paq);

    SpPluq::new(p, q, l, u, s)
}

// Applies permutations (p, q) to `a` and partitions the result into four blocks at row/col r:
//
//   paq = [[a0 | a1],   a0: r×r,     a1: r×(n-r)
//          [a2 | a3]]   a2: (m-r)×r, a3: (m-r)×(n-r)
fn split_by_pqr<R>(a: &SpMat<R>, p: &PermOwned, q: &PermOwned, r: usize) -> [SpMat<R>; 4]
where R: Ring, for<'x> &'x R: RingOps<R> {
    use std::cmp::Ordering::Less;

    let (m, n) = a.shape();
    let [mut a0, mut a1, mut a2, mut a3] = [vec![], vec![], vec![], vec![]];

    for (i, j, v) in a.iter() {
        let (pi, qj) = (p.at(i), q.at(j));
        let v = v.clone();
        match (pi.cmp(&r), qj.cmp(&r)) {
            (Less, Less) => a0.push((pi,     qj,     v)),
            (Less, _   ) => a1.push((pi,     qj - r, v)),
            (_   , Less) => a2.push((pi - r, qj,     v)),
            (_   , _   ) => a3.push((pi - r, qj - r, v)),
        }
    }
    [
        SpMat::from_entries((r,     r    ), a0),
        SpMat::from_entries((r,     n - r), a1),
        SpMat::from_entries((m - r, r    ), a2),
        SpMat::from_entries((m - r, n - r), a3),
    ]
}

// Builds (l, u, s) from the four blocks paq = [[a0|a1],[a2|a3]].
// Satisfies p*A*q = l*u + [[0,0],[0,s]].
//
// Rows: (u0,u1,r0,r1) = (a0,a1,a2,a3).  l1 = r0*u0^{-1},  l = [I_r;l1],  u = [u0|u1],  s = r1-l1*u1.
// Cols: (l0,l1,r0,r1) = (a0,a2,a1,a3).  u1 = l0^{-1}*r0,  l = [l0;l1],  u = [I_r|u1],  s = r1-l1*u1.
fn build_lus<R>(piv_type: PivotType, paq: [SpMat<R>; 4]) -> (SpMat<R>, SpMat<R>, SpMat<R>)
where R: Ring, for<'x> &'x R: RingOps<R> {
    let [a0, a1, a2, a3] = paq;
    let (m, n) = (a0.nrows() + a2.nrows(), a0.ncols() + a1.ncols());
    let r = a0.ncols();

    if r == 0 {
        return (SpMat::zero((m, 0)), SpMat::zero((0, n)), a3);
    }

    match piv_type {
        PivotType::Rows => {
            let [u0, u1] = [a0, a1];
            let [r0, r1] = [a2, a3];
            let l1 = solve_triangular_left(TriangularType::Upper, &u0, &r0);
            let l = SpMat::id(r).stack(&l1);
            let u = u0.concat(&u1);
            let s = r1 - &l1 * &u1;
            (l, u, s)
        },
        PivotType::Cols => {
            let [l0, l1] = [a0, a2];
            let [r0, r1] = [a1, a3];
            let u1 = solve_triangular(TriangularType::Lower, &l0, &r0);
            let l = l0.stack(&l1);
            let u = SpMat::id(r).concat(&u1);
            let s = r1 - &l1 * &u1;
            (l, u, s)
        }
    }
}

/// Computes a full PLUQ decomposition of `a`.
pub fn pluq<R>(a: &SpMat<R>, config: PivotFinderConfig) -> SpPluq<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let piv_type = config.piv_type;
    let pp1 = pre_pluq(a, config);
    if pp1.s.is_zero() { 
        return pp1
    }
    
    let pp2 = dense_pluq_in(&pp1.s, piv_type);

    debug!("merge pluq: {} + {}", pp1.rank(), pp2.rank());

    merge_pluq(pp1, pp2)
}

fn dense_pluq_in<R>(s: &SpMat<R>, piv_type: PivotType) -> SpPluq<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let transpose = piv_type == PivotType::Rows;
    let (ms, ns) = s.shape();
    let (row_idx, col_idx, mat) = extract_dense(s, transpose);
    let (m0, n0) = (row_idx.len(), col_idx.len());

    let raw = dense_pluq(&mat);
    let dp = if transpose { raw.transpose() } else { raw };

    let p2 = extend_perm(&dp.p, &row_idx, ms);
    let q2 = extend_perm(&dp.q, &col_idx, ns);

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

    let row_perm = perm_for_indices(s.nrows(), row_idx.iter());
    let col_perm = perm_for_indices(s.ncols(), col_idx.iter());

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
fn merge_pluq<R>(pp1: SpPluq<R>, pp2: SpPluq<R>) -> SpPluq<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (m, n) = (pp1.l.nrows(), pp1.u.ncols());
    let r1 = pp1.rank();
    let r2 = pp2.rank();

    assert_eq!(pp1.s.shape(), (m - r1, n - r1));
    assert_eq!(pp2.l.nrows(), m - r1);
    assert_eq!(pp2.u.ncols(), n - r1);

    // Fast path: when pp2 contributes no new pivots, the merged result equals
    // pp1 (its schur complement is unchanged: pp2.s has the same shape as
    // pp1.s and represents the same residual when r2 == 0).
    if r2 == 0 {
        return pp1;
    }

    let p = merge_perm(&pp1.p, &pp2.p);
    let q = merge_perm(&pp1.q, &pp2.q);

    let l = {
        let [l0, l1] = pp1.l.divide_at_row(r1);
        let l1 = l1.permute_rows(pp2.p.view());
        let zero_tr = SpMat::zero((r1, r2));
        SpMat::combine_blocks([
            &l0, &zero_tr, 
            &l1, &pp2.l
        ])
    };

    let u = {
        let [u0, u1] = pp1.u.divide_at_col(r1);
        let u1 = u1.permute_cols(pp2.q.view());
        let zero_bl = SpMat::zero((r2, r1));
        SpMat::combine_blocks([
            &u0, &u1, 
            &zero_bl, &pp2.u
        ])
    };

    let s = pp2.s;

    SpPluq::new(p, q, l, u, s)
}

/// Solves `a * x = y` over a field using sparse PLUQ.
///
/// Returns `Some(x)` if a solution exists, `None` otherwise.
pub fn solve_pluq<R>(a: &SpMat<R>, y: &[R]) -> Option<Vec<R>>
where R: Field, for<'x> &'x R: FieldOps<R> {
    assert_eq!(y.len(), a.nrows());

    debug!("solve pluq, a: {:?}", a.shape());

    let pp = pluq(a, PivotFinderConfig {
        piv_type: PivotType::Rows,
        ..Default::default()
    });

    let yp = perm_apply(pp.p.view(), y);
    let xq = solve_lu(&pp.l, &pp.u, &yp)?;
    let x = perm_apply(pp.q.inv(), &xq);

    Some(x)
}

// Solves `L * U * x = y` and returns `x` of length `n = u.ncols()` with
// entries beyond `r = l.ncols()` set to zero (free variables = 0).
//
// Requires the top r × r block of L to be unit lower triangular and the top
// r × r block of U to be invertible upper triangular.
//
// Returns `None` when `solve_l(l, y, true)` detects an inconsistent residual.
// When `l` is square (`l.nrows() == r`) the residual is empty and the call
// always succeeds.
fn solve_lu<R>(l: &SpMat<R>, u: &SpMat<R>, y: &[R]) -> Option<Vec<R>>
where R: Field, for<'x> &'x R: FieldOps<R> {
    assert_eq!(l.ncols(), u.nrows());
    assert_eq!(y.len(), l.nrows());

    let z = solve_l(l, y, true)?;
    let x = solve_u(u, &z);

    Some(x)
}

// Solves `l[0..r, 0..r] * z = y[0..r]` by forward substitution, where
// `r = l.ncols()`. The top r × r block of L must be lower triangular with
// non-zero diagonal.
//
// If `check_consistency` is true and `r < y.len()`, also verifies the residual
// `y[r..] - l[r.., :] * z` is zero, returning `None` when it isn't. When `r ==
// y.len()` the residual is trivially empty so the check is skipped.
fn solve_l<R>(l: &SpMat<R>, y: &[R], check_consistency: bool) -> Option<Vec<R>>
where R: Field, for<'x> &'x R: FieldOps<R> {
    assert_eq!(l.nrows(), y.len());
    let r = l.ncols();

    let x = if r == y.len() { 
        let y = SpVec::from(y.to_vec());
        solve_triangular_vec(TriangularType::Lower, l, &y).to_dense()
    } else { 
        let l0 = l.submat(0..r, 0..r);
        let y0 = SpVec::from(y[..r].to_vec());
        let x = solve_triangular_vec(TriangularType::Lower, &l0, &y0).to_dense();

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
    assert_eq!(l.nrows(), y.len());
    assert_eq!(l.ncols(), x.len());
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
// `r = u.nrows()`, and returns `x` of length `n = u.ncols()` with entries
// beyond `r` set to zero. The top r × r block of U must be upper triangular
// with non-zero diagonal.
fn solve_u<R>(u: &SpMat<R>, y: &[R]) -> Vec<R>
where R: Field, for<'x> &'x R: FieldOps<R> {
    let (r, n) = u.shape();
    assert_eq!(y.len(), r);
    assert!(n >= r);

    let mut x = if n == r {
        solve_triangular_vec(TriangularType::Upper, u, &SpVec::from(y.to_vec())).to_dense()
    } else {
        let u0 = u.submat(0..r, 0..r);
        solve_triangular_vec(TriangularType::Upper, &u0, &SpVec::from(y.to_vec())).to_dense()
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
pub fn solve_pluq_incr<R>(a: &SpMat<R>, y: &[R], max_piv: usize, chunk: usize) -> Option<Vec<R>>
where R: Field, for<'x> &'x R: FieldOps<R> {
    assert_eq!(y.len(), a.nrows());

    debug!("solve pluq (incremental), a: {:?}", a.shape());

    let mut pp = pre_pluq(a, PivotFinderConfig {
        piv_type: PivotType::Rows,
        max_pivots: max_piv,
        ..Default::default()
    });
    let mut yp = perm_apply(pp.p.view(), y);

    let mut step = 1;
    let total_step = (a.nrows() - pp.rank()) / chunk;

    while pp.s.nrows() > 0 {
        debug!("solve pluq ({}/{})", step, total_step);
        debug!("  current rank: {}", pp.rank());

        let r_old = pp.rank();
        let (pp_next, r_next, c) = chunk_pluq(&pp.s, chunk);
        let p_next = pp_next.p.clone();
        
        debug!("merge pluq: {} + {}", pp.rank(), pp_next.rank());

        pp = merge_pluq(pp, pp_next);

        // Apply the chunk's row perm to the tail of yp so it stays in sync with pp.l.
        let yp_tail = perm_apply(p_next.view(), &yp[r_old..]);
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

    let xq = solve_lu(&pp.l, &pp.u, &yp)?;
    let x = perm_apply(pp.q.inv(), &xq);

    Some(x)
}

// Takes the top `min(chunk_size, s.nrows())` rows of `s`, runs `pluq` on them,
// and lifts the result to act on all of `s` via `extend_chunk_to_full`.
// Returns `(pp_chunk_full, r_chunk, c)`.
fn chunk_pluq<R>(s: &SpMat<R>, chunk_size: usize) -> (SpPluq<R>, usize, usize)
where R: Ring, for<'x> &'x R: RingOps<R> {
    let c = chunk_size.min(s.nrows());
    let [s_chunk, s_rest] = s.divide_at_row(c);
    let pp_chunk = pluq(&s_chunk, PivotFinderConfig {
        piv_type: PivotType::Rows,
        ..Default::default()
    });
    let r_chunk = pp_chunk.rank();
    let pp_full = extend_chunk_to_full(pp_chunk, &s_rest);
    (pp_full, r_chunk, c)
}

// Lifts a PLUQ of the top `c` rows of some matrix `s` (= `pp_chunk`) to a
// partial PLUQ acting on all rows of `s`, by absorbing the untouched bottom
// rows `s_rest` (shape (m_s - c, n_s)) into L (below the new pivots) and into
// the new schur complement.
//
// Sparse analog of `dense_pluq_in`, applied to a single chunk.
fn extend_chunk_to_full<R>(pp_chunk: SpPluq<R>, s_rest: &SpMat<R>) -> SpPluq<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (c, n_s) = (pp_chunk.l.nrows(), pp_chunk.u.ncols());
    let r_chunk = pp_chunk.rank();
    let m_rest = s_rest.nrows();
    let m_s = c + m_rest;

    assert_eq!(s_rest.ncols(), n_s);
    assert_eq!(pp_chunk.s.shape(), (c - r_chunk, n_s - r_chunk));

    let s_rest_q = s_rest.permute_cols(pp_chunk.q.view());
    let [s_rest_left, s_rest_right] = s_rest_q.divide_at_col(r_chunk);

    let [u_top, u_right] = pp_chunk.u.divide_at_col(r_chunk);

    let l_ext = solve_triangular_left(TriangularType::Upper, &u_top, &s_rest_left);
    let s_ext = s_rest_right - &l_ext * &u_right;

    let chunk_idx: Vec<usize> = (0..c).collect();
    let p = extend_perm(&pp_chunk.p, &chunk_idx, m_s);
    let q = pp_chunk.q;
    let l = pp_chunk.l.stack(&l_ext);
    let u = u_top.concat(&u_right);
    let s = pp_chunk.s.stack(&s_ext);

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
    let m = pp.l.nrows();
    assert!(r + k <= m);
    assert_eq!(yp.len(), m);
    assert_eq!(pp.p.dim(), m);

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
    pp.s = pp.s.submat_rows(k..pp.s.nrows());

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
    pp.p = PermOwned::new(new_p_at);
}

// Composes perm1 with perm2: the first `r` positions stay, the rest are
// shifted by `r` and remapped by perm2 (where `r = perm1.dim() - perm2.dim()`).
fn merge_perm(perm1: &PermOwned, perm2: &PermOwned) -> PermOwned {
    assert!(perm1.dim() >= perm2.dim());
    let n = perm1.dim();
    let r = n - perm2.dim();
    PermOwned::new((0..n).map(|i| {
        let j = perm1.at(i);
        if j < r { j } else { r + perm2.at(j - r) }
    }).collect())
}

// Lifts a compact permutation (acting on compact_idx elements of [0..full_n]) to the
// full index space.  compact_idx[k] maps to compact_perm.at(k) (within [0..mr]);
// all other indices map to consecutive positions starting at mr (in sorted order).
fn extend_perm(compact_perm: &PermOwned, compact_idx: &[usize], full_n: usize) -> PermOwned {
    let front = perm_for_indices(full_n, compact_idx.iter());
    let mr = compact_idx.len();
    PermOwned::new((0..full_n).map(|i| {
        let k = front.at(i);
        if k < mr { compact_perm.at(k) } else { k }
    }).collect())
}

// Applies permutation p to y: yp[p(i)] = y[i].
fn perm_apply<R>(p: PermView, y: &[R]) -> Vec<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let pinv = p.inv();
    (0..y.len()).map(|i| y[pinv.at(i)].clone()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::One;

    fn cfg(piv_type: PivotType) -> PivotFinderConfig {
        PivotFinderConfig { piv_type, ..Default::default() }
    }

    fn sample() -> SpMat<i32> {
        SpMat::from_dense_data((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0,
        ])
    }

    // ---- split_by_pqr ----

    #[test]
    fn test_split() {
        use sprs::PermOwned;
        // a = [[1,2],[3,4]], r=1, identity perms → paq = a, partition at row/col 1:
        // a0=[[1]], a1=[[2]], a2=[[3]], a3=[[4]]
        let a = sp((2, 2), [r(1), r(2), r(3), r(4)]);
        let p = PermOwned::new(vec![0, 1]);
        let q = PermOwned::new(vec![0, 1]);
        let [a0, a1, a2, a3] = split_by_pqr(&a, &p, &q, 1);
        assert_eq!(a0, sp((1, 1), [r(1)]));
        assert_eq!(a1, sp((1, 1), [r(2)]));
        assert_eq!(a2, sp((1, 1), [r(3)]));
        assert_eq!(a3, sp((1, 1), [r(4)]));
    }

    // ---- build_lus ----

    #[test]
    fn test_build_cols() {
        // paq = [[1,2],[3,4]], r=1: a0=[[1]], a1=[[2]], a2=[[3]], a3=[[4]]
        // l0=a0=[[1]], l1=a2=[[3]], r0=a1=[[2]], r1=a3=[[4]]
        // u1 = [[1]]^{-1}*[[2]] = [[2]], l = [[1],[3]], u = [[1,2]], s = [[4]]-[[3]]*[[2]] = [[-2]]
        let paq = [sp((1,1),[r(1)]), sp((1,1),[r(2)]), sp((1,1),[r(3)]), sp((1,1),[r(4)])];
        let (l, u, s) = build_lus(PivotType::Cols, paq);
        assert_eq!(l, sp((2, 1), [r(1), r(3)]));
        assert_eq!(u, sp((1, 2), [r(1), r(2)]));
        assert_eq!(s, sp((1, 1), [r(-2)]));
    }

    #[test]
    fn test_build_rows() {
        // paq = [[1,2],[3,4]], r=1: a0=[[1]], a1=[[2]], a2=[[3]], a3=[[4]]
        // u0=a0=[[1]], u1=a1=[[2]], r0=a2=[[3]], r1=a3=[[4]]
        // l1 = [[3]]*[[1]]^{-1} = [[3]], l = [[1],[3]], u = [[1,2]], s = [[4]]-[[3]]*[[2]] = [[-2]]
        let paq = [sp((1,1),[r(1)]), sp((1,1),[r(2)]), sp((1,1),[r(3)]), sp((1,1),[r(4)])];
        let (l, u, s) = build_lus(PivotType::Rows, paq);
        assert_eq!(l, sp((2, 1), [r(1), r(3)]));
        assert_eq!(u, sp((1, 2), [r(1), r(2)]));
        assert_eq!(s, sp((1, 1), [r(-2)]));
    }

    // ---- pre_pluq ----

    fn check_pre_pluq_rows(a: &SpMat<i32>) {
        let pp = pre_pluq(a, cfg(PivotType::Rows));
        let (m, n) = a.shape();
        let r = pp.rank();

        assert_eq!(pp.l.shape(), (m, r));
        assert_eq!(pp.u.shape(), (r, n));
        assert_eq!(pp.s.shape(), (m - r, n - r));

        let paq = a.permute(pp.p.view(), pp.q.view());
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

        let paq = a.permute(pp.p.view(), pp.q.view());
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
        assert_eq!(pp.s, a.permute(pp.p.view(), pp.q.view()));
    }

    #[test]
    fn test_pre_pluq_square_full_rank() {
        let a = SpMat::from_dense_data((3, 3), [1, 0, 0, 0, 1, 0, 0, 0, 1]);
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

        let paq = a.permute(pp.p.view(), pp.q.view());
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

        let paq = a.permute(pp.p.view(), pp.q.view());
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
        let s = SpMat::from_dense_data((3, 3), [1i32, 0, 2, 0, 0, 0, 3, 0, 4]);
        let (row_idx, col_idx, mat) = extract_dense(&s, false);
        assert_eq!(row_idx, vec![0usize, 2]);
        assert_eq!(col_idx, vec![0usize, 2]);
        assert_eq!(mat, crate::dense::Mat::from_data((2, 2), [1i32, 2, 3, 4]));
    }

    #[test]
    fn test_extract_dense_transpose() {
        // Same S, but with transpose=true.  S0^T (2×2) = [[1,3],[2,4]].
        let s = SpMat::from_dense_data((3, 3), [1i32, 0, 2, 0, 0, 0, 3, 0, 4]);
        let (row_idx, col_idx, mat) = extract_dense(&s, true);
        assert_eq!(row_idx, vec![0usize, 2]);
        assert_eq!(col_idx, vec![0usize, 2]);
        assert_eq!(mat, crate::dense::Mat::from_data((2, 2), [1i32, 3, 2, 4]));
    }

    // ---- dense_pluq_in ----

    fn check_dense_pluq_in(s: &SpMat<i32>, piv_type: PivotType) {
        let pp = dense_pluq_in(s, piv_type);
        let (ms, ns) = s.shape();
        let r = pp.rank();

        assert_eq!(pp.l.shape(), (ms, r));
        assert_eq!(pp.u.shape(), (r, ns));
        assert_eq!(pp.s.shape(), (ms - r, ns - r));

        let psq = s.permute(pp.p.view(), pp.q.view());
        let rem = SpMat::from_entries((ms, ns),
            pp.s.iter_nz().map(|(i, j, v)| (i + r, j + r, v.clone()))
        );
        assert_eq!(psq, &pp.l * &pp.u + &rem, "p*s*q != l*u + rest");
    }

    #[test]
    fn test_dense_pluq_in_cols_with_zero_row_and_col() {
        // S has a zero row (row 1) and a zero col (col 1).
        let s = SpMat::from_dense_data((3, 3), [1i32, 0, 2, 0, 0, 0, 3, 0, 4]);
        check_dense_pluq_in(&s, PivotType::Cols);
    }

    #[test]
    fn test_dense_pluq_in_rows_with_zero_row_and_col() {
        let s = SpMat::from_dense_data((3, 3), [1i32, 0, 2, 0, 0, 0, 3, 0, 4]);
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
        let s = SpMat::from_dense_data((3, 3), [1i32,2,3,4,5,6,7,8,9]);
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
    
    fn sp(shape: (usize, usize), data: impl IntoIterator<Item = R>) -> SpMat<R> {
        SpMat::from_dense_data(shape, data)
    }

    #[test]
    fn test_solve_l_square() {
        // l = [[2, 0], [3, 4]], y = [4, 11]
        // l[0..2,0..2]*x = [4,11] → x = [2, 5/4]
        let l = sp((2, 2), [r(2), r(0), r(3), r(4)]);
        let y = vec![r(4), r(11)];
        let x = solve_l(&l, &y, true);
        assert_eq!(x, Some(vec![r(2), r(5)/r(4)]));
    }

    #[test]
    fn test_solve_l_rectangular_consistent() {
        // l = [[1,0],[2,1],[3,4]] (3×2 lower triangular with extra row), y = [1,2,3].
        // Forward sub on top 2×2: z = [1, 2 - 2*1] = [1, 0].
        // Residual at row 2: 3 - 3*1 - 4*0 = 0 → consistent.
        let l = sp((3, 2), [r(1), r(0), r(2), r(1), r(3), r(4)]);
        let y = vec![r(1), r(2), r(3)];
        assert_eq!(solve_l(&l, &y, true), Some(vec![r(1), r(0)]));
    }

    #[test]
    fn test_solve_l_rectangular_inconsistent() {
        // Same l as above but y = [1,2,4]. Residual at row 2: 4 - 3 - 0 = 1 ≠ 0 → None.
        let l = sp((3, 2), [r(1), r(0), r(2), r(1), r(3), r(4)]);
        let y = vec![r(1), r(2), r(4)];
        assert_eq!(solve_l(&l, &y, true), None);
    }

    #[test]
    fn test_solve_l_no_check() {
        // Same inconsistent input as above; with check=false the residual is ignored
        // and Some(z) is still returned (z is the forward-sub solution on the top).
        let l = sp((3, 2), [r(1), r(0), r(2), r(1), r(3), r(4)]);
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
        let l = sp((3, 2), [r(1), r(0), r(2), r(1), r(3), r(4)]);
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
        let l = sp((4, 2), [r(1), r(0), r(2), r(1), r(3), r(4), r(5), r(6)]);
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
        let u = sp((2, 2), [r(1), r(2), r(0), r(3)]);
        let y = vec![r(4), r(6)];
        assert_eq!(solve_u(&u, &y), vec![r(0), r(2)]);
    }

    #[test]
    fn test_solve_u_rectangular() {
        // u = [[1, 2, 5, 6], [0, 3, 7, 8]], y = [4, 6].
        // Top 2×2 same as above → x[..2] = [0, 2]; trailing entries are zeros.
        let u = sp((2, 4), [r(1), r(2), r(5), r(6), r(0), r(3), r(7), r(8)]);
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

    fn solve_check(a: &SpMat<R>, y: &[R]) -> Vec<R> {
        let x = solve_pluq(a, y).expect("expected a solution");
        let (m, _n) = a.shape();
        let mut ax = vec![r(0); m];
        for (i, j, v) in a.iter_nz() { ax[i] = ax[i].clone() + v * &x[j]; }
        for i in 0..m {
            assert_eq!(ax[i], y[i], "row {i}: (A*x)[{i}] != y[{i}]");
        }
        x
    }

    #[test]
    fn test_solve_square() {
        let a = sp((2, 2), [r(1), r(2), r(3), r(4)]);
        solve_check(&a, &[r(5), r(6)]);
    }

    #[test]
    fn test_solve_overdetermined_consistent() {
        let a = sp((3, 2), [r(1), r(0), r(0), r(1), r(1), r(1)]);
        solve_check(&a, &[r(2), r(3), r(5)]);
    }

    #[test]
    fn test_solve_overdetermined_inconsistent() {
        let a = sp((3, 2), [r(1), r(0), r(0), r(1), r(1), r(1)]);
        assert!(solve_pluq(&a, &[r(1), r(1), r(0)]).is_none());
    }

    #[test]
    fn test_solve_underdetermined() {
        let a = sp((2, 3), [r(1), r(0), r(2), r(0), r(1), r(3)]);
        solve_check(&a, &[r(4), r(5)]);
    }

    #[test]
    fn test_solve_zero_rhs() {
        let a = sp((2, 2), [r(1), r(2), r(3), r(4)]);
        let x = solve_check(&a, &[r(0), r(0)]);
        assert_eq!(x, vec![r(0), r(0)]);
    }

    #[test]
    fn test_solve_no_solution() {
        let a = sp((2, 2), [r(1), r(2), r(2), r(4)]);
        assert!(solve_pluq(&a, &[r(1), r(0)]).is_none());
    }

    #[test]
    fn test_solve_identity() {
        let a: SpMat<R> = SpMat::from_entries((4, 4), (0..4).map(|k| (k, k, r(1))));
        let y = [r(1), r(2), r(3), r(4)];
        let x = solve_check(&a, &y);
        assert_eq!(x, y);
    }

    // ---- extend_chunk_to_full ----

    fn check_extend_chunk(s: &SpMat<i32>, c: usize) {
        let (m, n) = s.shape();
        assert!(c <= m);
        let [s_top, s_rest] = s.divide_at_row(c);
        let pp_chunk = pluq(&s_top, cfg(PivotType::Rows));
        let pp = extend_chunk_to_full(pp_chunk, &s_rest);
        let r = pp.rank();

        assert_eq!(pp.l.shape(), (m, r));
        assert_eq!(pp.u.shape(), (r, n));
        assert_eq!(pp.s.shape(), (m - r, n - r));

        let psq = s.permute(pp.p.view(), pp.q.view());
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
        let (pp, r_chunk, c) = chunk_pluq(&s, 3);
        let r = pp.rank();

        assert_eq!(c, 3);
        assert_eq!(r, r_chunk);
        assert_eq!(pp.l.shape(), (m, r));
        assert_eq!(pp.u.shape(), (r, n));
        assert_eq!(pp.s.shape(), (m - r, n - r));

        let psq = s.permute(pp.p.view(), pp.q.view());
        let rem = SpMat::from_entries((m, n),
            pp.s.iter_nz().map(|(i, j, v)| (i + r, j + r, v.clone()))
        );
        assert_eq!(psq, &pp.l * &pp.u + &rem);
    }

    #[test]
    fn test_chunk_pluq_oversize() {
        let s = sample();
        let m = s.nrows();
        let (_, _, c) = chunk_pluq(&s, 100);
        assert_eq!(c, m);
    }

    // ---- solve_pluq_incr integration tests ----

    fn solve_incr_check(a: &SpMat<R>, y: &[R], max_piv: usize, chunk: usize) -> Vec<R> {
        let x = solve_pluq_incr(a, y, max_piv, chunk).expect("expected a solution");
        let (m, _) = a.shape();
        let mut ax = vec![r(0); m];
        for (i, j, v) in a.iter_nz() { ax[i] = ax[i].clone() + v * &x[j]; }
        for i in 0..m {
            assert_eq!(ax[i], y[i], "row {i}: (A*x)[{i}] != y[{i}] (max_piv={max_piv}, chunk={chunk})");
        }
        x
    }

    #[test]
    fn test_solve_incr_square() {
        let a = sp((2, 2), [r(1), r(2), r(3), r(4)]);
        // exercise different (max_piv, chunk) combinations
        for (mp, ch) in [(0, 1), (0, 2), (1, 1), (usize::MAX, 1), (usize::MAX, 100)] {
            solve_incr_check(&a, &[r(5), r(6)], mp, ch);
        }
    }

    #[test]
    fn test_solve_incr_overdetermined_consistent() {
        let a = sp((3, 2), [r(1), r(0), r(0), r(1), r(1), r(1)]);
        solve_incr_check(&a, &[r(2), r(3), r(5)], 0, 2);
        solve_incr_check(&a, &[r(2), r(3), r(5)], 1, 1);
    }

    #[test]
    fn test_solve_incr_overdetermined_inconsistent() {
        let a = sp((3, 2), [r(1), r(0), r(0), r(1), r(1), r(1)]);
        for (mp, ch) in [(0, 1), (0, 3), (usize::MAX, 1)] {
            assert!(solve_pluq_incr(&a, &[r(1), r(1), r(0)], mp, ch).is_none());
        }
    }

    #[test]
    fn test_solve_incr_underdetermined() {
        let a = sp((2, 3), [r(1), r(0), r(2), r(0), r(1), r(3)]);
        solve_incr_check(&a, &[r(4), r(5)], 0, 1);
        solve_incr_check(&a, &[r(4), r(5)], 1, 1);
    }

    #[test]
    fn test_solve_incr_zero_rhs() {
        let a = sp((2, 2), [r(1), r(2), r(3), r(4)]);
        let x = solve_incr_check(&a, &[r(0), r(0)], 0, 1);
        assert_eq!(x, vec![r(0), r(0)]);
    }

    #[test]
    fn test_solve_incr_no_solution() {
        let a = sp((2, 2), [r(1), r(2), r(2), r(4)]);
        for (mp, ch) in [(0, 1), (0, 2), (usize::MAX, 1)] {
            assert!(solve_pluq_incr(&a, &[r(1), r(0)], mp, ch).is_none());
        }
    }

    #[test]
    fn test_solve_incr_identity() {
        let a: SpMat<R> = SpMat::from_entries((4, 4), (0..4).map(|k| (k, k, r(1))));
        let y = [r(1), r(2), r(3), r(4)];
        let x = solve_incr_check(&a, &y, 0, 2);
        assert_eq!(x, y);
    }

    #[test]
    fn test_solve_incr_zero_matrix_zero_rhs() {
        let a = SpMat::<R>::zero((3, 4));
        // Ax = 0 with A=0 has any x as a solution; expect all-zero free vars.
        let x = solve_pluq_incr(&a, &vec![r(0); 3], 0, 1).expect("zero rhs is consistent");
        assert_eq!(x, vec![r(0); 4]);
    }

    #[test]
    fn test_solve_incr_zero_matrix_nonzero_rhs() {
        let a = SpMat::<R>::zero((3, 4));
        assert!(solve_pluq_incr(&a, &[r(1), r(0), r(0)], 0, 1).is_none());
    }

    #[test]
    fn test_solve_incr_matches_solve_pluq() {
        // Random sparse system of moderate size; the two solvers should agree.
        let a: SpMat<R> = sp((6, 9), [
            r(1), r(0), r(0), r(0), r(0), r(1), r(0), r(0), r(1),
            r(0), r(1), r(1), r(1), r(0), r(1), r(0), r(1), r(0),
            r(0), r(0), r(1), r(1), r(0), r(0), r(0), r(1), r(1),
            r(0), r(1), r(0), r(0), r(1), r(0), r(0), r(0), r(0),
            r(0), r(0), r(1), r(0), r(0), r(0), r(0), r(0), r(0),
            r(0), r(1), r(0), r(0), r(0), r(1), r(0), r(1), r(0),
        ]);
        let y = vec![r(1), r(2), r(3), r(0), r(1), r(0)];

        // Pick y that's reachable: y = a * (1, 1, ..., 1) is guaranteed consistent.
        let mut y_consistent = vec![r(0); 6];
        for (i, _, v) in a.iter_nz() { y_consistent[i] = y_consistent[i].clone() + v.clone(); }

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
        use sprs::PermOwned;
        // perm1 (size 5) = [2, 0, 3, 1, 4]; r = 2; perm2 (size 3) = [1, 2, 0].
        // For each i in 0..5, let j = perm1.at(i):
        //   i=0: j=2 ≥ r → r + perm2.at(0) = 2 + 1 = 3
        //   i=1: j=0 < r → 0
        //   i=2: j=3 ≥ r → r + perm2.at(1) = 2 + 2 = 4
        //   i=3: j=1 < r → 1
        //   i=4: j=4 ≥ r → r + perm2.at(2) = 2 + 0 = 2
        let perm1 = PermOwned::new(vec![2, 0, 3, 1, 4]);
        let perm2 = PermOwned::new(vec![1, 2, 0]);
        let p = merge_perm(&perm1, &perm2);
        for (i, expected) in [3, 0, 4, 1, 2].iter().enumerate() {
            assert_eq!(p.at(i), *expected, "mismatch at i={i}");
        }
    }

    // ---- extend_perm ----

    #[test]
    fn test_extend_perm() {
        use sprs::PermOwned;
        // compact_idx = [1, 3] in full space of size 5.
        // compact_perm swaps the two: at(0)=1, at(1)=0.
        // Expected:
        //   i=1 (compact_idx[0]) -> compact_perm.at(0) = 1
        //   i=3 (compact_idx[1]) -> compact_perm.at(1) = 0
        //   rest = [0,2,4] -> positions [2,3,4]
        //     i=0 -> 2,  i=2 -> 3,  i=4 -> 4
        let cp = PermOwned::new(vec![1, 0]);
        let idx = vec![1usize, 3];
        let p = extend_perm(&cp, &idx, 5);
        assert_eq!(p.at(0), 2);
        assert_eq!(p.at(1), 1);
        assert_eq!(p.at(2), 3);
        assert_eq!(p.at(3), 0);
        assert_eq!(p.at(4), 4);
    }

    #[test]
    fn test_extend_perm_identity() {
        use sprs::PermOwned;
        // compact_idx = [0, 2, 5] with identity compact_perm.
        // extend_perm should equal perm_for_indices(7, [0,2,5]).
        let cp = PermOwned::identity(3);
        let idx = vec![0usize, 2, 5];
        let p = extend_perm(&cp, &idx, 7);
        let expected = perm_for_indices(7, idx.iter());
        for i in 0..7 {
            assert_eq!(p.at(i), expected.at(i), "mismatch at i={i}");
        }
    }

    // ---- perm_apply ----

    #[test]
    fn test_perm_apply() {
        use sprs::PermOwned;
        let p = PermOwned::new(vec![1, 2, 0]); // 0→1, 1→2, 2→0
        let y = vec![r(10), r(20), r(30)];
        let yp = perm_apply(p.view(), &y);
        // yp[p(0)=1]=10, yp[p(1)=2]=20, yp[p(2)=0]=30
        assert_eq!(yp, vec![r(30), r(10), r(20)]);
    }
}
