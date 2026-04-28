// Implemented with the help of Claude Code.

use sprs::PermOwned;
use yui_core::{Ring, RingOps, Field, FieldOps};
use crate::MatTrait;
use crate::dense::Mat;

/// Result of a PLUQ decomposition satisfying `p_mat * A * q_mat = L * U + s`:
///   - `p`: row permutation (pivot rows first)
///   - `q`: column permutation (pivot columns first)
///   - `l`: `m × rank`, unit lower triangular — elimination multipliers
///   - `u`: `rank × n`, upper echelon — the reduced pivot rows
///   - `s`: `(m - rank) × (n - rank)`, Schur complement — zero when `R` is a field
pub struct Pluq<R> {
    pub p: PermOwned,
    pub q: PermOwned,
    pub l: Mat<R>,
    pub u: Mat<R>,
    pub s: Mat<R>,
}

impl<R> Pluq<R> {
    pub fn rank(&self) -> usize { self.l.ncols() }
}

/// Computes a PLUQ decomposition of `a` over a ring by Gaussian elimination.
///
/// Only unit elements are used as pivots, so every elimination step is exact.
/// Satisfies `p_mat * a * q_mat = l * u + s`, where `s` is the Schur complement
/// and is zero when `R` is a field (every non-zero element is a unit).
pub fn pluq<R>(a: &Mat<R>) -> Pluq<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (m, n) = a.shape();
    let mut work = a.clone();
    let mut row_of: Vec<usize> = (0..m).collect();

    let (pivot_cols, l) = reduce(&mut work, &mut row_of);
    let rank = pivot_cols.len();
    let p = row_perm(&row_of, m);
    let cols = col_order(&pivot_cols, n);
    let q = col_perm(&cols, n);
    let u = build_u(&work, &cols, rank);
    let s = build_s(&work, &cols, rank);

    Pluq { p, q, l, u, s }
}

/// Solves `A * x = y` over a field using PLUQ decomposition.
///
/// Returns `Some(x)` if a solution exists, `None` otherwise.
pub fn solve_pluq<R>(a: &Mat<R>, y: &[R]) -> Option<Vec<R>>
where R: Field, for<'x> &'x R: FieldOps<R> {
    assert_eq!(y.len(), a.nrows());
    let Pluq { p, q, l, u, .. } = pluq(a);
    let yp = apply_perm(&p, y);
    let z  = forward_sub(&l, &yp);
    if !check_consistent(&l, &yp, &z) { return None; }
    let xp = back_sub(&u, &z);
    Some((0..xp.len()).map(|j| xp[q.at(j)].clone()).collect())
}

// Applies permutation p to y: result[k] = y[p^{-1}(k)], i.e., result[p(i)] = y[i].
fn apply_perm<R: Clone>(p: &PermOwned, y: &[R]) -> Vec<R> {
    let pinv = p.inv();
    (0..y.len()).map(|k| y[pinv.at(k)].clone()).collect()
}

// Solves L * z = yp[0..rank] by forward substitution (L is unit lower triangular).
fn forward_sub<R>(l: &Mat<R>, yp: &[R]) -> Vec<R>
where R: Field, for<'x> &'x R: FieldOps<R> {
    (0..l.ncols()).fold(vec![], |mut z, k| {
        let val = (0..k).fold(yp[k].clone(), |v, j| v - &l[(k, j)] * &z[j]);
        z.push(val);
        z
    })
}

// Returns true if L[i,:] * z == yp[i] for every non-pivot row i >= rank.
fn check_consistent<R>(l: &Mat<R>, yp: &[R], z: &[R]) -> bool
where R: Ring + PartialEq, for<'x> &'x R: RingOps<R> {
    let rank = z.len();
    (rank..yp.len()).all(|i| {
        (0..rank).fold(R::zero(), |acc, j| acc + &l[(i, j)] * &z[j]) == yp[i]
    })
}

// Solves U * xp = z by back substitution; free variables xp[rank..n] stay zero.
fn back_sub<R>(u: &Mat<R>, z: &[R]) -> Vec<R>
where R: Field, for<'x> &'x R: FieldOps<R> {
    let (rank, n) = (z.len(), u.ncols());
    (0..rank).rev().fold(vec![R::zero(); n], |mut xp, k| {
        xp[k] = (k + 1..rank).fold(z[k].clone(), |v, j| v - &u[(k, j)] * &xp[j])
            * u[(k, k)].inv().unwrap();
        xp
    })
}

// Gaussian elimination in place, eliminating only below each pivot.
// Only unit elements are accepted as pivots so that `inv()` is always valid.
// Builds L column by column: each pivot step contributes one column
// (unit on the diagonal, multipliers below, zeros above).
// Any row swaps are mirrored in the already-built L columns.
//
// Returns (pivot_col_indices, L).
fn reduce<R>(work: &mut Mat<R>, row_of: &mut Vec<usize>) -> (Vec<usize>, Mat<R>)
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (m, n) = work.shape();
    let mut pivot_cols = Vec::new();
    let mut l_cols: Vec<Vec<R>> = Vec::new();
    let mut r = 0;

    for j in 0..n {
        if r >= m { break; }
        let Some(pivot_pos) = (r..m).find(|&i| work[(i, j)].is_unit()) else { continue; };
        if pivot_pos != r {
            work.swap_rows(r, pivot_pos);
            row_of.swap(r, pivot_pos);
            l_cols.iter_mut().for_each(|col| col.swap(r, pivot_pos));
        }
        let l_col = build_l_col(work, r, j);
        eliminate_below(work, &l_col, r);
        l_cols.push(l_col);
        pivot_cols.push(j);
        r += 1;
    }

    let rank = pivot_cols.len();
    let l = Mat::from_generator((m, rank), |i, k| l_cols[k][i].clone());
    (pivot_cols, l)
}

// Builds column r of L: 1 on the diagonal, multipliers (entry * pivot_inv) below, 0 above.
fn build_l_col<R>(work: &Mat<R>, r: usize, j: usize) -> Vec<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    use std::cmp::Ordering::*;
    let pivot_inv = work[(r, j)].inv().unwrap();
    (0..work.nrows()).map(|i| match i.cmp(&r) {
        Less    => R::zero(),
        Equal   => R::one(),
        Greater => work[(i, j)].clone() * pivot_inv.clone(),
    }).collect()
}

// Subtracts `l_col[i] * row_r` from each row `i > r`, zeroing out the pivot column below.
fn eliminate_below<R>(work: &mut Mat<R>, l_col: &[R], r: usize)
where R: Ring, for<'x> &'x R: RingOps<R> {
    let m = work.nrows();
    (r + 1..m)
        .filter(|&i| !l_col[i].is_zero())
        .for_each(|i| work.add_row_to(r, i, &-l_col[i].clone()));
}

// Builds the row permutation: p.at(orig) = current position of that row.
fn row_perm(row_of: &[usize], m: usize) -> PermOwned {
    PermOwned::new(row_of.iter().enumerate().fold(vec![0usize; m], |mut v, (pos, &orig)| {
        v[orig] = pos; v
    }))
}

// Returns the full column reordering: pivot columns first, non-pivot columns last.
fn col_order(pivot_cols: &[usize], n: usize) -> Vec<usize> {
    use std::collections::HashSet;
    let pivot_set: HashSet<usize> = pivot_cols.iter().cloned().collect();
    pivot_cols.iter().cloned().chain((0..n).filter(|j| !pivot_set.contains(j))).collect()
}

// Builds the column permutation from the full ordered column list.
fn col_perm(cols: &[usize], n: usize) -> PermOwned {
    PermOwned::new(cols.iter().enumerate().fold(vec![0usize; n], |mut v, (new_j, &old_j)| {
        v[old_j] = new_j; v
    }))
}

// Extracts U: the first `rank` rows of the reduced matrix with columns reordered by `cols`.
fn build_u<R>(work: &Mat<R>, cols: &[usize], rank: usize) -> Mat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    Mat::from_generator((rank, work.ncols()), |i, j| work[(i, cols[j])].clone())
}

// Builds the Schur complement s: (m-rank)×(n-rank), the bottom-right non-pivot block.
// Satisfies p*A*q = L*U + [[0,0],[0,s]].
fn build_s<R>(work: &Mat<R>, cols: &[usize], rank: usize) -> Mat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (m, n) = work.shape();
    Mat::from_generator((m - rank, n - rank), |i, j| work[(i + rank, cols[j + rank])].clone())
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use yui_core::num::Ratio;

    use super::*;

    type R = Ratio<i64>;

    fn r(n: i64) -> R { R::from(n) }
    fn rf(n: i64, d: i64) -> R { R::new(n, d) }

    fn sample() -> Mat<R> {
        Mat::from_data((3, 4), [
            r(1), r(2), r(3), r(4),
            r(2), r(4), r(5), r(6),
            r(3), r(6), r(7), r(8),
        ])
    }

    // Applies permutations to compute p_mat * a * q_mat as a plain matrix.
    fn apply_perms(a: &Mat<R>, p: &PermOwned, q: &PermOwned) -> Mat<R> {
        let (m, n) = a.shape();
        let mut out = Mat::zero((m, n));
        for i in 0..m {
            for j in 0..n {
                out[(p.at(i), q.at(j))] = a[(i, j)].clone();
            }
        }
        out
    }

    // Checks structural invariants and the main decomposition identity.
    // Returns the Pluq so callers can assert field-specific properties (s.is_zero, rank, …).
    fn check(a: &Mat<R>) -> Pluq<R> {
        let (m, n) = a.shape();
        let pp = pluq(a);
        let rank = pp.rank();

        assert_eq!(pp.l.shape(), (m, rank),         "L shape");
        assert_eq!(pp.u.shape(), (rank, n),          "U shape");
        assert_eq!(pp.s.shape(), (m - rank, n - rank), "s shape");

        // L is unit lower triangular: 1s on diagonal, 0 above diagonal
        for k in 0..rank {
            assert_eq!(pp.l[(k, k)], r(1), "L[{k},{k}] should be 1");
            for i in 0..k {
                assert_eq!(pp.l[(i, k)], r(0), "L[{i},{k}] should be 0 (above diagonal)");
            }
        }

        // U is upper echelon: 0 below diagonal in the rank×rank block
        for k in 0..rank {
            for i in (k + 1)..rank {
                assert_eq!(pp.u[(i, k)], r(0), "U[{i},{k}] should be 0 (below diagonal)");
            }
        }

        // Main invariant: p_mat * A * q_mat = L * U + [[0,0],[0,s]]
        let paq = apply_perms(a, &pp.p, &pp.q);
        let rem_full = Mat::from_generator((m, n), |i, j| {
            if i >= rank && j >= rank { pp.s[(i - rank, j - rank)].clone() } else { R::zero() }
        });
        assert_eq!(paq, &pp.l * &pp.u + &rem_full, "p*A*q should equal L*U + s");

        pp
    }

    #[test]
    fn test_sample() {
        let pp = check(&sample());
        assert_eq!(pp.rank(), 2);
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_zero() {
        let pp = check(&Mat::<R>::zero((3, 4)));
        assert_eq!(pp.rank(), 0);
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_identity() {
        let pp = check(&Mat::id(3));
        assert_eq!(pp.rank(), 3);
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_full_row_rank() {
        let a = Mat::from_data((2, 3), [
            r(1), r(0), r(2),
            r(0), r(1), r(3),
        ]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 2);
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_full_col_rank() {
        let a = Mat::from_data((3, 2), [
            r(1), r(2),
            r(3), r(4),
            r(5), r(6),
        ]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 2);
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_rank_deficient_cols() {
        // Column 2 = 2 * column 0
        let a = Mat::from_data((3, 3), [
            r(1), r(0), r(2),
            r(2), r(1), r(4),
            r(3), r(2), r(6),
        ]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 2);
        assert_eq!(pp.q.at(0), 0);
        assert_eq!(pp.q.at(1), 1);
        assert_eq!(pp.q.at(2), 2); // col 2 is non-pivot
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_pivot_not_in_first_col() {
        // First column is all zeros
        let a = Mat::from_data((2, 3), [
            r(0), r(1), r(2),
            r(0), r(3), r(4),
        ]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 2);
        assert!(pp.q.at(0) >= pp.rank(), "col 0 is non-pivot");
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_row_swap() {
        // First row has no unit in col 0
        let a = Mat::from_data((3, 3), [
            r(0), r(1), r(2),
            r(1), r(0), r(3),
            r(2), r(1), r(4),
        ]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 3);
        assert_eq!(pp.p.at(1), 0, "original row 1 should move to position 0");
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_fractions() {
        let a = Mat::from_data((2, 2), [
            rf(1, 2), rf(1, 3),
            rf(1, 4), rf(1, 5),
        ]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 2);
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_single_row() {
        let a = Mat::from_data((1, 4), [r(0), r(2), r(0), r(3)]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 1);
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_single_col() {
        let a = Mat::from_data((3, 1), [r(2), r(0), r(4)]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 1);
        assert!(pp.s.is_zero());
    }

    // Over a ring (Z), only ±1 are units.  Check that s is non-zero when
    // the matrix has non-unit entries that cannot be fully eliminated.
    #[test]
    fn test_ring_nonzero_rem() {
        // Only the (1,1) entry is ±1; the rest are non-units in Z.
        let a = Mat::<i32>::from_data((2, 2), [2, 3, 1, 4]);
        let pp = pluq(&a);

        // rank 1: only the unit entry (1) at row 1, col 0 becomes a pivot
        assert_eq!(pp.rank(), 1);
        assert_eq!(pp.l.shape(), (2, 1));
        assert_eq!(pp.u.shape(), (1, 2));
        assert_eq!(pp.s.shape(), (1, 1)); // (m-rank, n-rank) = (1, 1)

        // Invariant holds over Z
        let paq: Mat<i32> = {
            let (m, n) = a.shape();
            let mut out = Mat::zero((m, n));
            for i in 0..m { for j in 0..n { out[(pp.p.at(i), pp.q.at(j))] = a[(i, j)]; } }
            out
        };
        let rem_full = Mat::from_generator((2, 2), |i, j| {
            if i >= 1 && j >= 1 { pp.s[(i - 1, j - 1)] } else { 0 }
        });
        assert_eq!(paq, &pp.l * &pp.u + &rem_full);

        // s is non-zero (2 is not a unit in Z)
        assert!(!pp.s.is_zero());
    }

    // ---- solve_pluq tests ----

    fn solve_check(a: &Mat<R>, y: &[R]) -> Vec<R> {
        let x = solve_pluq(a, y).expect("expected a solution");
        // verify A * x = y
        let (m, n) = a.shape();
        assert_eq!(x.len(), n);
        for i in 0..m {
            let ax_i: R = (0..n).fold(R::zero(), |acc, j| acc + &a[(i, j)] * &x[j]);
            assert_eq!(ax_i, y[i], "row {i}: (A*x)[{i}] != y[{i}]");
        }
        x
    }

    #[test]
    fn test_solve_square_full_rank() {
        // 2×2 invertible matrix
        let a = Mat::from_data((2, 2), [r(1), r(2), r(3), r(4)]);
        let y = vec![r(5), r(6)];
        solve_check(&a, &y);
    }

    #[test]
    fn test_solve_overdetermined_consistent() {
        // 3×2 matrix, consistent y
        let a = Mat::from_data((3, 2), [r(1), r(0), r(0), r(1), r(1), r(1)]);
        let y = vec![r(2), r(3), r(5)]; // y = a * [2, 3]
        solve_check(&a, &y);
    }

    #[test]
    fn test_solve_overdetermined_inconsistent() {
        let a = Mat::from_data((3, 2), [r(1), r(0), r(0), r(1), r(1), r(1)]);
        let y = vec![r(1), r(1), r(0)]; // 1+1 != 0, inconsistent
        assert!(solve_pluq(&a, &y).is_none());
    }

    #[test]
    fn test_solve_underdetermined() {
        // 2×3 matrix, rank 2; infinitely many solutions — we just get one
        let a = Mat::from_data((2, 3), [r(1), r(0), r(2), r(0), r(1), r(3)]);
        let y = vec![r(4), r(5)];
        solve_check(&a, &y);
    }

    #[test]
    fn test_solve_zero_rhs() {
        let a = Mat::from_data((2, 2), [r(1), r(2), r(3), r(4)]);
        let y = vec![r(0), r(0)];
        let x = solve_check(&a, &y);
        assert_eq!(x, vec![r(0), r(0)]);
    }

    #[test]
    fn test_solve_no_solution_rank_deficient() {
        // rank-1 matrix; y not in column space
        let a = Mat::from_data((2, 2), [r(1), r(2), r(2), r(4)]);
        let y = vec![r(1), r(0)]; // not in column space
        assert!(solve_pluq(&a, &y).is_none());
    }
}
