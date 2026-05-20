// Implemented with the help of Claude Code.

use log::debug;
use nalgebra::Scalar;
use yui_core::{Ring, RingOps, Field, FieldOps};
use crate::{MatTrait, Perm};
use crate::dense::Mat;

/// Result of a PLUQ decomposition satisfying `p_mat * A * q_mat = L * U + s`:
///   - `p`: row permutation (pivot rows first)
///   - `q`: column permutation (pivot columns first)
///   - `l`: `m × rank`, lower triangular — pivot values on diagonal
///   - `u`: `rank × n`, unit upper triangular — elimination multipliers
///   - `s`: `(m - rank) × (n - rank)`, Schur complement — zero when `R` is a field
pub struct Pluq<R> {
    pub p: Perm,
    pub q: Perm,
    pub l: Mat<R>,
    pub u: Mat<R>,
    pub s: Mat<R>,
}

impl<R> Pluq<R> {
    pub fn rank(&self) -> usize { self.l.n_cols() }
}

impl<R: Scalar> Pluq<R> {
    /// Returns the PLUQ decomposition of `A^T`.
    /// If `self` satisfies `p * A * q = l * u + rest` (Cols convention),
    /// then `self.transpose()` satisfies `q * A^T * p = u^T * l^T + rest^T` (Rows convention).
    pub fn transpose(self) -> Self {
        let Pluq { p, q, l, u, s } = self;
        Pluq { p: q, q: p, l: u.transpose(), u: l.transpose(), s: s.transpose() }
    }
}

/// Computes a PLUQ decomposition of `a` over a ring by Gaussian elimination.
///
/// Only unit elements are used as pivots, so every elimination step is exact.
/// Satisfies `p_mat * a * q_mat = l * u + s`, where `s` is the Schur complement
/// and is zero when `R` is a field (every non-zero element is a unit).
///
/// The elimination uses column operations (zeroing out each pivot row to the right),
/// which aligns with `DMatrix`'s column-major storage and avoids strided row access.
pub fn pluq<R>(a: &Mat<R>) -> Pluq<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    debug!("compute dense pluq: {:?}", a.shape());

    let (m, n) = a.shape();
    let mut work = a.clone();
    let mut col_of: Vec<usize> = (0..n).collect();

    let (pivot_rows, u) = reduce(&mut work, &mut col_of);
    let rank = pivot_rows.len();
    let p = Perm::forward_indices(m, pivot_rows.iter().copied());
    let q = Perm::from_indices(col_of).inv();
    let p_inv = p.inv();
    let l = build_l(&work, &p_inv, rank);
    let s = build_s(&work, &p_inv, rank);

    Pluq { p, q, l, u, s }
}

/// Solves `A * x = y` over a field using PLUQ decomposition.
/// Uses column-based elimination; see [`pluq`] for details.
///
/// Returns `Some(x)` if a solution exists, `None` otherwise.
pub fn solve_pluq<R>(a: &Mat<R>, y: &[R]) -> Option<Vec<R>>
where R: Field, for<'x> &'x R: FieldOps<R> {
    debug!("dense solve: {:?}", a.shape());

    assert_eq!(y.len(), a.n_rows());
    let Pluq { p, q, l, u, .. } = pluq(a);
    let yp = p.apply_to(y.to_vec());

    debug!("forward sub: {:?}", l.shape());

    let z  = forward_sub(&l, &yp);
    if !check_consistent(&l, &yp, &z) { return None; }

    debug!("back sub: {:?}", u.shape());

    let xp = back_sub(&u, &z);
    Some((0..xp.len()).map(|j| xp[q.at(j)].clone()).collect())
}

// Solves L * z = yp[0..rank] by forward substitution (L is lower triangular, pivot values on diagonal).
fn forward_sub<R>(l: &Mat<R>, yp: &[R]) -> Vec<R>
where R: Field, for<'x> &'x R: FieldOps<R> {
    (0..l.n_cols()).fold(vec![], |mut z, k| {
        let pivot_inv = l[(k, k)].inv().unwrap();
        let val = (0..k).fold(yp[k].clone(), |v, j| v - &l[(k, j)] * &z[j]) * pivot_inv;
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

// Solves U * xp = z by back substitution (U is unit upper triangular); free variables xp[rank..n] stay zero.
fn back_sub<R>(u: &Mat<R>, z: &[R]) -> Vec<R>
where R: Field, for<'x> &'x R: FieldOps<R> {
    let (rank, n) = (z.len(), u.n_cols());
    (0..rank).rev().fold(vec![R::zero(); n], |mut xp, k| {
        xp[k] = (k + 1..rank).fold(z[k].clone(), |v, j| v - &u[(k, j)] * &xp[j]);
        xp
    })
}

// Gaussian elimination in place, eliminating only right of each pivot.
// Only unit elements are accepted as pivots so that `inv()` is always valid.
// Builds U row by row: each pivot step contributes one row
// (unit on the diagonal, multipliers to the right, zeros to the left).
// Any column swaps are mirrored in the already-built U rows.
//
// Returns (pivot_row_indices, U).
fn reduce<R>(work: &mut Mat<R>, col_of: &mut Vec<usize>) -> (Vec<usize>, Mat<R>)
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (m, n) = work.shape();
    let mut pivot_rows = Vec::new();
    let mut u_rows: Vec<Vec<R>> = Vec::new();
    let mut c = 0;

    for i in 0..m {
        if c >= n { break; }
        let Some(pivot_pos) = (c..n).find(|&j| work[(i, j)].is_unit()) else { continue; };
        if pivot_pos != c {
            work.swap_cols(c, pivot_pos);
            col_of.swap(c, pivot_pos);
            u_rows.iter_mut().for_each(|row| row.swap(c, pivot_pos));
        }
        let u_row = build_u_row(work, i, c);
        eliminate_right(work, &u_row, c);
        u_rows.push(u_row);
        pivot_rows.push(i);
        c += 1;
    }

    let rank = pivot_rows.len();
    let u = Mat::generate((rank, n), |k, j| u_rows[k][j].clone());
    (pivot_rows, u)
}

// Builds row i of U: 1 on the diagonal, multipliers (entry * pivot_inv) to the right, 0 to the left.
fn build_u_row<R>(work: &Mat<R>, i: usize, c: usize) -> Vec<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    use std::cmp::Ordering::*;
    let pivot_inv = work[(i, c)].inv().unwrap();
    (0..work.n_cols()).map(|j| match j.cmp(&c) {
        Less    => R::zero(),
        Equal   => R::one(),
        Greater => work[(i, j)].clone() * pivot_inv.clone(),
    }).collect()
}

// Subtracts `u_row[j] * col_c` from each column `j > c`, zeroing out the pivot row to the right.
fn eliminate_right<R>(work: &mut Mat<R>, u_row: &[R], c: usize)
where R: Ring, for<'x> &'x R: RingOps<R> {
    let n = work.n_cols();
    (c + 1..n)
        .filter(|&j| !u_row[j].is_zero())
        .for_each(|j| work.add_col_to(c, j, &-u_row[j].clone()));
}

// Extracts L: the first `rank` columns of the reduced matrix with rows reordered by `p_inv`.
fn build_l<R>(work: &Mat<R>, p_inv: &Perm, rank: usize) -> Mat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    Mat::generate((work.n_rows(), rank), |i, k| work[(p_inv.at(i), k)].clone())
}

// Builds the Schur complement s: (m-rank)×(n-rank), the bottom-right non-pivot block.
// Satisfies p*A*q = L*U + [[0,0],[0,s]].
fn build_s<R>(work: &Mat<R>, p_inv: &Perm, rank: usize) -> Mat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (m, n) = work.shape();
    Mat::generate((m - rank, n - rank), |i, j| work[(p_inv.at(i + rank), rank + j)].clone())
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use yui_core::num::Ratio;

    use super::*;

    // ---- Pluq::transpose tests ----

    #[test]
    fn test_pluq_transpose() {
        // a = [[1,2,3],[4,5,6]], rank 2.
        // pluq(a) factors a; pluq(a).transpose() should factor a^T.
        type R = Ratio<i64>;
        let r = |n: i64| R::from(n);

        let a = Mat::from_row_major((2, 3), [r(1),r(2),r(3),r(4),r(5),r(6)]);
        let at = a.transpose(); // 3×2

        let dp = pluq(&a).transpose();
        let rank = dp.rank();

        let (m, n) = at.shape(); // (3, 2)
        assert_eq!(dp.l.shape(), (m, rank));
        assert_eq!(dp.u.shape(), (rank, n));

        // p * A^T * q = l * u + rest
        let paq = apply_perms(&at, &dp.p, &dp.q);
        let rem_full = Mat::generate((m, n), |i, j| {
            if i >= rank && j >= rank { dp.s[(i - rank, j - rank)].clone() } else { R::zero() }
        });
        assert_eq!(paq, &dp.l * &dp.u + &rem_full);
    }

    type R = Ratio<i64>;

    fn r(n: i64) -> R { R::from(n) }
    fn rf(n: i64, d: i64) -> R { R::new(n, d) }

    fn sample() -> Mat<R> {
        Mat::from_row_major((3, 4), [
            r(1), r(2), r(3), r(4),
            r(2), r(4), r(5), r(6),
            r(3), r(6), r(7), r(8),
        ])
    }

    // Applies permutations to compute p_mat * a * q_mat as a plain matrix.
    fn apply_perms(a: &Mat<R>, p: &Perm, q: &Perm) -> Mat<R> {
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

        // L is lower triangular: 0 above diagonal
        for k in 0..rank {
            for i in 0..k {
                assert_eq!(pp.l[(i, k)], r(0), "L[{i},{k}] should be 0 (above diagonal)");
            }
        }

        // U is unit upper triangular: 1s on diagonal, 0 below diagonal in the rank×rank block
        for k in 0..rank {
            assert_eq!(pp.u[(k, k)], r(1), "U[{k},{k}] should be 1");
            for i in (k + 1)..rank {
                assert_eq!(pp.u[(i, k)], r(0), "U[{i},{k}] should be 0 (below diagonal)");
            }
        }

        // Main invariant: p_mat * A * q_mat = L * U + [[0,0],[0,s]]
        let paq = apply_perms(a, &pp.p, &pp.q);
        let rem_full = Mat::generate((m, n), |i, j| {
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
        let a = Mat::from_row_major((2, 3), [
            r(1), r(0), r(2),
            r(0), r(1), r(3),
        ]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 2);
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_full_col_rank() {
        let a = Mat::from_row_major((3, 2), [
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
        let a = Mat::from_row_major((3, 3), [
            r(1), r(0), r(2),
            r(2), r(1), r(4),
            r(3), r(2), r(6),
        ]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 2);
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_pivot_not_in_first_col() {
        // First column is all zeros
        let a = Mat::from_row_major((2, 3), [
            r(0), r(1), r(2),
            r(0), r(3), r(4),
        ]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 2);
        assert!(pp.q.at(0) >= pp.rank(), "col 0 is non-pivot");
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_col_swap() {
        // First column has no unit in row 0
        let a = Mat::from_row_major((3, 3), [
            r(0), r(1), r(2),
            r(1), r(0), r(3),
            r(2), r(1), r(4),
        ]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 3);
        assert_eq!(pp.q.at(1), 0, "original col 1 should move to position 0");
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_fractions() {
        let a = Mat::from_row_major((2, 2), [
            rf(1, 2), rf(1, 3),
            rf(1, 4), rf(1, 5),
        ]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 2);
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_single_row() {
        let a = Mat::from_row_major((1, 4), [r(0), r(2), r(0), r(3)]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 1);
        assert!(pp.s.is_zero());
    }

    #[test]
    fn test_single_col() {
        let a = Mat::from_row_major((3, 1), [r(2), r(0), r(4)]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 1);
        assert!(pp.s.is_zero());
    }

    // Over a ring (Z), only ±1 are units.  Check that s is non-zero when
    // the matrix has non-unit entries that cannot be fully eliminated.
    #[test]
    fn test_ring_nonzero_rem() {
        // Row 0 has no units; row 1 col 0 has unit 1 → rank 1.
        let a = Mat::<i32>::from_row_major((2, 2), [2, 3, 1, 4]);
        let pp = pluq(&a);

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
        let rem_full = Mat::generate((2, 2), |i, j| {
            if i >= 1 && j >= 1 { pp.s[(i - 1, j - 1)] } else { 0 }
        });
        assert_eq!(paq, &pp.l * &pp.u + &rem_full);

        // s is non-zero (remaining entry is not a unit in Z)
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
        let a = Mat::from_row_major((2, 2), [r(1), r(2), r(3), r(4)]);
        let y = vec![r(5), r(6)];
        solve_check(&a, &y);
    }

    #[test]
    fn test_solve_overdetermined_consistent() {
        // 3×2 matrix, consistent y
        let a = Mat::from_row_major((3, 2), [r(1), r(0), r(0), r(1), r(1), r(1)]);
        let y = vec![r(2), r(3), r(5)]; // y = a * [2, 3]
        solve_check(&a, &y);
    }

    #[test]
    fn test_solve_overdetermined_inconsistent() {
        let a = Mat::from_row_major((3, 2), [r(1), r(0), r(0), r(1), r(1), r(1)]);
        let y = vec![r(1), r(1), r(0)]; // 1+1 != 0, inconsistent
        assert!(solve_pluq(&a, &y).is_none());
    }

    #[test]
    fn test_solve_underdetermined() {
        // 2×3 matrix, rank 2; infinitely many solutions — we just get one
        let a = Mat::from_row_major((2, 3), [r(1), r(0), r(2), r(0), r(1), r(3)]);
        let y = vec![r(4), r(5)];
        solve_check(&a, &y);
    }

    #[test]
    fn test_solve_zero_rhs() {
        let a = Mat::from_row_major((2, 2), [r(1), r(2), r(3), r(4)]);
        let y = vec![r(0), r(0)];
        let x = solve_check(&a, &y);
        assert_eq!(x, vec![r(0), r(0)]);
    }

    #[test]
    fn test_solve_no_solution_rank_deficient() {
        // rank-1 matrix; y not in column space
        let a = Mat::from_row_major((2, 2), [r(1), r(2), r(2), r(4)]);
        let y = vec![r(1), r(0)]; // not in column space
        assert!(solve_pluq(&a, &y).is_none());
    }
}
