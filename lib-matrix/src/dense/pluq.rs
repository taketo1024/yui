// Implemented with the help of Claude Code. 

use sprs::PermOwned;
use yui_core::{Ring, RingOps, Field, FieldOps};
use crate::MatTrait;
use crate::dense::Mat;

/// Result of a PLUQ decomposition satisfying `p_mat * A * q_mat = L * U + rem`:
///   - `p`: row permutation (pivot rows first)
///   - `q`: column permutation (pivot columns first)
///   - `l`: `m × rank`, unit lower triangular — elimination multipliers
///   - `u`: `rank × n`, upper echelon — the reduced pivot rows
///   - `rem`: `m × n`, remainder — zero when `R` is a field
pub struct Pluq<R> {
    pub p: PermOwned,
    pub q: PermOwned,
    pub l: Mat<R>,
    pub u: Mat<R>,
    pub rem: Mat<R>,
}

impl<R> Pluq<R> {
    pub fn rank(&self) -> usize { self.l.ncols() }
}

/// Computes a PLUQ decomposition of `a` over a ring by Gaussian elimination.
///
/// Only unit elements are used as pivots, so every elimination step is exact.
/// Satisfies `p_mat * a * q_mat = l * u + rem`, where `rem` is zero when `R`
/// is a field (every non-zero element is a unit).
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
    let u = build_u(&work, &cols, rank, n);
    let rem = build_rem(&work, &cols, rank, m, n);

    Pluq { p, q, l, u, rem }
}

/// Solves `A * x = y` over a field using PLUQ decomposition.
///
/// Returns `Some(x)` if a solution exists, `None` otherwise.
pub fn solve_pluq<R>(a: &Mat<R>, y: &[R]) -> Option<Vec<R>>
where R: Field, for<'x> &'x R: FieldOps<R> {
    let (m, n) = a.shape();
    assert_eq!(y.len(), m);

    let pp = pluq(a);
    let rank = pp.rank();
    let Pluq { p, q, l, u, .. } = pp;

    // y' = P * y
    let mut yp = vec![R::zero(); m];
    for i in 0..m {
        yp[p.at(i)] = y[i].clone();
    }

    // Forward substitution: solve L * z = y' (L is unit lower triangular, m × rank)
    let mut z = vec![R::zero(); rank];
    for k in 0..rank {
        let mut val = yp[k].clone();
        for j in 0..k {
            val = val - &l[(k, j)] * &z[j];
        }
        z[k] = val; // L[k,k] = 1
    }

    // Consistency check for non-pivot rows
    for i in rank..m {
        let mut lhs = R::zero();
        for j in 0..rank {
            lhs = lhs + &l[(i, j)] * &z[j];
        }
        if lhs != yp[i] {
            return None;
        }
    }

    // Back substitution: solve U * x' = z, free variables x'[rank..n] = 0
    let mut xp = vec![R::zero(); n];
    for k in (0..rank).rev() {
        let mut val = z[k].clone();
        for j in (k + 1)..rank {
            val = val - &u[(k, j)] * &xp[j];
        }
        xp[k] = val * u[(k, k)].inv().unwrap();
    }

    // Recover x: x = Q * x', i.e., x[j] = x'[q.at(j)]
    let mut x = vec![R::zero(); n];
    for j in 0..n {
        x[j] = xp[q.at(j)].clone();
    }

    Some(x)
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
            for col in l_cols.iter_mut() { col.swap(r, pivot_pos); }
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
    let m = work.nrows();
    let pivot_inv = work[(r, j)].inv().unwrap();
    let mut col = vec![R::zero(); m];
    col[r] = R::one();
    for i in (r + 1)..m {
        col[i] = work[(i, j)].clone() * pivot_inv.clone();
    }
    col
}

// Subtracts `l_col[i] * row_r` from each row `i > r`, zeroing out the pivot column below.
fn eliminate_below<R>(work: &mut Mat<R>, l_col: &[R], r: usize)
where R: Ring, for<'x> &'x R: RingOps<R> {
    let m = work.nrows();
    for i in (r + 1)..m {
        if l_col[i].is_zero() { continue; }
        work.add_row_to(r, i, &-l_col[i].clone());
    }
}

// Builds the row permutation: p.at(orig) = current position of that row.
fn row_perm(row_of: &[usize], m: usize) -> PermOwned {
    let mut fwd = vec![0usize; m];
    for (pos, &orig) in row_of.iter().enumerate() {
        fwd[orig] = pos;
    }
    PermOwned::new(fwd)
}

// Returns the full column reordering: pivot columns first, non-pivot columns last.
fn col_order(pivot_cols: &[usize], n: usize) -> Vec<usize> {
    let mut is_pivot = vec![false; n];
    for &j in pivot_cols { is_pivot[j] = true; }
    let non_pivot: Vec<usize> = (0..n).filter(|&j| !is_pivot[j]).collect();
    pivot_cols.iter().chain(non_pivot.iter()).cloned().collect()
}

// Builds the column permutation from the full ordered column list.
fn col_perm(cols: &[usize], n: usize) -> PermOwned {
    let mut fwd = vec![0usize; n];
    for (new_j, &old_j) in cols.iter().enumerate() {
        fwd[old_j] = new_j;
    }
    PermOwned::new(fwd)
}

// Extracts U: the first `rank` rows of the reduced matrix with columns reordered by `cols`.
fn build_u<R>(work: &Mat<R>, cols: &[usize], rank: usize, n: usize) -> Mat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    Mat::from_generator((rank, n), |i, j| work[(i, cols[j])].clone())
}

// Builds rem: m×n, with zero rows for the pivot rows and the column-reordered
// remaining rows for the non-pivot rows.  Satisfies p*A*q = L*U + rem.
fn build_rem<R>(work: &Mat<R>, cols: &[usize], rank: usize, m: usize, n: usize) -> Mat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    Mat::from_generator((m, n), |i, j| {
        if i < rank { R::zero() }
        else { work[(i, cols[j])].clone() }
    })
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
    // Returns the Pluq so callers can assert field-specific properties (rem.is_zero, rank, …).
    fn check(a: &Mat<R>) -> Pluq<R> {
        let (m, n) = a.shape();
        let pp = pluq(a);
        let rank = pp.rank();

        assert_eq!(pp.l.shape(),   (m, rank), "L shape");
        assert_eq!(pp.u.shape(),   (rank, n), "U shape");
        assert_eq!(pp.rem.shape(), (m, n),    "rem shape");

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

        // rem has zero rows in the pivot positions
        for k in 0..rank {
            for j in 0..n {
                assert_eq!(pp.rem[(k, j)], r(0), "rem[{k},{j}] should be 0 (pivot row)");
            }
        }

        // Main invariant: p_mat * A * q_mat = L * U + rem
        let paq = apply_perms(a, &pp.p, &pp.q);
        assert_eq!(paq, &pp.l * &pp.u + &pp.rem, "p*A*q should equal L*U + rem");

        pp
    }

    #[test]
    fn test_sample() {
        let pp = check(&sample());
        assert_eq!(pp.rank(), 2);
        assert!(pp.rem.is_zero());
    }

    #[test]
    fn test_zero() {
        let pp = check(&Mat::<R>::zero((3, 4)));
        assert_eq!(pp.rank(), 0);
        assert!(pp.rem.is_zero());
    }

    #[test]
    fn test_identity() {
        let pp = check(&Mat::id(3));
        assert_eq!(pp.rank(), 3);
        assert!(pp.rem.is_zero());
    }

    #[test]
    fn test_full_row_rank() {
        let a = Mat::from_data((2, 3), [
            r(1), r(0), r(2),
            r(0), r(1), r(3),
        ]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 2);
        assert!(pp.rem.is_zero());
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
        assert!(pp.rem.is_zero());
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
        assert!(pp.rem.is_zero());
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
        assert!(pp.rem.is_zero());
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
        assert!(pp.rem.is_zero());
    }

    #[test]
    fn test_fractions() {
        let a = Mat::from_data((2, 2), [
            rf(1, 2), rf(1, 3),
            rf(1, 4), rf(1, 5),
        ]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 2);
        assert!(pp.rem.is_zero());
    }

    #[test]
    fn test_single_row() {
        let a = Mat::from_data((1, 4), [r(0), r(2), r(0), r(3)]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 1);
        assert!(pp.rem.is_zero());
    }

    #[test]
    fn test_single_col() {
        let a = Mat::from_data((3, 1), [r(2), r(0), r(4)]);
        let pp = check(&a);
        assert_eq!(pp.rank(), 1);
        assert!(pp.rem.is_zero());
    }

    // Over a ring (Z), only ±1 are units.  Check that rem is non-zero when
    // the matrix has non-unit entries that cannot be fully eliminated.
    #[test]
    fn test_ring_nonzero_rem() {
        // Only the (1,1) entry is ±1; the rest are non-units in Z.
        let a = Mat::<i32>::from_data((2, 2), [2, 3, 1, 4]);
        let pp = pluq(&a);

        // rank 1: only the unit entry (1) at row 1, col 0 becomes a pivot
        assert_eq!(pp.rank(), 1);
        assert_eq!(pp.l.shape(),   (2, 1));
        assert_eq!(pp.u.shape(),   (1, 2));
        assert_eq!(pp.rem.shape(), (2, 2));

        // Invariant holds over Z
        let paq: Mat<i32> = {
            let (m, n) = a.shape();
            let mut out = Mat::zero((m, n));
            for i in 0..m { for j in 0..n { out[(pp.p.at(i), pp.q.at(j))] = a[(i, j)]; } }
            out
        };
        assert_eq!(paq, &pp.l * &pp.u + &pp.rem);

        // rem is non-zero (2 is not a unit in Z)
        assert!(!pp.rem.is_zero());
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
