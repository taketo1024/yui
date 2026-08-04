// Partial sparse PLUQ decomposition.
// Implemented with the help of Claude Code. 
//
// Reference:
//   "Parallel Sparse PLUQ Factorization modulo p", Bouillaguet–Delaplace–Voge.
//   https://hal.inria.fr/hal-01646133/document

use sprs::PermOwned;
use yui_core::{Ring, RingOps};

use crate::MatTrait;
use super::SpMat;
use super::pivot::{PivotCondition, PivotType, find_pivots, perms_by_pivots};

/// Result of a partial PLUQ decomposition.
///
/// For `PivotType::Rows`:  `p * A * q = piv.stack(&rem)`
///   - `piv` has `rank` rows; each column `k < rank` has its pivot at row `k`.
///   - `rem` has `m - rank` rows (the non-pivot rows, permuted).
///
/// For `PivotType::Cols`:  `p * A * q = piv.concat(&rem)`
///   - `piv` has `rank` cols; each row `k < rank` has its pivot at column `k`.
///   - `rem` has `n - rank` cols (the non-pivot columns, permuted).
pub struct PartialPluq<R> {
    pub p: PermOwned,
    pub q: PermOwned,
    pub rank: usize,
    pub piv: SpMat<R>,
    pub rem: SpMat<R>,
}

/// Computes a partial PLUQ decomposition of `a` over a field.
///
/// When `recursive` is false (single-shot), calls `PivotFinder` once and splits
/// `p * a * q` at the pivot rank without performing any arithmetic.
///
/// When `recursive` is true, at each step extracts the non-pivot sub-block
/// (rows and columns beyond the current rank in the permuted space) and repeats
/// until `find_pivots` returns no new pivots, composing the permutations throughout.
pub fn pre_pluq<R>(a: &SpMat<R>, piv_type: PivotType, piv_cond: PivotCondition, recursive: bool) -> PartialPluq<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (m, n) = a.shape();

    let mut cur_p = PermOwned::identity(m);
    let mut cur_q = PermOwned::identity(n);
    let mut total_r = 0;
    let mut cur_mat: Option<SpMat<R>> = None;

    loop {
        let mat: &SpMat<R> = cur_mat.as_ref().unwrap_or(a);
        let pivots = find_pivots(mat, piv_type, piv_cond);
        
        if pivots.is_empty() { break; }

        let r = pivots.len();
        let (p, q) = perms_by_pivots(mat, &pivots);

        cur_p = compose_perms(&cur_p, total_r, &p, m);
        cur_q = compose_perms(&cur_q, total_r, &q, n);
        total_r += r;
        
        let (cm, cn) = mat.shape();

        if !recursive || r >= cm || r >= cn { break; }

        cur_mat = Some(permuted_sub(mat, &p, &q, r));
    }

    let (piv, rem) = split(a, piv_type, &cur_p, &cur_q, total_r);

    PartialPluq { p: cur_p, q: cur_q, rank: total_r, piv, rem }
}

fn split<R>(a: &SpMat<R>, piv_type: PivotType, p: &PermOwned, q: &PermOwned, r: usize) -> (SpMat<R>, SpMat<R>)
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (m, n) = a.shape();
    match piv_type {
        PivotType::Rows => {
            // p * a * q = [piv; rem]  (piv: r×n, rem: (m-r)×n)
            let mut piv_ent = vec![];
            let mut rem_ent = vec![];
            for (i, j, v) in a.iter() {
                let pi = p.at(i);
                let qj = q.at(j);
                if pi < r {
                    piv_ent.push((pi, qj, v.clone()));
                } else {
                    rem_ent.push((pi - r, qj, v.clone()));
                }
            }
            (
                SpMat::from_entries((r, n), piv_ent),
                SpMat::from_entries((m - r, n), rem_ent),
            )
        }
        PivotType::Cols => {
            // p * a * q = [piv | rem]  (piv: m×r, rem: m×(n-r))
            let mut piv_ent = vec![];
            let mut rem_ent = vec![];
            for (i, j, v) in a.iter() {
                let pi = p.at(i);
                let qj = q.at(j);
                if qj < r {
                    piv_ent.push((pi, qj, v.clone()));
                } else {
                    rem_ent.push((pi, qj - r, v.clone()));
                }
            }
            (
                SpMat::from_entries((m, r), piv_ent),
                SpMat::from_entries((m, n - r), rem_ent),
            )
        }
    }
}

// Applies (p, q) to `mat` and extracts the sub-block landing at positions >= r,
// all in a single pass — equivalent to mat.permute(p,q).submat(r..m, r..n).
fn permuted_sub<R>(mat: &SpMat<R>, p: &PermOwned, q: &PermOwned, r: usize) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (m, n) = mat.shape();
    mat.extract((m - r, n - r), |i, j| {
        let pi = p.at(i);
        let qj = q.at(j);
        if pi >= r && qj >= r { Some((pi - r, qj - r)) } else { None }
    })
}

// Composes two permutations: outer (size `total`) then inner (size `total - split`)
// applied to the tail.  new[i] = outer[i] if < split, else split + inner[outer[i] - split].
fn compose_perms(outer: &PermOwned, split: usize, inner: &PermOwned, total: usize) -> PermOwned {
    let fwd: Vec<usize> = (0..total).map(|i| {
        let oi = outer.at(i);
        if oi < split { oi } else { split + inner.at(oi - split) }
    }).collect();
    PermOwned::new(fwd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::One;

    // The 6×9 matrix used in pivot.rs tests — rank 5 for Rows, rank 6 for Cols.
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

    fn check_rows(a: &SpMat<i32>) {
        let pp = pre_pluq(a, PivotType::Rows, PivotCondition::AnyUnit, false);
        let (m, n) = a.shape();
        let r = pp.rank;

        assert_eq!(pp.piv.shape(), (r, n));
        assert_eq!(pp.rem.shape(), (m - r, n));

        // p * a * q = piv.stack(&rem)
        let stacked = pp.piv.stack(&pp.rem);
        assert_eq!(stacked, a.permute(pp.p.view(), pp.q.view()));

        // The first r×r block has pivots on the diagonal and zeros below it
        // (ensured by the top-sort in PivotFinder::result).
        let b = pp.piv.clone().into_dense();
        for k in 0..r {
            assert!(b[(k, k)].is_one(), "b[{k},{k}] should be a pivot (=1)");
        }
        for j in 0..r {
            for i in j + 1..r {
                assert_eq!(b[(i, j)], 0, "b[{i},{j}] should be zero (below diagonal)");
            }
        }
    }

    fn check_cols(a: &SpMat<i32>) {
        let pp = pre_pluq(a, PivotType::Cols, PivotCondition::AnyUnit, false);
        let (m, n) = a.shape();
        let r = pp.rank;

        assert_eq!(pp.piv.shape(), (m, r));
        assert_eq!(pp.rem.shape(), (m, n - r));

        // p * a * q = piv.concat(&rem)
        let concatenated = pp.piv.concat(&pp.rem);
        assert_eq!(concatenated, a.permute(pp.p.view(), pp.q.view()));

        // The first r×r block has pivots on the diagonal and zeros above it
        let b = pp.piv.clone().into_dense();
        for k in 0..r {
            assert!(b[(k, k)].is_one(), "b[{k},{k}] should be a pivot (=1)");
        }
        for i in 0..r {
            for j in i + 1..r {
                assert_eq!(b[(i, j)], 0, "b[{i},{j}] should be zero (above diagonal)");
            }
        }
    }

    // Recursive variants: verify shape, reconstruction, and diagonal pivots.
    // (The block-triangular sub-structure within each pass is an internal detail.)
    fn check_rows_recursive(a: &SpMat<i32>) {
        let r0 = pre_pluq(a, PivotType::Rows, PivotCondition::AnyUnit, false).rank;
        let pp = pre_pluq(a, PivotType::Rows, PivotCondition::AnyUnit, true);
        let (m, n) = a.shape();
        let r = pp.rank;

        assert!(r >= r0);
        assert_eq!(pp.piv.shape(), (r, n));
        assert_eq!(pp.rem.shape(), (m - r, n));

        let stacked = pp.piv.stack(&pp.rem);
        assert_eq!(stacked, a.permute(pp.p.view(), pp.q.view()));

        let b = pp.piv.clone().into_dense();
        for k in 0..r {
            assert!(b[(k, k)].is_one(), "b[{k},{k}] should be a pivot (=1)");
        }
    }

    fn check_cols_recursive(a: &SpMat<i32>) {
        let r0 = pre_pluq(a, PivotType::Cols, PivotCondition::AnyUnit, false).rank;
        let pp = pre_pluq(a, PivotType::Cols, PivotCondition::AnyUnit, true);
        let (m, n) = a.shape();
        let r = pp.rank;

        assert!(r >= r0);
        assert_eq!(pp.piv.shape(), (m, r));
        assert_eq!(pp.rem.shape(), (m, n - r));

        let concatenated = pp.piv.concat(&pp.rem);
        assert_eq!(concatenated, a.permute(pp.p.view(), pp.q.view()));

        let b = pp.piv.clone().into_dense();
        for k in 0..r {
            assert!(b[(k, k)].is_one(), "b[{k},{k}] should be a pivot (=1)");
        }
    }

    #[test]
    fn test_rows() {
        check_rows(&sample());
    }

    #[test]
    fn test_cols() {
        check_cols(&sample());
    }

    #[test]
    fn test_zero() {
        let a = SpMat::<i32>::zero((4, 5));
        let pp = pre_pluq(&a, PivotType::Rows, PivotCondition::AnyUnit, false);
        assert_eq!(pp.rank, 0);
        assert_eq!(pp.piv.shape(), (0, 5));
        assert_eq!(pp.rem.shape(), (4, 5));
        assert_eq!(pp.rem, a.permute(pp.p.view(), pp.q.view()));
    }

    #[test]
    fn test_square_full_rank() {
        let a = SpMat::from_dense_data((3, 3), [
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
        ]);
        let pp = pre_pluq(&a, PivotType::Rows, PivotCondition::AnyUnit, false);
        assert_eq!(pp.rank, 3);
        assert_eq!(pp.rem.shape(), (0, 3));
    }

    #[test]
    fn test_rank_rows() {
        assert_eq!(pre_pluq(&sample(), PivotType::Rows, PivotCondition::AnyUnit, false).rank, 5);
    }

    #[test]
    fn test_rank_cols() {
        assert_eq!(pre_pluq(&sample(), PivotType::Cols, PivotCondition::AnyUnit, false).rank, 6);
    }

    #[test]
    fn test_rand_rows() {
        let a = SpMat::<i32>::rand((40, 60), 0.1);
        check_rows(&a);
    }

    #[test]
    fn test_rand_cols() {
        let a = SpMat::<i32>::rand((40, 60), 0.1);
        check_cols(&a);
    }

    #[test]
    fn test_recursive_rows() {
        check_rows_recursive(&sample());
    }

    #[test]
    fn test_recursive_cols() {
        check_cols_recursive(&sample());
    }

    #[test]
    fn test_rand_recursive_rows() {
        let a = SpMat::<i32>::rand((40, 60), 0.1);
        check_rows_recursive(&a);
    }

    #[test]
    fn test_rand_recursive_cols() {
        let a = SpMat::<i32>::rand((40, 60), 0.1);
        check_cols_recursive(&a);
    }
}