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

/// Result of a partial PLUQ decomposition satisfying `p * A * q = l * u + rem`.
///
/// For `PivotType::Rows`:
///   - `l`: `m × rank`, upper identity block `[I_rank; 0]`
///   - `u`: `rank × n`, upper-echelon pivot rows
///   - `rem`: `m × n`, zero in the first `rank` rows
///
/// For `PivotType::Cols`:
///   - `l`: `m × rank`, lower-echelon pivot columns
///   - `u`: `rank × n`, left identity block `[I_rank | 0]`
///   - `rem`: `m × n`, zero in the first `rank` columns
pub struct PartialPluq<R> {
    pub p: PermOwned,
    pub q: PermOwned,
    pub l: SpMat<R>,
    pub u: SpMat<R>,
    pub rem: SpMat<R>,
}

impl<R> PartialPluq<R> {
    pub fn rank(&self) -> usize { self.l.ncols() }
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

    let (l, u, rem) = split(a, piv_type, &cur_p, &cur_q, total_r);

    PartialPluq { p: cur_p, q: cur_q, l, u, rem }
}

// Sparse identity block: `shape`-sized matrix with 1s at (k, k) for k < r.
fn eye<R>(shape: (usize, usize), r: usize) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    SpMat::from_entries(shape, (0..r).map(|k| (k, k, R::one())))
}

// Splits permuted `a` into `(l, u, rem)` satisfying `p * a * q = l * u + rem`.
//
// Rows: l = [I_r; 0] (m×r), u = pivot rows (r×n), rem zeros in rows 0..r.
// Cols: l = pivot cols (m×r), u = [I_r | 0] (r×n), rem zeros in cols 0..r.
fn split<R>(a: &SpMat<R>, piv_type: PivotType, p: &PermOwned, q: &PermOwned, r: usize) -> (SpMat<R>, SpMat<R>, SpMat<R>)
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (m, n) = a.shape();
    match piv_type {
        PivotType::Rows => {
            // u = pivot rows (r×n); rem = non-pivot entries shifted to rows r..m
            let mut u_ent = vec![];
            let mut rem_ent = vec![];
            for (i, j, v) in a.iter() {
                let pi = p.at(i);
                let qj = q.at(j);
                if pi < r {
                    u_ent.push((pi, qj, v.clone()));
                } else {
                    rem_ent.push((pi, qj, v.clone()));
                }
            }
            let l = eye((m, r), r);
            let u = SpMat::from_entries((r, n), u_ent);
            let rem = SpMat::from_entries((m, n), rem_ent);
            (l, u, rem)
        }
        PivotType::Cols => {
            // l = pivot cols (m×r); rem = non-pivot entries shifted to cols r..n
            let mut l_ent = vec![];
            let mut rem_ent = vec![];
            for (i, j, v) in a.iter() {
                let pi = p.at(i);
                let qj = q.at(j);
                if qj < r {
                    l_ent.push((pi, qj, v.clone()));
                } else {
                    rem_ent.push((pi, qj, v.clone()));
                }
            }
            let l = SpMat::from_entries((m, r), l_ent);
            let u = eye((r, n), r);
            let rem = SpMat::from_entries((m, n), rem_ent);
            (l, u, rem)
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
        let r = pp.rank();

        assert_eq!(pp.l.shape(), (m, r));
        assert_eq!(pp.u.shape(), (r, n));
        assert_eq!(pp.rem.shape(), (m, n));

        // p * a * q = l * u + rem
        let paq = a.permute(pp.p.view(), pp.q.view());
        assert_eq!(paq, &pp.l * &pp.u + &pp.rem);

        // u's first r×r block: pivots on diagonal, zeros below
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

    fn check_cols(a: &SpMat<i32>) {
        let pp = pre_pluq(a, PivotType::Cols, PivotCondition::AnyUnit, false);
        let (m, n) = a.shape();
        let r = pp.rank();

        assert_eq!(pp.l.shape(), (m, r));
        assert_eq!(pp.u.shape(), (r, n));
        assert_eq!(pp.rem.shape(), (m, n));

        // p * a * q = l * u + rem
        let paq = a.permute(pp.p.view(), pp.q.view());
        assert_eq!(paq, &pp.l * &pp.u + &pp.rem);

        // l's first r×r block: pivots on diagonal, zeros above
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

    // Recursive variants: verify shape, reconstruction, and diagonal pivots.
    fn check_rows_recursive(a: &SpMat<i32>) {
        let r0 = pre_pluq(a, PivotType::Rows, PivotCondition::AnyUnit, false).rank();
        let pp = pre_pluq(a, PivotType::Rows, PivotCondition::AnyUnit, true);
        let (m, n) = a.shape();
        let r = pp.rank();

        assert!(r >= r0);
        assert_eq!(pp.l.shape(), (m, r));
        assert_eq!(pp.u.shape(), (r, n));
        assert_eq!(pp.rem.shape(), (m, n));

        let paq = a.permute(pp.p.view(), pp.q.view());
        assert_eq!(paq, &pp.l * &pp.u + &pp.rem);

        let b = pp.u.clone().into_dense();
        for k in 0..r {
            assert!(b[(k, k)].is_one(), "u[{k},{k}] should be a pivot (=1)");
        }
    }

    fn check_cols_recursive(a: &SpMat<i32>) {
        let r0 = pre_pluq(a, PivotType::Cols, PivotCondition::AnyUnit, false).rank();
        let pp = pre_pluq(a, PivotType::Cols, PivotCondition::AnyUnit, true);
        let (m, n) = a.shape();
        let r = pp.rank();

        assert!(r >= r0);
        assert_eq!(pp.l.shape(), (m, r));
        assert_eq!(pp.u.shape(), (r, n));
        assert_eq!(pp.rem.shape(), (m, n));

        let paq = a.permute(pp.p.view(), pp.q.view());
        assert_eq!(paq, &pp.l * &pp.u + &pp.rem);

        let b = pp.l.clone().into_dense();
        for k in 0..r {
            assert!(b[(k, k)].is_one(), "l[{k},{k}] should be a pivot (=1)");
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
        assert_eq!(pp.rank(), 0);
        assert_eq!(pp.l.shape(), (4, 0));
        assert_eq!(pp.u.shape(), (0, 5));
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
        assert_eq!(pp.rank(), 3);
        assert_eq!(pp.l.shape(), (3, 3));
        assert_eq!(pp.u.shape(), (3, 3));
        assert_eq!(pp.rem.shape(), (3, 3));
        assert!(pp.rem.is_zero());
    }

    #[test]
    fn test_rank_rows() {
        assert_eq!(pre_pluq(&sample(), PivotType::Rows, PivotCondition::AnyUnit, false).rank(), 5);
    }

    #[test]
    fn test_rank_cols() {
        assert_eq!(pre_pluq(&sample(), PivotType::Cols, PivotCondition::AnyUnit, false).rank(), 6);
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