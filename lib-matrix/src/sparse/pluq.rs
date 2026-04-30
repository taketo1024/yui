// Sparse PLUQ decomposition & linear solver.
// Implemented with the help of Claude Code.

use log::debug;
use sprs::PermOwned;
use sprs::PermView;
use yui_core::{Ring, RingOps, Field, FieldOps};

use crate::MatTrait;
use crate::dense::pluq::solve_pluq as dense_solve_pluq;
use super::SpMat;
use super::SpVec;
use super::pivot::{PivotCondition, PivotFinderConfig, PivotType, find_pivots, perms_by_pivots};
use super::triang::{TriangularType, solve_triangular, solve_triangular_left, solve_triangular_vec};
use super::util::perm_for_indices;

/// Result of a partial PLUQ decomposition.
///
/// Satisfies `p * A * q = l * u + s` where `s` is the
/// `(m - rank) × (n - rank)` Schur complement (bottom-right block).
pub struct PartialPluq<R> {
    pub p: PermOwned,
    pub q: PermOwned,
    pub l: SpMat<R>,
    pub u: SpMat<R>,
    pub s: SpMat<R>,
}

impl<R> PartialPluq<R> {
    pub fn rank(&self) -> usize { self.l.ncols() }
}

/// Computes a partial PLUQ decomposition of `a`.
pub fn pre_pluq<R>(a: &SpMat<R>, piv_type: PivotType, piv_cond: PivotCondition) -> PartialPluq<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    debug!("pre PLUQ: {:?}", a.shape());

    let config = PivotFinderConfig { piv_type, piv_cond, ..Default::default() };
    let pivots = find_pivots(a, config);
    let r = pivots.len();
    let (p, q) = perms_by_pivots(a, &pivots);
    let paq = split(a, &p, &q, r);
    let (l, u, s) = build(piv_type, paq);

    PartialPluq { p, q, l, u, s }
}

// Applies permutations (p, q) to `a` and partitions the result into four blocks at row/col r:
//
//   paq = [[a0 | a1],   a0: r×r,     a1: r×(n-r)
//          [a2 | a3]]   a2: (m-r)×r, a3: (m-r)×(n-r)
fn split<R>(a: &SpMat<R>, p: &PermOwned, q: &PermOwned, r: usize) -> [SpMat<R>; 4]
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
fn build<R>(piv_type: PivotType, paq: [SpMat<R>; 4]) -> (SpMat<R>, SpMat<R>, SpMat<R>)
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

/// Solves `a * x = y` over a field using sparse PLUQ.
///
/// Returns `Some(x)` if a solution exists, `None` otherwise.
///
/// Uses `pre_pluq` with column pivots. After the Schur extension:
///   `p * a * q = l * u + s`,  u = [I_r | U1].
/// Algorithm:
///   z      = L0^{-1} * yp[0..r]  (forward substitution)
///   yp_res = yp[r..m] - L1 * z
///   xq_res: solve S * xq_res = yp_res  (S = pp.s)
///   xq_top  = z - U1 * xq_res          (U1 = u[0..r, r..n])
pub fn solve_pluq<R>(a: &SpMat<R>, y: &[R]) -> Option<Vec<R>>
where R: Field, for<'x> &'x R: FieldOps<R> {
    let (m, n) = a.shape();
    assert_eq!(y.len(), m);

    debug!("sparse solve: {:?}", a.shape());

    let pp = pre_pluq(a, PivotType::Cols, PivotCondition::AnyUnit);
    let r = pp.rank();
    let yp = perm_apply(pp.p.view(), y);

    debug!("solve top: {:?}", pp.l.shape());

    let z = solve_top(&pp.l, &yp);
    let yp_res = compute_yp_res(&pp.l, &yp, &z);

    debug!("solve res: {:?}", pp.s.shape());
    
    let xq_res = solve_res(&pp.s, &yp_res)?;
    let u1 = pp.u.submat(0..r, r..n);
    let xq_top = back_sub_piv(&z, &u1, &xq_res);

    Some(reconstruct_x(&pp.q, r, &xq_top, &xq_res))
}

// Applies permutation p to y: yp[p(i)] = y[i].
fn perm_apply<R>(p: PermView, y: &[R]) -> Vec<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let pinv = p.inv();
    (0..y.len()).map(|i| y[pinv.at(i)].clone()).collect()
}

// Solves l[0..r, 0..r] * x = yp[0..r] by forward substitution.
// Requires l[0..r, 0..r] to be lower triangular with non-zero diagonal.
fn solve_top<R>(l: &SpMat<R>, yp: &[R]) -> Vec<R>
where R: Field, for<'x> &'x R: FieldOps<R> {
    let r = l.ncols();
    if r == 0 { return vec![]; }

    let l0 = l.submat(0..r, 0..r);
    let b = SpVec::from(yp[..r].to_vec());
    solve_triangular_vec(TriangularType::Lower, &l0, &b).to_dense()
}

// Computes yp[r..m] - l[r..m, :] * x_piv where r = x_piv.len().
fn compute_yp_res<R>(l: &SpMat<R>, yp: &[R], x_piv: &[R]) -> Vec<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let r = x_piv.len();
    let mut yp_res = yp[r..].to_vec();
    for (i, j, v) in l.iter_nz() {
        if i >= r {
            yp_res[i - r] = yp_res[i - r].clone() - v * &x_piv[j];
        }
    }
    yp_res
}

// Solves s * xq_res = yp_res.
// Returns None if any zero-row of s has nonzero yp_res (inconsistency).
// Compresses to a dense subsystem over non-zero rows/cols and delegates to dense_solve_pluq.
fn solve_res<R>(s: &SpMat<R>, yp_res: &[R]) -> Option<Vec<R>>
where R: Field, for<'x> &'x R: FieldOps<R> {
    use std::collections::BTreeSet;

    let (m, n) = s.shape();
    assert_eq!(yp_res.len(), m);

    let mut row_set: BTreeSet<usize> = BTreeSet::new();
    let mut col_set: BTreeSet<usize> = BTreeSet::new();
    for (i, j, _) in s.iter_nz() { 
        row_set.insert(i); 
        col_set.insert(j); 
    }

    let m0 = row_set.len();
    let n0 = col_set.len();
    let row_perm = perm_for_indices(m, row_set.iter());
    let col_perm = perm_for_indices(n, col_set.iter());

    // Consistency check: zero rows must have zero yp_res.
    if (0..m).any(|i| row_perm.at(i) >= m0 && !yp_res[i].is_zero()) {
        return None;
    }

    if m0 == 0 {
        return Some(vec![R::zero(); n]);
    }

    let s0 = s.extract((m0, n0), |i, j| {
        let (ri, ci) = (row_perm.at(i), col_perm.at(j));
        (ri < m0 && ci < n0).then_some((ri, ci))
    }).into_dense();

    let rhs = perm_apply(row_perm.view(), yp_res);
    let mut xq = dense_solve_pluq(&s0, &rhs[..m0])?;
    xq.resize(n, R::zero());

    Some(perm_apply(col_perm.inv(), &xq))
}

// Computes xq_top = z - U1 * xq_res.
// After Cols extension, u = [I_r | U1], so U1 = u[0..r, r..n].
fn back_sub_piv<R>(z: &[R], u1: &SpMat<R>, xq_res: &[R]) -> Vec<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let mut xq_top = z.to_vec();
    for (i, j, v) in u1.iter_nz() {
        xq_top[i] = xq_top[i].clone() - v * &xq_res[j];
    }
    xq_top
}

// Reconstructs x from permuted solution: x[j] = x'[q(j)], x' = [x_piv | x_free].
fn reconstruct_x<R: Clone>(q: &PermOwned, r: usize, x_piv: &[R], x_free: &[R]) -> Vec<R> {
    let n = r + x_free.len();
    (0..n).map(|j| {
        let qj = q.at(j);
        if qj < r { x_piv[qj].clone() } else { x_free[qj - r].clone() }
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::One;

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
        let pp = pre_pluq(a, PivotType::Rows, PivotCondition::AnyUnit);
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

    fn check_cols(a: &SpMat<i32>) {
        let pp = pre_pluq(a, PivotType::Cols, PivotCondition::AnyUnit);
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
    fn test_rows() { check_rows(&sample()); }

    #[test]
    fn test_cols() { check_cols(&sample()); }

    #[test]
    fn test_zero() {
        let a = SpMat::<i32>::zero((4, 5));
        let pp = pre_pluq(&a, PivotType::Rows, PivotCondition::AnyUnit);
        assert_eq!(pp.rank(), 0);
        assert_eq!(pp.l.shape(), (4, 0));
        assert_eq!(pp.u.shape(), (0, 5));
        assert_eq!(pp.s.shape(), (4, 5)); // (m-r, n-r) = (4, 5) when r=0
        assert_eq!(pp.s, a.permute(pp.p.view(), pp.q.view()));
    }

    #[test]
    fn test_square_full_rank() {
        let a = SpMat::from_dense_data((3, 3), [1, 0, 0, 0, 1, 0, 0, 0, 1]);
        let pp = pre_pluq(&a, PivotType::Rows, PivotCondition::AnyUnit);
        assert_eq!(pp.rank(), 3);
        assert_eq!(pp.s.shape(), (0, 0)); // full rank: Schur complement is empty
    }

    #[test]
    fn test_rank_rows() {
        assert_eq!(pre_pluq(&sample(), PivotType::Rows, PivotCondition::AnyUnit).rank(), 5);
    }

    #[test]
    fn test_rank_cols() {
        assert_eq!(pre_pluq(&sample(), PivotType::Cols, PivotCondition::AnyUnit).rank(), 6);
    }

    #[test]
    fn test_rand_rows() { check_rows(&SpMat::<i32>::rand((40, 60), 0.1)); }

    #[test]
    fn test_rand_cols() { check_cols(&SpMat::<i32>::rand((40, 60), 0.1)); }

    // ---- helper unit tests ----

    #[test]
    fn test_split() {
        use sprs::PermOwned;
        // a = [[1,2],[3,4]], r=1, identity perms → paq = a, partition at row/col 1:
        // a0=[[1]], a1=[[2]], a2=[[3]], a3=[[4]]
        let a = sp((2, 2), [r(1), r(2), r(3), r(4)]);
        let p = PermOwned::new(vec![0, 1]);
        let q = PermOwned::new(vec![0, 1]);
        let [a0, a1, a2, a3] = split(&a, &p, &q, 1);
        assert_eq!(a0, sp((1, 1), [r(1)]));
        assert_eq!(a1, sp((1, 1), [r(2)]));
        assert_eq!(a2, sp((1, 1), [r(3)]));
        assert_eq!(a3, sp((1, 1), [r(4)]));
    }

    #[test]
    fn test_build_cols() {
        // paq = [[1,2],[3,4]], r=1: a0=[[1]], a1=[[2]], a2=[[3]], a3=[[4]]
        // l0=a0=[[1]], l1=a2=[[3]], r0=a1=[[2]], r1=a3=[[4]]
        // u1 = [[1]]^{-1}*[[2]] = [[2]], l = [[1],[3]], u = [[1,2]], s = [[4]]-[[3]]*[[2]] = [[-2]]
        let paq = [sp((1,1),[r(1)]), sp((1,1),[r(2)]), sp((1,1),[r(3)]), sp((1,1),[r(4)])];
        let (l, u, s) = build(PivotType::Cols, paq);
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
        let (l, u, s) = build(PivotType::Rows, paq);
        assert_eq!(l, sp((2, 1), [r(1), r(3)]));
        assert_eq!(u, sp((1, 2), [r(1), r(2)]));
        assert_eq!(s, sp((1, 1), [r(-2)]));
    }

    use yui_core::num::Ratio;
    type R = Ratio<i64>;
    fn r(n: i64) -> R { R::from(n) }

    fn sp(shape: (usize, usize), data: impl IntoIterator<Item = R>) -> SpMat<R> {
        SpMat::from_dense_data(shape, data)
    }

    #[test]
    fn test_perm_apply() {
        use sprs::PermOwned;
        let p = PermOwned::new(vec![1, 2, 0]); // 0→1, 1→2, 2→0
        let y = vec![r(10), r(20), r(30)];
        let yp = perm_apply(p.view(), &y);
        // yp[p(0)=1]=10, yp[p(1)=2]=20, yp[p(2)=0]=30
        assert_eq!(yp, vec![r(30), r(10), r(20)]);
    }

    #[test]
    fn test_solve_top() {
        // l = [[2, 0], [3, 4]], yp = [4, 11]
        // l[0..2,0..2]*x = [4,11] → x = [2, 5/4]
        let l = sp((2, 2), [r(2), r(0), r(3), r(4)]);
        let yp = vec![r(4), r(11)];
        let x = solve_top(&l, &yp);
        assert_eq!(x, vec![r(2), r(5)/r(4)]);
    }

    #[test]
    fn test_compute_yp_res() {
        // l = [[1,0],[2,1],[3,0]], x_piv = [3,2], yp = [*,*,7,8]
        // yp_res = yp[2..] - l[2..,:]*x_piv = [7-9, 8-0] = [-2, 8]
        let l = sp((3, 2), [r(1), r(0), r(2), r(1), r(3), r(0)]);
        let yp = vec![r(0), r(0), r(7), r(8)]; // only [2..] matters
        let x_piv = vec![r(3), r(2)];
        let res = compute_yp_res(&l, &yp, &x_piv);
        assert_eq!(res, vec![r(7) - r(9), r(8)]);
    }

    #[test]
    fn test_solve_res_trivial() {
        // all-zero s, zero yp_res → x_free = [0, 0]
        let rem = SpMat::<R>::zero((2, 2));
        let yp_res = vec![r(0), r(0)];
        assert_eq!(solve_res(&rem, &yp_res), Some(vec![r(0), r(0)]));
    }

    #[test]
    fn test_solve_res_inconsistent() {
        // zero rem, nonzero yp_res → None
        let rem = SpMat::<R>::zero((2, 2));
        let yp_res = vec![r(1), r(0)];
        assert!(solve_res(&rem, &yp_res).is_none());
    }

    #[test]
    fn test_solve_res_nontrivial() {
        // rem = [[2, 0], [0, 3]], yp_res = [4, 6] → x_free = [2, 2]
        let rem = sp((2, 2), [r(2), r(0), r(0), r(3)]);
        let yp_res = vec![r(4), r(6)];
        let x = solve_res(&rem, &yp_res).unwrap();
        // check rem * x = yp_res
        let ax0 = r(2) * x[0].clone();
        let ax1 = r(3) * x[1].clone();
        assert_eq!(ax0, r(4));
        assert_eq!(ax1, r(6));
    }

    #[test]
    fn test_reconstruct_x() {
        use sprs::PermOwned;
        // q = [1, 0] (swap), r=1, x_piv=[10], x_free=[20]
        // j=0: q(0)=1 >= r=1 → x_free[0]=20
        // j=1: q(1)=0 < r=1  → x_piv[0]=10
        let q = PermOwned::new(vec![1, 0]);
        let x = reconstruct_x(&q, 1, &[r(10)], &[r(20)]);
        assert_eq!(x, vec![r(20), r(10)]);
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
}
