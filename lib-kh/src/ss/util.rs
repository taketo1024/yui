//! The `H = 1` divisibility machinery shared by `ss` and `ssi`.
//!
//! Over `F[H]` the `H`-divisibility `d_H` of a Lee class is read off from a family of linear
//! systems `d(w) = z` over `F`, posed in q-truncations of the Bar-Natan specialization. Level `m`
//! (window `q < q₀ + 2m`) is solvable exactly when `d_H ≥ m`, so the first failure gives `d_H`.
//! No homology, SNF or basis tracking is involved.
//!
//! Reference: `research/h1-divisibility.tex`.

use itertools::Itertools;
use log::debug;

use yui_core::lc::Lc;
use yui_core::abst::{EucRing, EucRingOps, Field, FieldOps};
use yui_core::lc::LcKey;
use yui_matrix::MatTrait;
use yui_matrix::sparse::{SpMat, SpVec};
use yui_matrix::sparse::pluq::solve_pluq;

use crate::util::FastPoly;

type P<F> = FastPoly<'H', F>;

/// Which pipeline computes the divisibility.
/// - `V1`: from the bigraded homology over `F[H]`; simple, memory-heavy.
/// - `V2`: the `H = 1` specialization below — q-truncated `F` solves, no homology.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum SsVersion {
    V1,
    #[default]
    V2,
}

/// The canon cycle's homogeneous q-degree `q₀ = w − r` (`+1` reduced), where `w` is the writhe
/// and `r` the number of Seifert circles.
pub(super) fn canon_q_deg(w: i32, r: usize, reduced: bool) -> isize {
    w as isize - r as isize + if reduced { 1 } else { 0 }
}

/// Climb `m` while the level-`m` system `d(w) = z` in `C̄ / F_{q₀+2m}` solves; the first failure
/// gives `d_H = m`. `None` means the built window is too tight to trust — either the canon
/// representative was truncated away, or the climb hit the built ceiling still solvable.
pub(super) fn max_solvable_level<F>(
    d: &SpMat<P<F>>, v: &SpVec<P<F>>, src_q: &[isize], tgt_q: &[isize], q0: isize, at: isize,
) -> Option<i32>
where F: Field, for<'x> &'x F: FieldOps<F> {
    // beyond this the window contains everything the build kept.
    let q_top = Iterator::chain(src_q.iter(), tgt_q.iter()).max().copied().unwrap_or(q0);

    let mut m = 0;
    loop {
        let q_hi = q0 + 2 * (m + 1);
        if q_hi > q_top + 2 {
            return None; // reached the built ceiling still solvable — widen for the true failure
        }

        let (a, y) = truncated_system(d, v, src_q, tgt_q, q_hi);
        debug!("h = {at}, level {}: solve system of size {:?}.", m + 1, a.shape());

        if solve_pluq(&a, &y).is_some() {
            m += 1;
        } else {
            debug!("h = {at}: d = {m}.");
            return Some(m as i32);
        }
    }
}

// The level system: rows/cols of `d` restricted to the q-window `q < q_hi` (the full quotient
// `C̄ / F_{q_hi}`), entries evaluated at `H = 1`.
fn truncated_system<F>(
    d: &SpMat<P<F>>, v: &SpVec<P<F>>,
    src_q: &[isize], tgt_q: &[isize], q_hi: isize,
) -> (SpMat<F>, SpVec<F>)
where F: Field, for<'x> &'x F: FieldOps<F> {
    let in_win = |q: isize| q < q_hi;
    let (n_rows, row_map) = window_index(tgt_q, &in_win);
    let (n_cols, col_map) = window_index(src_q, &in_win);

    let a = SpMat::from_entries((n_rows, n_cols), d.iter_nz().filter_map(|(i, j, p)| {
        let i1 = row_map[i]?;
        let j1 = col_map[j]?;
        debug_assert!(
            2 * (p.deg() as isize) == tgt_q[i] - src_q[j],
            "the differential is not q-homogeneous"
        );
        let e = eval_at_one(p);
        (!e.is_zero()).then_some((i1, j1, e))
    }));

    let y = SpVec::from_entries(n_rows, v.iter_nz().filter_map(|(i, p)| {
        let i1 = row_map[i]?;
        let e = eval_at_one(p);
        (!e.is_zero()).then_some((i1, e))
    }));

    (a, y)
}

fn window_index(qs: &[isize], in_win: impl Fn(isize) -> bool) -> (usize, Vec<Option<usize>>) {
    let mut n = 0;
    let map = qs.iter().map(|&q|
        in_win(q).then(|| {
            let i = n;
            n += 1;
            i
        })
    ).collect();
    (n, map)
}

fn eval_at_one<F>(p: &P<F>) -> F
where F: Field, for<'x> &'x F: FieldOps<F> {
    p.coeff().clone()
}

/// Every coefficient of `z` is a monomial `λ·H^k` (by the `FastPoly` representation); check that
/// `q(x) − 2k = q₀` for every term.
pub(super) fn assert_homogeneous<X, F>(z: &Lc<X, P<F>>, q_deg: impl Fn(&X) -> isize, q0: isize)
where X: LcKey, F: Field, for<'x> &'x F: FieldOps<F> {
    let bad = z.iter().filter(|(x, a)|
        q_deg(x) - 2 * (a.deg() as isize) != q0
    ).map(|(x, a)|
        format!("  q(x) = {}, coeff = {} (deg {})", q_deg(x), a, a.deg())
    ).collect_vec();

    assert!(
        bad.is_empty(),
        "canon cycle is not q-homogeneous at q0 = {q0}:\n{}", bad.iter().join("\n")
    );
}

/// The largest `k` with `c^k` dividing every entry of `v`; `None` if `v` is zero.
pub fn div_vec<R>(v: &SpVec<R>, c: &R) -> Option<i32>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    v.iter().filter_map(|(_, a)|
        div(a, c)
    ).min()
}

fn div<R>(a: &R, c: &R) -> Option<i32>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    if a.is_zero() { return None }

    let mut a = a.clone();
    let mut k = 0;

    while (&a % c).is_zero() {
        a /= c;
        k += 1;
    }

    Some(k)
}
