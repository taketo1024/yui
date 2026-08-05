//! `ssi` via specialization at `H = 1`: the divisibility `d_H` of each equivariant Lee class
//! is read off from a family of linear systems `d(w) = z` over `𝔽₂`, posed in q-truncations of
//! the Bar-Natan specialization of `CKhI(D)`. No homology, SNF or basis tracking is involved.
//!
//! Reference: `research/h1-divisibility.tex`. Solvability at level `m` (window `q < q₀ + 2m`)
//! is equivalent to `d_H([z]) ≥ m`, the `H`-divisibility modulo torsion — exactly the quantity
//! defining `ssi`.

use itertools::Itertools;
use num_traits::Zero;
use log::info;

use yui_core::num::FF2;
use yui_matrix::MatTrait;
use yui_matrix::sparse::{SpMat, SpVec};
use yui_matrix::sparse::pluq::solve_pluq;
use yui_link::{InvLink, Link};

use crate::kh::KhComplex;
use crate::tng::builder::SymBuildConfig;
use crate::util::FastPoly;
use crate::khi::{KhIChain, KhIComplex};

type F = FF2;
type P = FastPoly<'H', F>;

// The `H`-divisibilities `(d0, d1)` via `H = 1` solves. `expected` (the guessed s-value) sets the
// initial high-q build cut at `q = expected − 1` — only the high end is cut, so `d_H` (mod torsion)
// is unchanged. A too-tight window (the reduction lifts the canon representative above it) is
// widened by 2 and retried. `None` = full (un-windowed) build.
pub(crate) fn ssi_divisibility_v2(l: &InvLink, reduced: bool, config: SymBuildConfig, expected: Option<isize>) -> (i32, i32) {
    let q0 = canon_q_deg(l, reduced);

    let Some(s) = expected else {
        return divisibility_in_window(l, reduced, config, q0)
            .expect("the full (un-windowed) build must not truncate the canon cycle");
    };

    let q_hi0 = s - 1; // the top divisibility generator for the guessed s sits at q = s − 1.
    let mut q_hi = q_hi0;
    loop {
        info!("q-window: q0 = {q0}, cut above {q_hi}.");
        let cfg = SymBuildConfig { q_range: Some((isize::MIN + 1) ..= q_hi), ..config.clone() };

        if let Some(ds) = divisibility_in_window(l, reduced, cfg, q0) {
            return ds;
        }
        q_hi += 2; // one divisibility level = q-degree 2 (deg H = −2); raise the cut minimally.
        info!("canon cycle above window; raising cut to {q_hi}.");
        assert!(q_hi <= q_hi0 + 128, "q-window widening runaway");
    }
}

// The canon cycle's homogeneous q-degree `q₀ = w − r` (`+1` reduced), where `w` = writhe and
// `r` = number of Seifert circles.
fn canon_q_deg(l: &InvLink, reduced: bool) -> isize {
    let w = l.writhe() as isize;
    let r = l.seifert_circles().len() as isize;
    w - r + if reduced { 1 } else { 0 }
}

// `None` signals the q-window was too tight (a canon cycle's low-q representative was truncated) —
// the caller should widen and retry.
fn divisibility_in_window(l: &InvLink, reduced: bool, config: SymBuildConfig, q0: isize) -> Option<(i32, i32)> {
    let r = if reduced { 1 } else { 2 };
    let h = P::variable();
    let t = P::zero();

    // Build over `(a−1)..=1`: the solves only involve `C[−1] → C[0] → C[1]`, so nothing above 1
    // is built; the cheap low degrees stay (truncating the bottom makes the build frontier dense).
    // Only KhI degrees `−1..=1` are converted to the raw complex — the rest is never expanded.
    let requested = config.h_range.clone().unwrap_or(-(Link::MAX_CROSSING as isize) ..= 1);
    let range = KhComplex::<P>::clamp_h_range(l.inner(), reduced, requested);
    let (a, b) = (*range.start(), *range.end());
    assert!(a <= 0 && b >= 1, "ssi h-range must include 0 and 1, got {a}..={b}");
    let config = SymBuildConfig { h_range: Some((a - 1)..=b), ..config };
    let kc = KhIComplex::<P>::new_windowed(l, &h, &t, reduced, config, -1..=1);

    let zs = kc.canon_cycles(); // sorted by h-degree: B classes at 0, then Q classes at 1
    assert_eq!(zs.len(), 2 * r);

    // index gives the h-degree (0 for the first `r`, else 1) — robust to a cycle truncated to zero.
    let mut ds = [vec![], vec![]];
    for (i, z) in zs.iter().enumerate() {
        let i0 = if i < r { 0 } else { 1 };
        ds[i0 as usize].push(max_solvable_level(&kc, z, i0, q0)?);
    }

    let (d0, d1) = if reduced {
        (ds[0][0], ds[1][0])
    } else {
        assert_eq!(ds[0][0], ds[0][1]);
        assert_eq!(ds[1][0], ds[1][1]);
        (ds[0][0], ds[1][0])
    };

    assert!(d0 <= d1); // s̲ ≤ s̄.

    Some((d0, d1))
}

// Climb `m` while the level-`m` system `d(w) = z` in `C̄ / F_{q₀+2m}` solves; the first failure gives
// `d_H = m`. Returns `None` when the built window is too tight to trust: the canon representative's
// lowest generator sits above `q₀` (truncated), or the climb reaches the built ceiling still solvable.
fn max_solvable_level(kc: &KhIComplex<P>, z: &KhIChain<P>, i0: isize, q0: isize) -> Option<i32> {
    // an empty cycle means the window cut the whole representative (it sits entirely above `q_hi`) —
    // signal widen. A nonempty cycle is homogeneous at `q₀` (every term has `q(x) − 2·deg_H = q₀`),
    // so the truncated low part is enough to run the solve.
    if z.is_zero() {
        return None;
    }
    assert_homogeneous(kc, z, q0);

    let src_q = q_degs(kc, i0 - 1);
    let tgt_q = q_degs(kc, i0);
    let d = kc.inner().d_matrix(i0 - 1);
    let v = kc.inner()[i0].vectorize(z);

    // beyond this the window contains everything the build kept.
    let q_top = Iterator::chain(src_q.iter(), tgt_q.iter()).max().copied().unwrap_or(q0);

    let mut m = 0;
    loop {
        let q_hi = q0 + 2 * (m + 1);
        if q_hi > q_top + 2 {
            return None; // reached the built ceiling still solvable — widen for the true failure
        }

        let (a, y) = truncated_system(&d, &v, &src_q, &tgt_q, q_hi);
        info!("h = {i0}, level {}: solve system of size {:?}.", m + 1, a.shape());

        if solve_pluq(&a, &y).is_some() {
            m += 1;
        } else {
            info!("h = {i0}: d = {m}.");
            return Some(m as i32);
        }
    }
}

// The level system: rows/cols of `d^{i0-1}` restricted to the q-window `q < q_hi`
// (the full quotient `C̄ / F_{q_hi}`), entries evaluated at `H = 1`.
fn truncated_system(
    d: &SpMat<P>, v: &SpVec<P>,
    src_q: &[isize], tgt_q: &[isize], q_hi: isize,
) -> (SpMat<F>, SpVec<F>) {
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

fn q_degs(kc: &KhIComplex<P>, i: isize) -> Vec<isize> {
    kc.inner()[i].raw_generators().iter().map(|x| kc.q_deg_of(x)).collect()
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

fn eval_at_one(p: &P) -> F {
    *p.coeff()
}

// every coefficient of `z` is a monomial `λ·H^k` (by the `FastPoly` representation); check `q(x) − 2k = q₀`.
fn assert_homogeneous(kc: &KhIComplex<P>, z: &KhIChain<P>, q0: isize) {
    let bad = z.iter().filter(|(x, a)|
        kc.q_deg_of(*x) - 2 * (a.deg() as isize) != q0
    ).map(|(x, a)|
        format!("  q(x) = {}, coeff = {} (deg {})", kc.q_deg_of(x), a, a.deg())
    ).collect_vec();

    assert!(
        bad.is_empty(),
        "canon cycle is not q-homogeneous at q0 = {q0}:\n{}", bad.iter().join("\n")
    );
}
