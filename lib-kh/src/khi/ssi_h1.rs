//! `ssi` via specialization at `H = 1`: the divisibility `d_H` of each equivariant Lee class
//! is read off from a family of linear systems `d(w) = z` over `F`, posed in q-truncations of
//! the Bar-Natan specialization of `CKhI(D)`. No homology, SNF or basis tracking is involved.
//!
//! Reference: `research/h1-divisibility.tex`. Solvability at level `m` (window `q < q₀ + 2m`)
//! is equivalent to `d_H([z]) ≥ m`, the `H`-divisibility modulo torsion — exactly the quantity
//! defining `ssi`.

use itertools::Itertools;
use num_traits::Zero;
use log::info;

use yui_core::{Field, FieldOps, MathType};
use yui_core::poly::Poly;
use yui_matrix::MatTrait;
use yui_matrix::sparse::{SpMat, SpVec};
use yui_matrix::sparse::pluq::solve_pluq;
use yui_link::{InvLink, Link};

use crate::kh::KhComplex;
use crate::tng::builder::SymBuildConfig;
use crate::khi::{KhIChain, KhIComplex};

type PolyH<F> = Poly<'H', F>;

/// `ssi` from the `H = 1` linear systems. The complex is built over `F[H]` (q-homogeneous)
/// exactly as in `ssi_invariant`; only the solves happen over `F`. Char-2 only (the cone).
pub fn ssi_invariant_h1<F>(l: &InvLink, reduced: bool, config: SymBuildConfig) -> (i32, i32)
where F: Field, for<'x> &'x F: FieldOps<F> {
    assert!(l.is_knot());
    info!("compute ssi via H=1 solves over {}.", PolyH::<F>::math_symbol());

    let q0 = canon_q_deg(l, reduced);
    let (d0, d1) = ssi_divisibility_h1::<F>(l, reduced, config, q0)
        .expect("full (un-windowed) build must not truncate the canon cycle");
    finish_ssi(l, d0, d1)
}

/// `ssi_invariant_h1` with the build's high-q cut set from a target divisibility level: only
/// generators of q-degree `≤ q₀ + 2·(level + 2)` are built (`q₀ = w − r`). Only the high end is cut,
/// so `d_H` (mod torsion) is unchanged, and the solve still checks every level from `q₀`. If the
/// window is too tight — the reduction lifts the canon representative above it (truncated to zero) —
/// the build is retried one level wider.
///
/// `start_from` is the initial level. `None` uses the `ssi ≈ 0` heuristic `−(q₀+1)/2 − 1`: slice
/// knots have `ssi` near 0, hence `d_H ≈ −(q₀+1)/2` (large, since `q₀` is very negative), so this
/// starts near the answer and skips the doomed small windows. This is the memory/time lever for
/// large diagrams (Wh-doubles).
pub fn ssi_invariant_h1_windowed<F>(l: &InvLink, reduced: bool, config: SymBuildConfig, start_from: Option<i32>) -> (i32, i32)
where F: Field, for<'x> &'x F: FieldOps<F> {
    assert!(l.is_knot());
    info!("compute ssi via windowed H=1 solves over {}.", PolyH::<F>::math_symbol());

    let q0 = canon_q_deg(l, reduced);
    let md0 = start_from.unwrap_or(((-(q0 + 1) / 2) as i32 - 1).max(0));
    let mut md = md0;
    loop {
        let q_hi = q0 + 2 * (md as isize + 2);
        info!("q-window: q0 = {q0}, cut above {q_hi} (level = {md}).");
        let cfg = SymBuildConfig { q_range: Some((isize::MIN + 1) ..= q_hi), ..config.clone() };

        if let Some((d0, d1)) = ssi_divisibility_h1::<F>(l, reduced, cfg, q0) {
            return finish_ssi(l, d0, d1);
        }
        md += 1; // one divisibility level = q-degree 2 (deg H = −2); widen by the minimal step.
        info!("canon cycle above window; widening to level = {md}.");
        assert!(md <= md0 + 64, "q-window widening runaway");
    }
}

fn finish_ssi(l: &InvLink, d0: i32, d1: i32) -> (i32, i32) {
    let w = l.writhe();
    let r = l.seifert_circles().len() as i32;
    let (ss0, ss1) = (2 * d0 + w - r + 1, 2 * d1 + w - r + 1);
    info!("w = {w}, r = {r}, d0 = {d0}, d1 = {d1}; ssi = ({ss0}, {ss1}).");
    (ss0, ss1)
}

// The canon cycle's homogeneous q-degree `q₀`, computed from the diagram: `w − r` (`+1` reduced),
// where `w` = writhe and `r` = number of Seifert circles.
fn canon_q_deg(l: &InvLink, reduced: bool) -> isize {
    let w = l.writhe() as isize;
    let r = l.seifert_circles().len() as isize;
    w - r + if reduced { 1 } else { 0 }
}

// `None` signals the q-window was too tight (a canon cycle's low-q representative was truncated) —
// the caller should widen and retry. `q0 = w − r` is the diagram-computed homogeneous degree.
fn ssi_divisibility_h1<F>(l: &InvLink, reduced: bool, config: SymBuildConfig, q0: isize) -> Option<(i32, i32)>
where F: Field, for<'x> &'x F: FieldOps<F> {
    let r = if reduced { 1 } else { 2 };
    let h = PolyH::<F>::variable();
    let t = PolyH::<F>::zero();

    // same windowing as `ssi_divisibility`: cover the canon degrees 0, 1, one wider on both ends.
    let requested = config.h_range.clone().unwrap_or(-(Link::MAX_CROSSING as isize) ..= 1);
    let range = KhComplex::<PolyH<F>>::clamp_h_range(l.inner(), reduced, requested);
    let (a, b) = (*range.start(), *range.end());
    assert!(a <= 0 && b >= 1, "ssi h-range must include 0 and 1, got {a}..={b}");
    let config = SymBuildConfig { h_range: Some((a - 1)..=(b + 1)), ..config };
    let kc = KhIComplex::<PolyH<F>>::new_with_config(l, &h, &t, reduced, config);

    let zs = kc.canon_cycles(); // sorted by h-degree: B classes at 0, then Q classes at 1
    assert_eq!(zs.len(), 2 * r);

    // index gives the h-degree (0 for the first `r`, else 1) — robust to a cycle truncated to zero.
    let mut ds = [vec![], vec![]];
    for (i, z) in zs.iter().enumerate() {
        let i0 = if i < r { 0 } else { 1 };
        ds[i0 as usize].push(max_solvable_level::<F>(&kc, z, i0, q0)?);
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
fn max_solvable_level<F>(kc: &KhIComplex<PolyH<F>>, z: &KhIChain<PolyH<F>>, i0: isize, q0: isize) -> Option<i32>
where F: Field, for<'x> &'x F: FieldOps<F> {
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
fn truncated_system<F>(
    d: &SpMat<PolyH<F>>, v: &SpVec<PolyH<F>>,
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
            p.nterms() == 1 && 2 * (p.lead_deg() as isize) == tgt_q[i] - src_q[j],
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

fn q_degs<F>(kc: &KhIComplex<PolyH<F>>, i: isize) -> Vec<isize>
where F: Field, for<'x> &'x F: FieldOps<F> {
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

fn eval_at_one<F>(p: &PolyH<F>) -> F
where F: Field, for<'x> &'x F: FieldOps<F> {
    F::sum(p.iter().map(|(_, a)| a.clone()))
}

// every coefficient of `z` must be a monomial `λ·H^k` with `q(x) − 2k = q₀`.
fn assert_homogeneous<F>(kc: &KhIComplex<PolyH<F>>, z: &KhIChain<PolyH<F>>, q0: isize)
where F: Field, for<'x> &'x F: FieldOps<F> {
    let bad = z.iter().filter(|(x, a)|
        a.nterms() != 1 || kc.q_deg_of(*x) - 2 * (a.lead_deg() as isize) != q0
    ).map(|(x, a)|
        format!("  q(x) = {}, coeff = {} (nterms {}, lead_deg {})", kc.q_deg_of(x), a, a.nterms(), a.lead_deg())
    ).collect_vec();

    assert!(
        bad.is_empty(),
        "canon cycle is not q-homogeneous at q0 = {q0}:\n{}", bad.iter().join("\n")
    );
}

#[cfg(test)]
mod tests {
    use yui_core::num::FF2;
    use super::*;

    type F = FF2;

    fn ssi_h1(l: &InvLink, reduced: bool) -> (i32, i32) {
        ssi_invariant_h1::<F>(l, reduced, SymBuildConfig::default())
    }

    macro_rules! test {
        ($test:ident, $name:literal, $expected:expr) => {
            #[test]
            fn $test() -> Result<(), Box<dyn std::error::Error>> {
                let l = InvLink::load($name)?;
                assert_eq!(ssi_h1(&l, false), $expected);
                Ok(())
            }
        }
    }

    test!(k3_1, "3_1", (2, 2));
    test!(k4_1, "4_1", (0, 0));
    test!(k5_1, "5_1", (4, 4));
    test!(k6_2a, "6_2a", (2, 2));
    test!(k6_3, "6_3", (0, 0));
    test!(k7_6a, "7_6a", (-2, -2));

    // the q-windowed build (high-q cut at the true divisibility) reproduces the full ssi.
    #[test]
    fn windowed_matches_full() -> Result<(), Box<dyn std::error::Error>> {
        for name in ["3_1", "4_1", "5_1", "5_2a", "6_2a", "6_3", "7_6a"] {
            let l = InvLink::load(name)?;
            let full = ssi_invariant_h1::<F>(&l, false, SymBuildConfig::default());

            // tight window at the true divisibility d = (ss − w + r − 1) / 2 (max over the two classes).
            let (w, r) = (l.writhe(), l.seifert_circles().len() as i32);
            let d_max = [(full.0 - w + r - 1) / 2, (full.1 - w + r - 1) / 2].into_iter().max().unwrap();
            let win = ssi_invariant_h1_windowed::<F>(&l, false, SymBuildConfig::default(), Some(d_max));

            assert_eq!(win, full, "{name} (d_max = {d_max})");
        }
        Ok(())
    }

    // Gate for the windowed path: the 44-crossing Whitehead double of P(−3,3,−3). Two variants for
    // a runtime comparison: `no_full_deloop = false` deloops+eliminates at the final step (heavier
    // build, lighter solve); `true` defers that to the matrix-level expansion (lighter build, heavier
    // solve). Both must give (0, 2).
    fn run_wh_pretzel_3_windowed(no_full_deloop: bool) {
        use crate::tng::builder::{BuildMode, CutOption};
        let _ = env_logger::Builder::from_default_env().target(env_logger::Target::Stdout).try_init();

        let k = InvLink::sym_pretzel(-3, 3, -3);
        let w = k.whitehead_double(true, 0);
        let config = SymBuildConfig {
            mode: BuildMode::MinFill,
            cut: CutOption::Auto(2),
            max_elim_cost: Some(1 << 16),
            no_full_deloop,
            ..Default::default()
        };
        let ssi = ssi_invariant_h1_windowed::<F>(&w, false, config, None);
        println!("ssi(Wh+(P(-3,3,-3))) [windowed, no_full_deloop={no_full_deloop}] = {ssi:?}");
        assert_eq!(ssi, (0, 2));
    }

    #[test]
    #[ignore = "heavy: 44-crossing Whitehead double via windowed H=1"]
    fn ssi_wh_pretzel_3_windowed() {
        run_wh_pretzel_3_windowed(false);
    }

    #[test]
    #[ignore = "heavy: 44-crossing Whitehead double via windowed H=1 (no_full_deloop)"]
    fn ssi_wh_pretzel_3_windowed_nodeloop() {
        run_wh_pretzel_3_windowed(true);
    }

    // Minimal reproducer for the windowed canon-transport bug: Wh(6_2a) (30 crossings) truncates a
    // canon cycle to zero, while Wh(6_1a) (also 30) works. Root cause: the build-time q-filter moves
    // the canon cycle's representative onto high-q vertices that the window then drops. See notes.
    // Wh(6_2a) (30 crossings) needs the adaptive window widening: at max_div=3 the reduction pushes
    // the canon representative above the window (truncated to zero), so the driver widens until it
    // reappears. Should converge to ssi = (2, 2).
    #[test]
    #[ignore = "heavy-ish (~1 min): windowed Wh(6_2a) with adaptive widening"]
    fn wh_6_2a_windowed_adaptive() {
        use crate::tng::builder::{BuildMode, CutOption};
        let l = InvLink::load("6_2a").unwrap().whitehead_double(true, 0);
        let nx = l.inner().n_crossings();
        let cfg = SymBuildConfig { mode: BuildMode::MinFill, cut: CutOption::AtCrossings(vec![nx / 2]), max_elim_cost: Some(1 << 16), no_full_deloop: false, ..Default::default() };
        let ssi = ssi_invariant_h1_windowed::<F>(&l, false, cfg, None);
        assert_eq!(ssi, (2, 2));
    }

    #[test]
    fn k3_1_m() {
        let l = InvLink::test_data("3_1").mirror();
        assert_eq!(ssi_h1(&l, false), (-2, -2));
    }

    #[test]
    fn k3_1_red() {
        let l = InvLink::test_data("3_1");
        assert_eq!(ssi_h1(&l, true), (2, 2));
    }

    #[test]
    fn k9_46() {
        let l = InvLink::from_symmetric_pd_code(
            [[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]
        );
        assert_eq!(ssi_h1(&l, false), (0, 2));
    }

    // h1 must agree with the divisibility pipeline (`ssi_invariant`), incl. chunked builds.
    #[test]
    fn matches_ssi_invariant_chunked() -> Result<(), Box<dyn std::error::Error>> {
        use yui_core::poly::Poly;
        use crate::tng::builder::CutOption;
        use crate::khi::ssi_invariant;

        type P = Poly<'H', F>;
        let c = P::variable();

        for name in ["3_1", "4_1", "6_1a"] {
            let l = InvLink::load(name)?;
            let config = SymBuildConfig { cut: CutOption::Auto(3), ..Default::default() };
            let expected = ssi_invariant(&l, &c, false, config.clone());
            let h1 = ssi_invariant_h1::<F>(&l, false, config);
            assert_eq!(h1, expected, "{name}");
        }
        Ok(())
    }
}
