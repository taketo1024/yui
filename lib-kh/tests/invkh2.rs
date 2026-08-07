//! Computational results of
//!
//! > T. Sano, "Involutive Khovanov homology and equivariant knots II" (2026).
//! > `references/Sano2026_InvKh2.pdf`
//!
//! One test per claim, named after it, asserting the published value. The point is that the
//! paper can be checked against the program: pick a statement, find the test with that name,
//! read the expected `(s̲, s̄)` in the `assert_eq!`.
//!
//! Everything here is `#[ignore]`d — these are minutes-to-hours runs, not a test suite. Run one by
//! name, e.g.
//!
//! ```text
//! cargo test -r -p yui-kh --test invkh2 -- --ignored --nocapture prop_1_3_9_46a
//! RUST_LOG=info cargo test -r -p yui-kh --test invkh2 -- --ignored --nocapture prop_1_7_wh_p3
//! ```
//!
//! `--features big-link` is needed wherever a diagram exceeds 63 crossings (marked below).
//! To list the whole catalogue with its expected values:
//!
//! ```text
//! cargo test -p yui-kh --test invkh2 -- --list
//! ```
//!
//! Conventions. `(s̲, s̄)` is the equivariant Rasmussen invariant of §5, over `𝔽₂`. Since
//! `(s̲(−K), s̄(−K)) = (−s̄(K), −s̲(K))`, only one side of the chirality is computed; the helpers
//! below take the positive side. Knot names follow [Lam25] — `K_a` / `K_b` are the two
//! inversion classes of a knot with two strong inversions.

use yui_kh::tng::builder::{Strategy, CutOption, SymBuildConfig};
use yui_link::InvLink;

mod common;
use common::*;

// The production Wh⁺-pretzel runs, reproducible on any machine (the in-memory construction
// fixes the crossing order; a PD re-import would reorder under MinCut).
fn wh_pretzel_config(cut_at: usize) -> SymBuildConfig {
    SymBuildConfig {
        strategy: Strategy::MinFill,
        cut: CutOption::At(vec![cut_at]),
        max_elim_cost: Some(1 << 16),
        no_full_deloop: true,
        ..Default::default()
    }
}

// ---------------------------------------------------------------------------------------------
// Proposition 1.3 — within the dataset of [Lam25], these eight are the only strongly invertible
// knots with (s̲, s̄) ≠ (s, s). The paper's table, with `s` the 𝔽₂-Rasmussen invariant:
//
//     K          (s̲, s̄)     s          K          (s̲, s̄)     s
//     8_21b      ( 2,  4)    2          10_144b    ( 2,  4)    2
//     9_46a      (-2,  0)    0          10_160     ( 2,  4)    4
//     9_46b      (-2,  0)    0          10_162     (-4, -2)   -2
//     10_141b    ( 0,  2)    0          10_165     (-4, -2)   -2
//
// The `a`/`b` partners not listed here have (s̲, s̄) = (s, s): 8_21a = (2, 2),
// 10_141a = (0, 0), 10_144a = (2, 2). PD codes are the doubly-transvergent diagrams of
// [Lam25, Appendix B].
//
// Every test states the paper's own pair. `Side::Positive` computes whichever of `K`, `K*` has
// non-negative writhe — the cheaper build, since the window is `-n_neg ..= 1` — and converts back
// by Prop 1.3 of InvKh1.
// ---------------------------------------------------------------------------------------------

#[test]
#[ignore = "Prop 1.3: ssi(8_21b) = (2, 4), while s = 2"]
fn prop_1_3_8_21b() {
    init_logger();
    assert_ssi(&inv_8_21b(), SymBuildConfig::default(), Side::Positive, (2, 4));
}

#[test]
#[ignore = "Prop 1.3: ssi(9_46a) = (-2, 0), while s = 0"]
fn prop_1_3_9_46a() {
    init_logger();
    assert_ssi(&inv_9_46a(), SymBuildConfig::default(), Side::Positive, (-2, 0));
}

#[test]
#[ignore = "Prop 1.3: ssi(9_46b) = (-2, 0), while s = 0"]
fn prop_1_3_9_46b() {
    init_logger();
    assert_ssi(&inv_9_46b(), SymBuildConfig::default(), Side::Positive, (-2, 0));
}

#[test]
#[ignore = "Prop 1.3: ssi(10_141b) = (0, 2), while s = 0"]
fn prop_1_3_10_141b() {
    init_logger();
    assert_ssi(&inv_10_141b(), SymBuildConfig::default(), Side::Positive, (0, 2));
}

#[test]
#[ignore = "Prop 1.3: ssi(10_144b) = (2, 4), while s = 2"]
fn prop_1_3_10_144b() {
    init_logger();
    assert_ssi(&inv_10_144b(), SymBuildConfig::default(), Side::Positive, (2, 4));
}

#[test]
#[ignore = "Prop 1.3: ssi(10_160) = (2, 4), while s = 4"]
fn prop_1_3_10_160() {
    init_logger();
    assert_ssi(&inv_10_160(), SymBuildConfig::default(), Side::Positive, (2, 4));
}

#[test]
#[ignore = "Prop 1.3: ssi(10_162) = (-4, -2), while s = -2"]
fn prop_1_3_10_162() {
    init_logger();
    assert_ssi(&inv_10_162(), SymBuildConfig::default(), Side::Positive, (-4, -2));
}

#[test]
#[ignore = "Prop 1.3: ssi(10_165) = (-4, -2), while s = -2"]
fn prop_1_3_10_165() {
    init_logger();
    assert_ssi(&inv_10_165(), SymBuildConfig::default(), Side::Positive, (-4, -2));
}

// ---------------------------------------------------------------------------------------------
// Proposition 1.4 — the equivariant connected sum `K_a # −K_b := (K # −K, τ_a # −τ_b)` of the two
// inversion classes. Within the dataset of [Lam25], these three are the only such sums with
// non-trivial equivariant Rasmussen invariant, and all three equal (-2, 0).
//
// Each is a 16- to 20-crossing diagram, so these run in seconds to minutes.
// ---------------------------------------------------------------------------------------------

// `K_a # −K_b`, computed on the positive side of the chirality.
fn assert_conn_sum_ssi(a: InvLink, b: InvLink, expected: (i32, i32)) {
    let k = a.conn_sum(&b.mirror());
    println!("conn-sum: {} crossings, writhe {}", k.inner().n_crossings(), k.writhe());
    assert_ssi(&k, cut_config(CutOption::Auto(2)), Side::Positive, expected);
}

#[test]
#[ignore = "Prop 1.4: ssi(8_21a # -8_21b) = (-2, 0)"]
fn prop_1_4_8_21() {
    init_logger();
    assert_conn_sum_ssi(inv_8_21a(), inv_8_21b(), (-2, 0));
}

#[test]
#[ignore = "Prop 1.4: ssi(10_141a # -10_141b) = (-2, 0)"]
fn prop_1_4_10_141() {
    init_logger();
    assert_conn_sum_ssi(inv_10_141a(), inv_10_141b(), (-2, 0));
}

#[test]
#[ignore = "Prop 1.4: ssi(10_144a # -10_144b) = (-2, 0)"]
fn prop_1_4_10_144() {
    init_logger();
    assert_conn_sum_ssi(inv_10_144a(), inv_10_144b(), (-2, 0));
}

// ---------------------------------------------------------------------------------------------
// Proposition 1.5 — the `2K # (−2K)` construction: `K̃ := (K # K^r) # (−K # −K)` with the strong
// inversion `τ_f # (−τ # −τ)`, where `τ_f` is the flipping involution on `K # K^r` (Figure 1).
// The eight knots of Proposition 1.3 are exactly those for which `K̃` is a strongly invertible
// slice knot with non-trivial equivariant Rasmussen invariant, each equal to (-2, 0) or (0, 2).
//
// `dms_construction` builds `K̃` from the source knot at run time, so the diagram is derived here
// rather than pasted in. The results are 40 to 56 crossings; the last two take minutes even on the
// reduced pipeline.
// ---------------------------------------------------------------------------------------------

macro_rules! prop_1_5 {
    ($(#[$m:meta])* $name:ident, $expected:expr, $knot:ident) => {
        $(#[$m])*
        #[test]
        fn $name() {
            init_logger();
            let k = dms_construction($knot());
            println!("K~: {} crossings, writhe {}", k.n_crossings(), k.writhe());
            assert_ssi(&k, cut_config(CutOption::Auto(3)), Side::Positive, $expected);
        }
    }
}

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*9_46a # -2*9_46a) = (-2, 0)"]
    prop_1_5_9_46a, (-2, 0), inv_9_46a
);

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*9_46b # -2*9_46b) = (-2, 0)"]
    prop_1_5_9_46b, (-2, 0), inv_9_46b
);

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*8_21b # -2*8_21b) = (0, 2)"]
    prop_1_5_8_21b, (0, 2), inv_8_21b
);

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*10_160 # -2*10_160) = (-2, 0)"]
    prop_1_5_10_160, (-2, 0), inv_10_160
);

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*10_162 # -2*10_162) = (-2, 0); heavy, ~minutes"]
    prop_1_5_10_162, (-2, 0), inv_10_162
);

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*10_165 # -2*10_165) = (-2, 0); heavy, ~minutes"]
    prop_1_5_10_165, (-2, 0), inv_10_165
);

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*10_141b # -2*10_141b) = (0, 2); heavy, ~minutes (56 crossings)"]
    prop_1_5_10_141b, (0, 2), inv_10_141b
);

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*10_144b # -2*10_144b) = (0, 2); heavy, ~minutes (48 crossings)"]
    prop_1_5_10_144b, (0, 2), inv_10_144b
);

// ---------------------------------------------------------------------------------------------
// Theorem 4 — for every odd `p ≥ 3`, the strongly invertible pretzel knot `K = P(-p, p, -p)` has
// `(s̲(K), s̄(K)) = (0, 2)`. Proved by hand in §6 (no computer), so these are a check on the
// program, not evidence for the theorem. `p = 3` gives `m(9_46)`.
//
// `p = 3, 5` here; `p = 7` (21 crossings) is left out as it adds nothing beyond `p = 5`.
// ---------------------------------------------------------------------------------------------

#[test]
#[ignore = "Thm 4, p=3: ssi(P(-3, 3, -3)) = (0, 2)"]
fn thm_4_pretzel_3() {
    init_logger();
    let k = InvLink::sym_pretzel(-3, 3, -3);
    assert_eq!(ssi(&k), (0, 2));
}

#[test]
#[ignore = "Thm 4, p=5: ssi(P(-5, 5, -5)) = (0, 2)"]
fn thm_4_pretzel_5() {
    init_logger();
    let k = InvLink::sym_pretzel(-5, 5, -5);
    assert_eq!(ssi(&k), (0, 2));
}

// ---------------------------------------------------------------------------------------------
// Proposition 1.7 — the Whitehead doubles of the slice pretzel knots `P(-3, 3, -3)` and
// `P(-5, 5, -5)` have non-trivial equivariant Rasmussen invariant, both equal to (0, 2).
//
// With Theorem 4 and [Gut+23, Prop 2.1] + [CP21, Thm 1.2] / [Hay21, Thm 2.2], this proves
// Theorem 1: `Wh(D)` and `τWh(D)` are topologically but not smoothly isotopic — an exotic pair.
//
// Cost. `Wh(P(-3,3,-3))` is 44 crossings and takes about a minute. `Wh(P(-5,5,-5))` is 72
// crossings and needs `--features big-link`; the paper reports 17 hours on a 64-core machine with
// 3 TB of memory. Do not run it casually. Run logs: `logs/wh_pretzel_3.log`, `logs/wh-p5*.log`.
// ---------------------------------------------------------------------------------------------

#[test]
#[ignore = "Prop 1.7: ssi(Wh+(P(-3,3,-3))) = (0, 2); heavy, 44 crossings, ~1 min"]
fn prop_1_7_wh_p3() {
    init_logger();
    let k = InvLink::sym_pretzel(-3, 3, -3);
    let w = InvLink::whitehead_double(&k, true, 0);
    assert_ssi(&w, wh_pretzel_config(17), Side::Positive, (0, 2));
}

#[cfg(feature = "big-link")]
#[test]
#[ignore = "Prop 1.7: ssi(Wh+(P(-5,5,-5))) = (0, 2); 72 crossings — 17 h on 64 cores / 3 TB"]
fn prop_1_7_wh_p5() {
    init_logger();
    let k = InvLink::sym_pretzel(-5, 5, -5);
    let w = InvLink::whitehead_double(&k, true, 0);
    assert_ssi(&w, wh_pretzel_config(20), Side::Positive, (0, 2));
}

// ---------------------------------------------------------------------------------------------
// Knot data — the doubly-transvergent diagrams of [Lam25, Appendix B], with the standard strong
// inversion installed by `from_symmetric_pd_code`. `K_a` / `K_b` are the two inversion classes.
// ---------------------------------------------------------------------------------------------

fn inv_8_21a() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [1,12,20,11],[20,16,19,17],[19,3,18,4],[16,12,15,13],[15,7,14,8],[11,2,10,1],[10,6,9,7],[9,13,8,14],[6,2,5,3],[5,17,4,18]
    ])
}

fn inv_8_21b() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [1,17,20,18],[20,12,19,13],[16,7,15,6],[15,11,14,12],[14,18,13,19],[11,7,10,8],[10,2,9,3],[6,17,5,16],[5,1,4,2],[4,8,3,9]
    ])
}

fn inv_9_46a() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [20,13,19,12],[19,3,18,4],[16,1,15,20],[15,7,14,8],[12,17,11,16],[10,3,9,2],[9,13,8,14],[6,11,5,10],[5,17,4,18],[2,7,1,6]
    ])
}

fn inv_9_46b() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [3,9,4,8],[4,17,5,18],[6,11,7,12],[9,3,10,2],[10,15,11,16],[13,19,14,18],[14,7,15,8],[16,1,17,2],[19,13,20,12],[20,5,1,6]
    ])
}

fn inv_10_141b() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [1,25,28,26],[28,16,27,17],[22,9,21,8],[21,24,20,23],[20,11,19,10],[19,15,18,16],[18,26,17,27],[15,11,14,12],[14,2,13,3],[8,23,7,22],[7,10,6,9],[6,25,5,24],[5,1,4,2],[4,12,3,13]
    ])
}

fn inv_10_144b() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [1,21,24,22],[24,14,23,15],[20,7,19,6],[18,9,17,8],[17,13,16,14],[16,22,15,23],[13,9,12,10],[12,2,11,3],[8,19,7,18],[6,21,5,20],[5,1,4,2],[4,10,3,11]
    ])
}

fn inv_10_160() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [2,10,3,9],[5,17,6,16],[8,14,9,13],[10,4,11,3],[12,20,13,19],[14,1,15,2],[15,7,16,6],[17,5,18,4],[18,12,19,11],[20,7,1,8]
    ])
}

fn inv_10_162() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [18,4,19,3],[19,12,20,13],[1,8,2,9],[4,18,5,17],[5,10,6,11],[7,15,8,14],[9,2,10,3],[11,16,12,17],[13,20,14,1],[15,7,16,6]
    ])
}

fn inv_10_165() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [8,14,9,13],[9,2,10,3],[11,16,12,17],[14,8,15,7],[15,20,16,1],[17,4,18,5],[19,12,20,13],[1,6,2,7],[3,18,4,19],[5,10,6,11]
    ])
}

fn inv_10_141a() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [1,16,28,15],[28,3,27,2],[27,18,26,17],[26,22,25,23],[25,5,24,6],[22,18,21,19],[21,9,20,10],[15,2,14,1],[14,17,13,16],[13,4,12,3],[12,8,11,9],[11,19,10,20],[8,4,7,5],[7,23,6,24]
    ])
}

fn inv_10_144a() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [24,15,23,14],[23,19,22,20],[22,4,21,5],[19,15,18,16],[18,8,17,9],[14,1,13,24],[12,3,11,2],[11,7,10,8],[10,16,9,17],[7,3,6,4],[6,20,5,21],[2,13,1,12]
    ])
}

// `2K # (−2K)` as Dai–Mallick–Stoffregen define it (`references/DMS.pdf` §1.2,
// https://arxiv.org/abs/2201.01875): the flip sum `K # r(K)`, then the mirror of the flip sum as
// a single summand. `conn_sum` reads the splice edges off the axis, so the summands are the only
// input. The mirrored convention of the other paper gives `−ssi`.
pub fn dms_construction(k: InvLink) -> InvLink {
    let mk2 = InvLink::si_knot_from(k.inner().conn_sum(k.inner()).mirror()); // m(K # K) = m(K) # m(K)
    k.conn_sum(&k.reversed()).conn_sum(&mk2)
}
