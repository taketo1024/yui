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
// The `inv_flip_*` constructors at the bottom hold the *output* diagrams of that construction,
// baked from the campaign runs so the tests don't re-derive them (and so the crossing order —
// hence the build — is fixed).
// They are 80 to 112 crossings; the last four take minutes even on the reduced pipeline.
// Run logs: `logs/flip_sum_*.log` (kept out of git; see .gitignore).
// ---------------------------------------------------------------------------------------------

macro_rules! prop_1_5 {
    ($(#[$m:meta])* $name:ident, $expected:expr, $knot:ident) => {
        $(#[$m])*
        #[test]
        fn $name() {
            init_logger();
            assert_ssi(&$knot(), cut_config(CutOption::Auto(3)), Side::Positive, $expected);
        }
    }
}

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*9_46a # -2*9_46a) = (-2, 0)"]
    prop_1_5_9_46a, (-2, 0), inv_flip_9_46a
);

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*9_46b # -2*9_46b) = (-2, 0)"]
    prop_1_5_9_46b, (-2, 0), inv_flip_9_46b
);

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*8_21b # -2*8_21b) = (0, 2)"]
    prop_1_5_8_21b, (0, 2), inv_flip_8_21b
);

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*10_160 # -2*10_160) = (-2, 0)"]
    prop_1_5_10_160, (-2, 0), inv_flip_10_160
);

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*10_162 # -2*10_162) = (-2, 0); heavy, ~minutes"]
    prop_1_5_10_162, (-2, 0), inv_flip_10_162
);

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*10_165 # -2*10_165) = (-2, 0); heavy, ~minutes"]
    prop_1_5_10_165, (-2, 0), inv_flip_10_165
);

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*10_141b # -2*10_141b) = (0, 2); heavy, ~minutes (112 crossings)"]
    prop_1_5_10_141b, (0, 2), inv_flip_10_141b
);

prop_1_5!(
    #[ignore = "Prop 1.5: ssi(2*10_144b # -2*10_144b) = (0, 2); heavy, ~minutes (96 crossings)"]
    prop_1_5_10_144b, (0, 2), inv_flip_10_144b
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

/// `2K # (-2K)` output diagram of `9_46a` (40 crossings), baked from the campaign
/// run so the crossing order — hence the build — is fixed.
fn inv_flip_9_46a() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [1,6,2,7],[4,78,5,77],[5,10,6,11],[8,74,9,73],[9,2,10,3],[12,20,13,19],[13,28,14,29],[16,12,17,11],[17,24,18,25],[20,16,21,15],[22,30,23,29],
        [23,18,24,19],[26,22,27,21],[27,14,28,15],[30,26,31,25],[31,36,32,37],[34,48,35,47],[35,40,36,41],[38,44,39,43],[39,32,40,33],[41,46,42,47],
        [44,38,45,37],[45,50,46,51],[48,34,49,33],[49,42,50,43],[52,60,53,59],[53,68,54,69],[56,52,57,51],[57,64,58,65],[60,56,61,55],[62,70,63,69],
        [63,58,64,59],[66,62,67,61],[67,54,68,55],[70,66,71,65],[71,76,72,77],[74,8,75,7],[75,80,76,1],[78,4,79,3],[79,72,80,73]
    ])
}

/// `2K # (-2K)` output diagram of `9_46b` (40 crossings), baked from the campaign
/// run so the crossing order — hence the build — is fixed.
fn inv_flip_9_46b() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [3,9,4,8],[4,77,5,78],[6,71,7,72],[9,3,10,2],[10,75,11,76],[11,27,12,26],[12,19,13,20],[15,31,16,30],[17,25,18,24],[18,13,19,14],[21,17,22,16],
        [22,29,23,30],[25,21,26,20],[27,15,28,14],[28,23,29,24],[33,39,34,38],[34,47,35,48],[36,41,37,42],[39,33,40,32],[40,45,41,46],[43,49,44,48],
        [44,37,45,38],[46,31,47,32],[49,43,50,42],[50,35,51,36],[51,67,52,66],[52,59,53,60],[55,71,56,70],[57,65,58,64],[58,53,59,54],[61,57,62,56],
        [62,69,63,70],[65,61,66,60],[67,55,68,54],[68,63,69,64],[73,79,74,78],[74,7,75,8],[76,1,77,2],[79,73,80,72],[80,5,1,6]
    ])
}

/// `2K # (-2K)` output diagram of `8_21b` (40 crossings), baked from the campaign
/// run so the crossing order — hence the build — is fixed.
fn inv_flip_8_21b() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [3,9,4,8],[4,2,5,1],[5,76,6,77],[9,3,10,2],[10,8,11,7],[11,14,12,15],[12,19,13,20],[16,26,17,25],[17,20,18,21],[18,13,19,14],[21,24,22,25],
        [22,29,23,30],[26,16,27,15],[27,30,28,31],[28,23,29,24],[33,39,34,38],[34,32,35,31],[35,46,36,47],[39,33,40,32],[40,38,41,37],[43,49,44,48],
        [44,42,45,41],[45,36,46,37],[49,43,50,42],[50,48,51,47],[51,54,52,55],[52,59,53,60],[56,66,57,65],[57,60,58,61],[58,53,59,54],[61,64,62,65],
        [62,69,63,70],[66,56,67,55],[67,70,68,71],[68,63,69,64],[73,79,74,78],[74,72,75,71],[75,6,76,7],[79,73,80,72],[80,78,1,77]
    ])
}

/// `2K # (-2K)` output diagram of `10_160` (40 crossings), baked from the campaign
/// run so the crossing order — hence the build — is fixed.
fn inv_flip_10_160() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [1,9,2,8],[3,79,4,78],[6,76,7,75],[9,3,10,2],[12,19,13,20],[15,26,16,27],[18,23,19,24],[20,13,21,14],[22,29,23,30],[24,12,25,11],[25,16,26,17],
        [27,14,28,15],[28,21,29,22],[30,18,31,17],[31,44,32,45],[33,41,34,40],[34,48,35,47],[36,46,37,45],[37,50,38,51],[39,33,40,32],[41,49,42,48],
        [43,39,44,38],[46,36,47,35],[49,43,50,42],[52,59,53,60],[55,66,56,67],[58,63,59,64],[60,53,61,54],[62,69,63,70],[64,52,65,51],[65,56,66,57],
        [67,54,68,55],[68,61,69,62],[70,58,71,57],[71,4,72,5],[73,1,74,80],[74,8,75,7],[76,6,77,5],[77,10,78,11],[79,73,80,72]
    ])
}

/// `2K # (-2K)` output diagram of `10_162` (40 crossings), baked from the campaign
/// run so the crossing order — hence the build — is fixed.
fn inv_flip_10_162() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [2,9,3,10],[4,78,5,77],[6,1,7,2],[7,75,8,74],[10,3,11,4],[11,19,12,18],[14,27,15,28],[15,21,16,20],[17,24,18,25],[19,13,20,12],[21,27,22,26],
        [23,31,24,30],[25,16,26,17],[28,13,29,14],[29,23,30,22],[32,39,33,40],[33,49,34,48],[36,46,37,45],[38,31,39,32],[40,35,41,36],[42,49,43,50],
        [44,38,45,37],[46,41,47,42],[47,35,48,34],[50,43,51,44],[51,59,52,58],[54,67,55,68],[55,61,56,60],[57,64,58,65],[59,53,60,52],[61,67,62,66],
        [63,71,64,70],[65,56,66,57],[68,53,69,54],[69,63,70,62],[72,79,73,80],[73,9,74,8],[76,6,77,5],[78,71,79,72],[80,75,1,76]
    ])
}

/// `2K # (-2K)` output diagram of `10_165` (40 crossings), baked from the campaign
/// run so the crossing order — hence the build — is fixed.
fn inv_flip_10_165() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [1,6,2,7],[4,78,5,77],[5,10,6,11],[7,74,8,75],[9,2,10,3],[11,17,12,16],[13,29,14,28],[15,21,16,20],[18,23,19,24],[19,13,20,12],[21,27,22,26],
        [24,17,25,18],[25,31,26,30],[27,15,28,14],[29,23,30,22],[31,36,32,37],[33,48,34,49],[35,40,36,41],[38,44,39,43],[39,32,40,33],[41,46,42,47],
        [44,38,45,37],[45,50,46,51],[47,34,48,35],[49,42,50,43],[51,57,52,56],[53,69,54,68],[55,61,56,60],[58,63,59,64],[59,53,60,52],[61,67,62,66],
        [64,57,65,58],[65,71,66,70],[67,55,68,54],[69,63,70,62],[71,76,72,77],[73,8,74,9],[75,80,76,1],[78,4,79,3],[79,72,80,73]
    ])
}

/// `2K # (-2K)` output diagram of `10_141b` (56 crossings), baked from the campaign
/// run so the crossing order — hence the build — is fixed.
fn inv_flip_10_141b() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [3,13,4,12],[4,2,5,1],[5,108,6,109],[6,9,7,10],[7,106,8,107],[13,3,14,2],[14,12,15,11],[15,18,16,19],[16,27,17,28],[22,36,23,35],[23,21,24,20],
        [24,34,25,33],[25,28,26,29],[26,17,27,18],[29,32,30,33],[30,41,31,42],[36,22,37,21],[37,35,38,34],[38,20,39,19],[39,42,40,43],[40,31,41,32],
        [45,55,46,54],[46,44,47,43],[47,66,48,67],[48,51,49,52],[49,64,50,65],[55,45,56,44],[56,54,57,53],[59,69,60,68],[60,58,61,57],[61,52,62,53],
        [62,65,63,66],[63,50,64,51],[69,59,70,58],[70,68,71,67],[71,74,72,75],[72,83,73,84],[78,92,79,91],[79,77,80,76],[80,90,81,89],[81,84,82,85],
        [82,73,83,74],[85,88,86,89],[86,97,87,98],[92,78,93,77],[93,91,94,90],[94,76,95,75],[95,98,96,99],[96,87,97,88],[101,111,102,110],[102,100,103,99],
        [103,10,104,11],[104,107,105,108],[105,8,106,9],[111,101,112,100],[112,110,1,109]
    ])
}

/// `2K # (-2K)` output diagram of `10_144b` (48 crossings), baked from the campaign
/// run so the crossing order — hence the build — is fixed.
fn inv_flip_10_144b() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [3,11,4,10],[4,2,5,1],[5,92,6,93],[7,90,8,91],[11,3,12,2],[12,10,13,9],[13,16,14,17],[14,23,15,24],[18,32,19,31],[20,30,21,29],[21,24,22,25],
        [22,15,23,16],[25,28,26,29],[26,35,27,36],[30,20,31,19],[32,18,33,17],[33,36,34,37],[34,27,35,28],[39,47,40,46],[40,38,41,37],[41,56,42,57],
        [43,54,44,55],[47,39,48,38],[48,46,49,45],[51,59,52,58],[52,50,53,49],[53,44,54,45],[55,42,56,43],[59,51,60,50],[60,58,61,57],[61,64,62,65],
        [62,71,63,72],[66,80,67,79],[68,78,69,77],[69,72,70,73],[70,63,71,64],[73,76,74,77],[74,83,75,84],[78,68,79,67],[80,66,81,65],[81,84,82,85],
        [82,75,83,76],[87,95,88,94],[88,86,89,85],[89,8,90,9],[91,6,92,7],[95,87,96,86],[96,94,1,93]
    ])
}
