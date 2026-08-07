//! Computational results of
//!
//! > T. Sano, "Involutive Khovanov homology and equivariant knots" (2025).
//! > `references/Sano2025_InvKh.pdf`, arXiv:2404.08568
//!
//! One test per claim, named after it, asserting the published value. See `invkh2.rs` for the
//! sequel's results and `invkh3.rs` for the interlocks.
//!
//! Everything here is `#[ignore]`d and run on demand:
//!
//! ```text
//! cargo test -r -p yui-kh --test invkh1 -- --ignored --nocapture
//! cargo test -p yui-kh --test invkh1 -- --list          # the catalogue, with expected values
//! ```
//!
//! `(s̲, s̄)` is the equivariant Rasmussen invariant over `𝔽₂` (§3), and `s` is the ordinary
//! `𝔽₂`-Rasmussen invariant. PD codes are at the bottom of the file.

use yui_kh::tng::builder::SymBuildConfig;
use yui_link::InvLink;

mod common;
use common::*;

// ---------------------------------------------------------------------------------------------
// Proposition 1.2 — the three strongly invertible slice knots of [HS24]:
//
//     K = m(9_46), 15n_103488, 17nh_73    have    s̲(K) = 0 < 2 = s̄(K).
//
// Combined with Corollary 1.11 this recovers that each admits a pair of non-smoothly-isotopic
// slice disks. `17nh_73` is the positron knot, and is `J_0` of Theorem 1.
// ---------------------------------------------------------------------------------------------

#[test]
#[ignore = "Prop 1.2: ssi(m(9_46)) = (0, 2)"]
fn prop_1_2_m9_46() {
    init_logger();
    assert_ssi(&inv_9_46(), SymBuildConfig::default(), Side::Positive, (0, 2));
}

#[test]
#[ignore = "Prop 1.2: ssi(15n_103488) = (0, 2)"]
fn prop_1_2_15n_103488() {
    init_logger();
    assert_ssi(&inv_15n_103488(), SymBuildConfig::default(), Side::Positive, (0, 2));
}

#[test]
#[ignore = "Prop 1.2: ssi(17nh_73) = (0, 2); J_0 of Theorem 1"]
fn prop_1_2_17nh_73() {
    init_logger();
    assert_ssi(&inv_17nh_73(), SymBuildConfig::default(), Side::Positive, (0, 2));
}

// ---------------------------------------------------------------------------------------------
// Proposition 1.3 — behaviour under mirroring: for the mirror `K*` of `K`,
//
//     s̲(K*) = -s̄(K)        (equivalently the pair maps by (s̲, s̄) ↦ (-s̄, -s̲)).
//
// Checked as a relation rather than against a stored number: compute both sides and compare.
// This is what licenses computing every other family on its positive-writhe representative.
// ---------------------------------------------------------------------------------------------

fn assert_mirror_relation(l: &InvLink) {
    let (s_lo, s_hi) = ssi(l);
    let (m_lo, m_hi) = ssi(&l.mirror());
    assert_eq!(m_lo, -s_hi, "s̲(K*) = -s̄(K)");
    assert_eq!(m_hi, -s_lo, "s̄(K*) = -s̲(K)");
}

#[test]
#[ignore = "Prop 1.3: (s̲, s̄)(K*) = (-s̄, -s̲)(K) on m(9_46)"]
fn prop_1_3_mirror_9_46() {
    init_logger();
    assert_mirror_relation(&inv_9_46());
}

#[test]
#[ignore = "Prop 1.3: (s̲, s̄)(K*) = (-s̄, -s̲)(K) on 15n_103488"]
fn prop_1_3_mirror_15n_103488() {
    init_logger();
    assert_mirror_relation(&inv_15n_103488());
}

#[test]
#[ignore = "Prop 1.3: (s̲, s̄)(K*) = (-s̄, -s̲)(K) on the 7_6a control (ssi-trivial)"]
fn prop_1_3_mirror_7_6a() {
    init_logger();
    assert_mirror_relation(&inv_7_6a());
}

// ---------------------------------------------------------------------------------------------
// Proposition 1.4 — subadditivity under the equivariant connected sum:
//
//     s̲(K) + s̲(K') ≤ s̲(K # K') ≤ s̲(K) + s̄(K') ≤ s̄(K # K') ≤ s̄(K) + s̄(K').
//
// Checked as a chain of inequalities on computed values, so it tests the invariant's algebra
// rather than one more knot.
// ---------------------------------------------------------------------------------------------

fn assert_conn_sum_bounds(a: &InvLink, b: &InvLink) {
    let (a_lo, a_hi) = ssi(a);
    let (b_lo, b_hi) = ssi(b);
    let k = a.conn_sum(b);
    let (k_lo, k_hi) = ssi(&k);

    assert!(a_lo + b_lo <= k_lo, "s̲(K) + s̲(K') ≤ s̲(K # K'): {a_lo} + {b_lo} ≤ {k_lo}");
    assert!(k_lo <= a_lo + b_hi, "s̲(K # K') ≤ s̲(K) + s̄(K'): {k_lo} ≤ {a_lo} + {b_hi}");
    assert!(a_lo + b_hi <= k_hi, "s̲(K) + s̄(K') ≤ s̄(K # K'): {a_lo} + {b_hi} ≤ {k_hi}");
    assert!(k_hi <= a_hi + b_hi, "s̄(K # K') ≤ s̄(K) + s̄(K'): {k_hi} ≤ {a_hi} + {b_hi}");
}

#[test]
#[ignore = "Prop 1.4: conn-sum bounds on m(9_46) # m(9_46)"]
fn prop_1_4_bounds_9_46_9_46() {
    init_logger();
    assert_conn_sum_bounds(&inv_9_46(), &inv_9_46());
}

#[test]
#[ignore = "Prop 1.4: conn-sum bounds on m(9_46) # 7_6a"]
fn prop_1_4_bounds_9_46_7_6a() {
    init_logger();
    assert_conn_sum_bounds(&inv_9_46(), &inv_7_6a());
}

// ---------------------------------------------------------------------------------------------
// Proposition 1.7 — the positive `(p, q)`-torus knot has `s̲ = s̄ = (p-1)(q-1)` with respect to
// its unique inverting involution. Only `T_{2,3}` is reachable with the current constructors:
// `P(-1,-1,-1)` is the positive trefoil (`writhe(P(a,b,a)) = -(2a+b) = 3`), so the prediction is
// `(2, 2)`.
//
// TODO: a general `InvLink::torus(p, q)` from the braid `(σ₁⋯σ_{p-1})^q` would let the rest of
// the family be tested; the work is picking the transvergent diagram, not closing the braid.
// ---------------------------------------------------------------------------------------------

#[test]
#[ignore = "Prop 1.7: ssi(T_{2,3}) = (2, 2) = ((p-1)(q-1), (p-1)(q-1))"]
fn prop_1_7_torus_2_3() {
    init_logger();
    let k = InvLink::sym_pretzel(-1, -1, -1);
    assert_eq!(k.writhe(), 3, "P(-1,-1,-1) should be the positive trefoil");
    assert_ssi(&k, SymBuildConfig::default(), Side::Positive, (2, 2));
}

// ---------------------------------------------------------------------------------------------
// Knot data — symmetric PD codes, installing the standard strong inversion.
// ---------------------------------------------------------------------------------------------

/// `m(9_46)` (9 crossings) — the simplest of the [HS24] knots. Also `P(-3, 3, -3)`.
fn inv_9_46() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]
    ])
}

/// `15n_103488` (15 crossings).
fn inv_15n_103488() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [1,11,2,10],[2,20,3,19],[5,17,6,16],[6,25,7,26],[9,22,10,23],[12,30,13,29],[14,8,15,7],[15,27,16,26],[18,4,19,3],[20,11,21,12],[21,1,22,30],[23,4,24,5],[24,18,25,17],[27,8,28,9],[28,14,29,13]
    ])
}

/// `17nh_73`, the positron knot — `J_0` of Theorem 1. A non-negative 17-crossing diagram.
fn inv_17nh_73() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [1,27,2,26],[3,13,4,12],[4,31,5,32],[8,27,9,28],[9,1,10,34],[10,18,11,17],[13,7,14,6],[14,21,15,22],[18,26,19,25],[19,2,20,3],[20,8,21,7],[23,33,24,32],[24,11,25,12],[28,16,29,15],[29,23,30,22],[30,5,31,6],[33,16,34,17]
    ])
}

/// `7_6a` — an ssi-trivial control: adjacent `Iτ`-invariant order-1 torsion is present, but the
/// pairing splits, so `(s̲, s̄) = (s, s)`.
fn inv_7_6a() -> InvLink {
    InvLink::from_symmetric_pd_code([
        [2,13,3,14],[4,11,5,12],[6,4,7,3],[8,1,9,2],[10,5,11,6],[12,10,13,9],[14,7,1,8]
    ])
}
