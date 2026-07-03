//! Sakuma-style equivariant connected sums: for each 2-bridge knot `K` with two strong
//! inversions `τ₁ (a), τ₂ (b)`, form `J = (K, τ₁) # (−K, −τ₂)` — slice as a plain knot,
//! possibly not equivariantly slice — and compute `ssi(J)`.
//!
//! note: 2-bridge knots are invertible, so `−K` (the concordance inverse) is just the mirror.
//! Controls `(K, τ) # (−K, −τ)` are equivariantly slice and must give ssi = (0, 0).

use std::collections::BTreeMap;
use std::fs;

use log::info;
use yui_core::poly::Poly;
use yui_core::num::FF2;
use yui_link::{Edge, InvLink, Link};
use yui_link::misc::jones_polynomial;
use yui_kh::khi::ssi_invariants_via_cone;
use yui_kh::tng::builder::{BuildMode, CutOption, SymBuildConfig};

type P = Poly<'H', FF2>;
type Pd = Vec<[Edge; 4]>;

fn ssi(l: &InvLink) -> (i32, i32) {
    ssi_invariants_via_cone(l, &P::variable(), false, SymBuildConfig::default())
}

// heavy config for the ~70-crossing Whitehead doubles, mirroring the Wh⁺(9_46) run.
fn ssi_heavy(l: &InvLink) -> (i32, i32) {
    // doubly-truncated window 0..=1: the full bottom..=1 window is memory-infeasible at this size.
    let config = SymBuildConfig {
        mode: BuildMode::MinFill,
        cut: CutOption::Auto(2),
        h_range: Some(0..=1),
        ..Default::default()
    };
    ssi_invariants_via_cone(l, &P::variable(), false, config)
}

fn load_knots() -> BTreeMap<String, Pd> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../knots.json");
    let json = fs::read_to_string(path).expect("knots.json next to the runner crate");
    serde_json::from_str(&json).unwrap()
}

// `sakuma-si wh <base> [neg]`: ssi of the (positive- or negative-clasp) Whitehead double
// of J = (K_a, τ₁) # (−K_b, −τ₂).
fn run_wh(base: &str, positive: bool) {
    let knots = load_knots();
    let ka = InvLink::from_symmetric_pd_code(knots[&format!("{base}a")].clone());
    let kb = InvLink::from_symmetric_pd_code(knots[&format!("{base}b")].clone());

    let j = ka.conn_sum(&kb.mirror());
    let wh = j.whitehead_double(positive, 0);

    info!("Wh{}({base}a # -{base}b): {} crossings, writhe {}", if positive { "+" } else { "-" }, wh.n_nodes(), wh.writhe());

    let ssi = ssi_heavy(&wh);
    println!("{base}, Wh{}, {:?}", if positive { "+" } else { "-" }, ssi);
}

// `sakuma-si pd <base> [neg]`: print the PD code of Wh±(K_a # -K_b), verifying it round-trips
// (writhe + Jones) so it can be fed to `ykh khi "<pd>" ...` on another machine.
fn run_pd(base: &str, positive: bool) {
    let knots = load_knots();
    let ka = InvLink::from_symmetric_pd_code(knots[&format!("{base}a")].clone());
    let kb = InvLink::from_symmetric_pd_code(knots[&format!("{base}b")].clone());

    let j = ka.conn_sum(&kb.mirror());
    let wh = j.whitehead_double(positive, 0);
    let inner = wh.inner();

    let pd: Vec<[Edge; 4]> = inner.nodes().map(|x| *x.edges()).collect();

    let reloaded = Link::from_pd_code(pd.clone());
    assert_eq!(reloaded.writhe(), inner.writhe(), "writhe mismatch on reload");
    assert_eq!(jones_polynomial(&reloaded), jones_polynomial(inner), "Jones mismatch on reload");

    let body = pd.iter().map(|e| format!("[{},{},{},{}]", e[0], e[1], e[2], e[3])).collect::<Vec<_>>().join(",");
    eprintln!("Wh{}({base}a # -{base}b): {} crossings, round-trip OK (writhe/Jones)", if positive { "+" } else { "-" }, pd.len());
    println!("[{body}]");
}

fn main() {
    env_logger::init();

    let args: Vec<String> = std::env::args().collect();
    if args.len() >= 3 && args[1] == "wh" {
        let positive = args.get(3).map_or(true, |s| s != "neg");
        run_wh(&args[2], positive);
        return;
    }
    if args.len() >= 3 && args[1] == "pd" {
        let positive = args.get(3).map_or(true, |s| s != "neg");
        run_pd(&args[2], positive);
        return;
    }

    let knots = load_knots();

    let mut pairs: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for name in knots.keys() {
        if name.ends_with('a') || name.ends_with('b') {
            pairs.entry(name[..name.len() - 1].to_string()).or_default().push(name.clone());
        }
    }

    println!("knot, ssi(a), ssi(b), ssi(a # -b)");

    for (base, variants) in pairs.iter().filter(|(_, v)| v.len() == 2) {
        let ka = InvLink::from_symmetric_pd_code(knots[&variants[0]].clone());
        let kb = InvLink::from_symmetric_pd_code(knots[&variants[1]].clone());

        info!("--- {base}: {} x {} crossings ---", ka.n_nodes(), kb.n_nodes());

        let ssi_a = ssi(&ka);
        let ssi_b = ssi(&kb);
        let sum = ssi(&ka.conn_sum(&kb.mirror()));

        println!("{base}, {:?}, {:?}, {:?}", ssi_a, ssi_b, sum);
    }
}
