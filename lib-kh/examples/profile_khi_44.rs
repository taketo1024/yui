//! One-shot run of the equivariant builder on a 44-crossing strongly invertible
//! knot. Intentionally NOT in the regular bench harness — too slow.
//!
//! Env knobs: `H_RANGE=a..=b` restricts the build; `CHUNK_BOUND=n` builds in
//! ≤n-crossing chunks; `PAIR_PENALTY=f` tunes the off-axis chooser handicap;
//! `SELECTIVE=1` selective deloop; `PREPROCESS=0` skips the half-mirror;
//! `NODE=mincut` uses the min-cutwidth crossing order.

use std::time::Instant;

use yui_core::num::FF2;
use yui_link::InvLink;
use yui_kh::tng::builder::{SymTngBuilder, SymBuildConfig, BuildMode, NodeOrder};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
    env_logger::Builder::from_default_env()
        .format_timestamp_millis()
        .init();

    let pd: &[[u8; 4]] = &[
        [2,64,3,63],[4,62,5,61],[7,31,8,30],[8,79,9,80],[11,19,12,18],
        [12,47,13,48],[13,75,14,74],[14,35,15,36],[19,7,20,6],[20,59,21,60],
        [24,86,25,85],[26,84,27,83],[27,43,28,42],[28,67,29,68],[31,11,32,10],
        [32,55,33,56],[37,51,38,50],[38,15,39,16],[39,35,40,34],[40,75,41,76],
        [43,23,44,22],[45,61,46,60],[46,5,47,6],[51,37,52,36],[52,73,53,74],
        [53,49,54,48],[54,17,55,18],[57,81,58,80],[58,29,59,30],[62,4,63,3],
        [64,2,65,1],[65,45,66,44],[66,21,67,22],[69,77,70,76],[70,33,71,34],
        [71,17,72,16],[72,49,73,50],[77,57,78,56],[78,9,79,10],[81,69,82,68],
        [82,41,83,42],[84,26,85,25],[86,24,1,23],
    ];
    let l = InvLink::from_symmetric_pd_code(pd.iter().copied());

    let zero = FF2::default();

    let h_range = std::env::var("H_RANGE").ok().map(|s| {
        let (a, b) = s.split_once("..=").expect("H_RANGE must be `a..=b`");
        a.trim().parse::<isize>().unwrap() ..= b.trim().parse::<isize>().unwrap()
    });
    let pair_penalty_coeff = std::env::var("PAIR_PENALTY").ok().map_or(1.0, |s| s.parse::<f64>().unwrap());
    let chunk_bound = std::env::var("CHUNK_BOUND").ok().map(|s| s.parse::<usize>().unwrap());
    // ELIM=deferred → MinFill; else SELECTIVE=1 → Selective; else Greedy.
    let mode = match std::env::var("ELIM").ok().as_deref() {
        Some("deferred") | Some("d") => BuildMode::MinFill,
        _ if std::env::var("SELECTIVE").is_ok() => BuildMode::Selective,
        _ => BuildMode::Greedy,
    };
    let node = match std::env::var("NODE").ok().as_deref() {
        Some("mincut") | Some("min-cut") => NodeOrder::MinCut,
        _ => NodeOrder::LoopGreedy,
    };
    let preprocess = std::env::var("PREPROCESS").map_or(true, |s| s != "0");
    let cfg = SymBuildConfig { pair_penalty_coeff, chunk_bound, h_range: h_range.clone(), mode, node, preprocess, ..Default::default() };

    let t0 = Instant::now();
    let b = SymTngBuilder::<FF2>::from_inv_link(&l, &zero, &zero, false).with_config(cfg).run();
    let elapsed = t0.elapsed();
    let c = b.into_tng_complex();

    println!("\n=== summary ===");
    println!("knot: 44-crossing (strongly invertible target)");
    println!("h_range: {:?}, chunk_bound: {:?}", h_range, chunk_bound);
    println!("build total: {:?}, verts: {}", elapsed, c.n_verts());
}
