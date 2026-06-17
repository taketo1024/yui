//! Quick greedy-vs-selective deloop comparison on small knots.
//! Run: cargo run --release --example compare_selective -p yui-kh

use std::time::Instant;
use yui_link::Link;
use yui_kh::kh::KhComplex;
use yui_kh::tng::builder::{BuildConfig, BuildMode};

fn time_us(l: &Link, mode: BuildMode, reps: usize) -> f64 {
    let cfg = BuildConfig { mode, ..Default::default() };
    let build = || { let _ = KhComplex::<i32>::new_with_config(l, &0, &0, false, cfg.clone()); };
    for _ in 0..reps / 5 { build(); } // warm up (avoid cold-cache bias from mode order)
    let t0 = Instant::now();
    for _ in 0..reps { build(); }
    t0.elapsed().as_secs_f64() / reps as f64 * 1e6
}

fn main() {
    println!("{:12} {:>12} {:>12} {:>8}", "knot", "greedy(µs)", "selective(µs)", "ratio");
    for name in ["3_1", "5_1", "6_3", "8_19", "14n_19265"] {
        let l = Link::test_data(name);
        let reps = if name == "14n_19265" { 30 } else { 2000 };
        let g = time_us(&l, BuildMode::Greedy, reps);
        let s = time_us(&l, BuildMode::Selective, reps);
        println!("{name:12} {g:12.1} {s:12.1} {:7.2}x", s / g);
    }
}
