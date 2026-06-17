//! Small-knot timing sweep: greedy vs selective vs min-fill on sym (KhI) builds.
//! Selective is expected to (if ever) win here, on small non-chunked knots.
//! Run: cargo run --release --example survey_modes_small -p yui-kh

use std::time::Instant;

use yui_core::num::FF2;
use yui_link::InvLink;
use yui_kh::tng::builder::{SymTngBuilder, SymBuildConfig, BuildMode};

fn time_us(l: &InvLink, mode: BuildMode, reps: usize) -> f64 {
    let zero = FF2::default();
    let run = || {
        let cfg = SymBuildConfig { mode, preprocess: true, ..Default::default() };
        let _ = SymTngBuilder::<FF2>::from_inv_link(l, &zero, &zero, false)
            .with_config(cfg).run().into_tng_complex();
    };
    for _ in 0..(reps / 5).max(1) { run(); } // warmup
    let t0 = Instant::now();
    for _ in 0..reps { run(); }
    t0.elapsed().as_secs_f64() / reps as f64 * 1e6
}

fn main() {
    let k9_46 = InvLink::from_symmetric_pd_code(
        [[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]
    );
    let knots: [(&str, InvLink, usize); 4] = [
        ("3_1", InvLink::test_data("3_1"), 3000),
        ("4_1", InvLink::test_data("4_1"), 2000),
        ("6_3", InvLink::test_data("6_3"),  500),
        ("9_46", k9_46,                      80),
    ];

    println!("{:8} {:>12} {:>12} {:>12}   {:>10} {:>10}", "knot", "greedy", "selective", "minfill", "sel/grd", "mf/grd");
    for (name, l, reps) in &knots {
        let g = time_us(l, BuildMode::Greedy,    *reps);
        let s = time_us(l, BuildMode::Selective, *reps);
        let m = time_us(l, BuildMode::MinFill,   *reps);
        println!("{name:8} {g:11.1}µ {s:11.1}µ {m:11.1}µ   {:9.2}x {:9.2}x", s / g, m / g);
    }
}
