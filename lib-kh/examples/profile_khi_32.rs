//! One-shot run of the equivariant builder on a 32-crossing strongly invertible
//! knot. Intentionally NOT in the regular bench harness — too slow.
//!
//! Usage:
//! ```
//! RUST_LOG=info cargo run --release --example profile_khi_32
//! ```
//!
//! Env knobs: `H_RANGE=a..=b` restricts the build; `CHUNK_BOUND=n` builds in
//! ≤n-crossing chunks; `PAIR_PENALTY=f` tunes the off-axis chooser handicap.

use std::time::Instant;

use yui_core::num::FF2;
use yui_link::InvLink;
use yui_kh::tng::builder::{SymTngBuilder, SymBuildConfig};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
    env_logger::Builder::from_default_env()
        .format_timestamp_millis()
        .init();

    let pd: &[[u8; 4]] = &[
        [48,2,49,1],[2,48,3,47],[40,3,41,4],[39,47,40,46],[9,5,10,4],
        [10,45,11,46],[8,23,9,24],[41,25,42,24],[7,59,8,58],[42,57,43,58],
        [59,7,60,6],[60,43,61,44],[22,5,23,6],[21,45,22,44],[61,57,62,56],
        [62,25,63,26],[19,27,20,26],[20,55,21,56],[18,64,19,63],[64,18,1,17],
        [11,55,12,54],[12,27,13,28],[38,53,39,54],[37,29,38,28],[29,37,30,36],
        [30,13,31,14],[52,35,53,36],[51,15,52,14],[15,51,16,50],[16,31,17,32],
        [33,33,34,32],[34,49,35,50],
    ];
    let l = InvLink::from_symmetric_pd_code(pd.iter().copied());

    let zero = FF2::default();

    let h_range = std::env::var("H_RANGE").ok().map(|s| {
        let (a, b) = s.split_once("..=").expect("H_RANGE must be `a..=b`");
        a.trim().parse::<isize>().unwrap() ..= b.trim().parse::<isize>().unwrap()
    });
    let pair_penalty_coeff = std::env::var("PAIR_PENALTY").ok().map_or(1.0, |s| s.parse::<f64>().unwrap());
    let chunk_bound = std::env::var("CHUNK_BOUND").ok().map(|s| s.parse::<usize>().unwrap());

    let cfg = SymBuildConfig { pair_penalty_coeff, chunk_bound, h_range: h_range.clone(), ..Default::default() };

    let t0 = Instant::now();
    let b = SymTngBuilder::<FF2>::from_inv_link(&l, &zero, &zero, false).with_config(cfg).run();
    let elapsed = t0.elapsed();
    let c = b.into_tng_complex();

    println!("\n=== summary ===");
    println!("knot: 32 crossings (strongly invertible target)");
    println!("h_range: {:?}, chunk_bound: {:?}", h_range, chunk_bound);
    println!("build total: {:?}, verts: {}", elapsed, c.n_verts());
}
