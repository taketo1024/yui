//! One-shot run of `KhIComplex::new` on the 60-crossing strongly invertible
//! target knot. Intentionally NOT in the regular bench harness — too slow.
//!
//! Usage:
//! ```
//! RUST_LOG=info cargo run --release --example profile_khi_60
//! ```
//!
//! Profile with samply:
//! ```
//! cargo build --profile=profiling --example profile_khi_60
//! samply record target/profiling/examples/profile_khi_60
//! ```

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
        [1,17,2,16],[6,23,7,24],[7,115,8,114],[8,103,9,104],[9,3,10,2],
        [10,96,11,95],[14,47,15,48],[15,91,16,90],[18,113,19,114],[19,25,20,24],
        [20,5,21,6],[21,101,22,100],[26,112,27,111],[28,33,29,34],[29,73,30,72],
        [31,107,32,106],[34,72,35,71],[36,63,37,64],[37,43,38,42],[38,55,39,56],
        [39,83,40,82],[44,77,45,78],[45,61,46,60],[48,13,49,14],[49,93,50,92],
        [50,88,51,87],[56,41,57,42],[57,65,58,64],[58,85,59,86],[59,53,60,52],
        [61,77,62,76],[66,83,67,84],[67,55,68,54],[68,43,69,44],[69,63,70,62],
        [70,36,71,35],[74,107,75,108],[75,31,76,30],[78,53,79,54],[79,85,80,84],
        [80,65,81,66],[81,41,82,40],[86,52,87,51],[88,93,89,94],[89,13,90,12],
        [91,47,92,46],[94,12,95,11],[96,3,97,4],[97,103,98,102],[98,115,99,116],
        [99,23,100,22],[104,17,105,18],[105,1,106,120],[108,73,109,74],[109,33,110,32],
        [110,28,111,27],[116,101,117,102],[117,5,118,4],[118,25,119,26],[119,113,120,112],
    ];
    let l = InvLink::from_symmetric_pd_code(pd.iter().copied());

    let zero = FF2::default();

    // `H_RANGE=a..=b` restricts the build; `PREPROCESS_BOUND=n` caps the symmetric chunk.
    let h_range = std::env::var("H_RANGE").ok().map(|s| {
        let (a, b) = s.split_once("..=").expect("H_RANGE must be `a..=b`");
        a.trim().parse::<isize>().unwrap() ..= b.trim().parse::<isize>().unwrap()
    });
    let preprocess_bound = std::env::var("PREPROCESS_BOUND").ok().map(|s| s.parse::<usize>().unwrap());
    let pair_penalty_coeff = std::env::var("PAIR_PENALTY").ok().map_or(1.0, |s| s.parse::<f64>().unwrap());
    let chunk_bound = std::env::var("CHUNK_BOUND").ok().map(|s| s.parse::<usize>().unwrap());

    let cfg = SymBuildConfig { preprocess_bound, pair_penalty_coeff, chunk_bound, h_range: h_range.clone(), ..Default::default() };

    let t0 = Instant::now();
    let b = SymTngBuilder::<FF2>::from_inv_link(&l, &zero, &zero, false).with_config(cfg).run();
    let elapsed = t0.elapsed();
    let c = b.into_tng_complex();

    println!("\n=== summary ===");
    println!("knot: 60 crossings (strongly invertible target)");
    println!("h_range: {:?}, preprocess_bound: {:?}, chunk_bound: {:?}", h_range, preprocess_bound, chunk_bound);
    println!("build total: {:?}, verts: {}", elapsed, c.n_verts());
}
