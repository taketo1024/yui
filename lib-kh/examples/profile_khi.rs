//! One-shot run of `KhIComplex::new` on the user-supplied 18-crossing strongly
//! invertible knot, with `env_logger` enabled so the builder's per-crossing
//! `info!` events are captured with timestamps.
//!
//! Usage:
//! ```
//! RUST_LOG=info cargo run --release --example profile_khi
//! ```
//!
//! Or, to record a sampling profile with samply:
//! ```
//! cargo build --release --example profile_khi
//! samply record target/release/examples/profile_khi
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
        [1,27,2,26],[5,16,6,17],[6,32,7,31],[10,27,11,28],[11,1,12,36],
        [13,8,14,9],[14,20,15,19],[17,4,18,5],[18,24,19,23],[21,32,22,33],
        [22,16,23,15],[25,3,26,2],[28,9,29,10],[29,24,30,25],[30,4,31,3],
        [33,20,34,21],[34,8,35,7],[35,13,36,12],
    ];
    let l = InvLink::from_symmetric_pd_code(pd.iter().copied());

    let zero = FF2::default();

    // isolate the variable that changed: the SymTngBuilder run with preprocess
    // on vs off (cone assembly afterwards is identical either way).
    let build = |preprocess: bool| {
        let cfg = SymBuildConfig { preprocess, ..Default::default() };
        let t0 = Instant::now();
        let b = SymTngBuilder::<FF2>::from_inv_link(&l, &zero, &zero, false).with_config(cfg).run();
        let dt = t0.elapsed();
        let c = b.into_tng_complex();
        (dt, c.n_verts())
    };

    let (t_off, v_off) = build(false);
    let (t_on, v_on) = build(true);

    println!("\n=== summary (18 crossings) ===");
    println!("preprocess off: {:?}  (verts {})", t_off, v_off);
    println!("preprocess on : {:?}  (verts {})", t_on, v_on);
}
