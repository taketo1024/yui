//! One-shot run of `KhComplex::new` on the 44-crossing link from
//! Bar-Natan's "Tweaking JavaKh" page (3-twisted double of 11n12).
//! Intentionally NOT in the regular bench harness — each iteration is
//! minutes long, and criterion's statistical model wastes hours on this scale.
//!
//! Usage:
//! ```
//! RUST_LOG=info cargo run --release --example profile_kh_br44
//! ```
//!
//! Profile with samply:
//! ```
//! cargo build --profile=profiling --example profile_kh_br44
//! samply record target/profiling/examples/profile_kh_br44
//! ```

use std::time::Instant;

use num_bigint::BigInt;
use yui_link::Braid;
use yui_kh::kh::KhComplex;
use yui_kh::tng::builder::{BuildConfig, BuildMode};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
    env_logger::Builder::from_default_env()
        .format_timestamp_millis()
        .init();

    // 8-strand braid word from
    // <https://katlas.org/wiki/Tweaking_JavaKh>
    let word: Vec<i32> = vec![
        2, 1, 3, 2, 2, 1, 3, 2, -4, -5, -3, -4, 2, 1, 3, 2,
        -4, -5, -3, -4, -4, -5, -3, -4, 6, 5, 7, 6, 4, 3, 5, 4,
        -2, -3, -1, -2, 4, 3, 5, 4, 6, 5, 7, 6,
    ];
    let braid = Braid::new(8, Braid::from_iter(word.iter().copied()).elements().to_vec());
    let l = braid.closure();

    let cfg = BuildConfig { mode: BuildMode::Greedy, ..Default::default() };

    // `BigInt` because even `i128` overflows mid-run on this knot (the
    // earlier i128 attempt panicked at step 42/45 in apply_bilin's `r * s`).
    let t0 = Instant::now();
    let c = KhComplex::<BigInt>::new_with_config(&l, &BigInt::from(0), &BigInt::from(0), false, cfg);
    let elapsed = t0.elapsed();

    println!("\n=== summary ===");
    println!("knot: 44-crossing (3-twisted double of 11n12)");
    println!("KhComplex::new total: {:?}", elapsed);
    println!("h-range: {:?}", c.h_range());
}
