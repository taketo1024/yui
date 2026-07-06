//! `Wh⁺((8_21, τ₁) # (−8_21, −τ₂))` — the Sakuma-style equivariant connected sum of the two
//! strong inversions of 8_21 (the unique ≤ 9-crossing pair with non-trivial `ssi(a # −b) = (−2, 0)`),
//! Whitehead-doubled: a candidate for a strongly-invertible knot with exotic slice disks.
//!
//! The double has 72 crossings, so this requires `--features big-link`, and the run is heavy
//! (OOMed at ~105 GB on a 128 GB machine during the chunk build — intended for large-memory hosts):
//!
//! ```sh
//! RUST_LOG=debug cargo test -p yui-kh --release --features big-link --test wh_sakuma \
//!     -- --ignored --nocapture 2> wh_8_21.log
//! ```
#![cfg(feature = "big-link")]

use yui_core::poly::Poly;
use yui_core::num::FF2;
use yui_link::InvLink;
use yui_kh::khi::ssi_invariants_via_cone;
use yui_kh::tng::builder::{BuildMode, CutOption, SymBuildConfig};

type P = Poly<'H', FF2>;

// transvergent PDs of the two strong inversions of 8_21 (Sakuma's table).
fn k8_21a() -> InvLink {
    InvLink::from_symmetric_pd_code(
        [[3,8,4,9],[6,12,7,11],[7,2,8,3],[9,14,10,15],[12,1,13,2],[13,5,14,4],[15,10,16,11],[16,5,1,6]]
    )
}

fn k8_21b() -> InvLink {
    InvLink::from_symmetric_pd_code(
        [[1,6,2,7],[3,16,4,17],[7,12,8,13],[10,5,11,6],[11,3,12,2],[13,18,14,1],[14,9,15,10],[15,4,16,5],[17,9,18,8]]
    )
}

#[test]
#[ignore = "heavy: 72 crossings, needs a large-memory machine"]
fn ssi_wh_8_21_pos() {
    let _ = env_logger::try_init(); // RUST_LOG controls the level; logs go to stderr (use 2>&1 for stdout)

    let j = k8_21a().conn_sum(&k8_21b().mirror());
    let wh = j.whitehead_double(true, 0);

    println!("Wh+(8_21a # -8_21b): {} crossings, writhe {}", wh.n_nodes(), wh.writhe());

    // 72 crossings: prune the discarded negative tail with the tight ssi window (`h_range = 0..=1`
    // → cone `-2..=2`), emit the cone directly, and cap the heavy eliminations (survivors defer to
    // the matrix reduction). The dense central-degree merge is still the wall — watch `merge_edges`.
    let config = SymBuildConfig {
        mode: BuildMode::MinFill,
        cut: CutOption::Auto(2),
        h_range: Some(0 ..= 1),
        cone_direct: true,
        elim_max_cost: Some(4096),
        ..Default::default()
    };
    let ssi = ssi_invariants_via_cone(&wh, &P::variable(), false, config);

    println!("ssi(Wh+(8_21a # -8_21b)) = {:?}", ssi);
}
