// `ssi` computations on larger diagrams, with the values established by the ssi-corks
// experiments as oracles; the fast per-knot correctness tests live in `src/khi/ssi.rs`.
// The heavy ones are `#[ignore]`d — run with `--ignored --nocapture` (and
// `--features big-link` for Wh(P5)).

use yui_kh::ssi::{ssi_invariant_with, SsiVersion};
use yui_kh::tng::builder::{Strategy, CutOption, SymBuildConfig};
use yui_link::InvLink;

mod common;
use common::*;

// RUST_LOG-controlled logging to stdout (tests initialize no logger by default).
fn init_logger() {
    use env_logger::{Builder, Target};
    let _ = Builder::from_default_env().target(Target::Stdout).try_init();
}

// The production Wh⁺-pretzel runs, reproducible on any machine (the in-memory construction
// fixes the crossing order; a PD re-import would reorder under MinCut).
fn wh_pretzel_config(cut_at: usize) -> SymBuildConfig {
    SymBuildConfig {
        strategy: Strategy::MinFill,
        cut: CutOption::At(vec![cut_at]),
        max_elim_cost: Some(1 << 16),
        no_full_deloop: true,
        ..Default::default()
    }
}

#[test]
fn k9_46() {
    let l = inv(common::k9_46());
    for ver in [SsiVersion::V1, SsiVersion::V2] {
        let ssi = ssi_invariant_with(&l, false, SymBuildConfig::default(), None, ver);
        assert_eq!(ssi, (0, 2), "{ver:?}");
    }
}

#[test]
#[ignore = "slow: 30-crossing interlock"]
fn interlock_9_46() {
    let (pd, _) = common::interlock_9_46();
    let l = inv(pd);
    let config = SymBuildConfig { cut: CutOption::Auto(2), ..Default::default() };
    assert_eq!(ssi_invariant_with(&l, false, config, None, SsiVersion::V2), (0, 4));
}

#[test]
#[ignore = "heavy: 48-crossing interlock"]
fn interlock_17nh() {
    let (pd, _) = common::interlock_17nh();
    let l = inv(pd);
    let config = SymBuildConfig { cut: CutOption::Auto(2), ..Default::default() };
    assert_eq!(ssi_invariant_with(&l, false, config, None, SsiVersion::V2), (0, 4));
}

#[test]
#[ignore = "heavy: 44-crossing Whitehead double (~1 min)"]
fn wh_pretzel_3() {
    init_logger();
    let k = InvLink::sym_pretzel(-3, 3, -3);
    let w = InvLink::whitehead_double(&k, true, 0);
    let ssi = ssi_invariant_with(&w, false, wh_pretzel_config(17), Some(2), SsiVersion::V2);
    assert_eq!(ssi, (0, 2));
}

#[cfg(feature = "big-link")]
#[test]
#[ignore = "heavy: 72-crossing Whitehead double"]
fn wh_pretzel_5() {
    init_logger();
    let k = InvLink::sym_pretzel(-5, 5, -5);
    let w = InvLink::whitehead_double(&k, true, 0);
    let ssi = ssi_invariant_with(&w, false, wh_pretzel_config(20), Some(2), SsiVersion::V2);
    println!("ssi(Wh+(P(-5,5,-5))) = {ssi:?}");
}

#[test]
#[ignore = "slow: 15-crossing SI knot, both pipelines"]
fn k15n_103488() {
    let l = InvLink::from_symmetric_pd_code(
        [[1,11,2,10],[2,20,3,19],[5,17,6,16],[6,25,7,26],[9,22,10,23],[12,30,13,29],[14,8,15,7],[15,27,16,26],[18,4,19,3],[20,11,21,12],[21,1,22,30],[23,4,24,5],[24,18,25,17],[27,8,28,9],[28,14,29,13]]
    );

    for ver in [SsiVersion::V1, SsiVersion::V2] {
        let ssi = ssi_invariant_with(&l, false, SymBuildConfig::default(), None, ver);
        assert_eq!(ssi, (0, 2), "{ver:?}");
    }
}

#[test]
#[ignore = "slow: 17-crossing SI knot via the matrix pipeline"]
fn k17nh_73() {
    let l = InvLink::from_symmetric_pd_code(
        [[1,27,2,26],[19,2,20,3],[3,13,4,12],[4,31,5,32],[30,5,31,6],[13,7,14,6],[8,27,9,28],[9,1,10,34],[10,18,11,17],[24,11,25,12],[14,21,15,22],[28,16,29,15],[33,16,34,17],[18,26,19,25],[20,8,21,7],[29,23,30,22],[23,33,24,32]]
    );

    let ssi = ssi_invariant_with(&l, false, SymBuildConfig::default(), None, SsiVersion::V1);

    assert_eq!(ssi, (0, 2));
}
