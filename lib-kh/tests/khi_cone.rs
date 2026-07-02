// KhI via the cobordism-level cone (`ConeBuilder` + `from_cone`) on the SI-cork dataset;
// the ssi values established by the ssi-corks experiments serve as oracles.

use itertools::Itertools;
use num_traits::Zero;
use yui_core::poly::Poly;
use yui_core::num::FF2;
use yui_link::InvLink;
use yui_kh::khi::{KhIHomology, ssi_invariants};
use yui_kh::tng::builder::{CutOption, SymBuildConfig};
use yui_kh::util::calc::div_vec;

mod common;
use common::*;

type P = Poly<'H', FF2>;

// ssi (s̲, s̄) from the canon classes via H-divisibility — cf. `cone_canon_ssi_matches_matrix`.
fn ssi_with(l: &InvLink, config: SymBuildConfig, cone: bool) -> (i32, i32) {
    let (c, t) = (P::variable(), P::zero());
    let config = SymBuildConfig { h_range: Some(isize::MIN + 1 ..= 1), ..config };
    let kh = if cone {
        KhIHomology::from_cone(l, &c, &t, false, config)
    } else {
        KhIHomology::new_with_config(l, &c, &t, false, config)
    };

    let zs = kh.canon_cycles();
    assert_eq!(zs.len(), 4);

    let ds = zs.iter().map(|z| {
        let h = kh.h_deg_of_chain(z);
        div_vec(&kh[h].vectorize_euc(z).subvec(0..2), &c).expect("invalid divisibility")
    }).collect_vec();

    let (w, r) = (l.writhe(), l.seifert_circles().len() as i32);
    (2 * ds[0] + w - r + 1, 2 * ds[2] + w - r + 1)
}

fn ssi_via_cone(l: &InvLink, config: SymBuildConfig) -> (i32, i32) {
    ssi_with(l, config, true)
}

#[test]
fn ssi_9_46_cone() {
    let l = inv(k9_46());
    let ssi = ssi_via_cone(&l, SymBuildConfig::default());
    assert_eq!(ssi, (0, 2), "oracle from the ssi-corks experiments");
    assert_eq!(ssi, ssi_invariants(&l, &P::variable(), false), "matrix-cone cross-check");
}

#[test]
#[ignore = "slow: 30-crossing interlock"]
fn ssi_interlock_9_46_cone() {
    let (pd, cut) = interlock_9_46();
    let l = inv(pd);
    let config = SymBuildConfig { cut: CutOption::Manual(cut), ..Default::default() };
    assert_eq!(ssi_with(&l, config, true), (0, 4), "oracle from the ssi-corks experiments");
}

#[test]
#[ignore = "slow: 30-crossing interlock"]
fn ssi_interlock_9_46_matrix() {
    let (pd, cut) = interlock_9_46();
    let l = inv(pd);
    let config = SymBuildConfig { cut: CutOption::Manual(cut), ..Default::default() };
    assert_eq!(ssi_with(&l, config, false), (0, 4), "oracle from the ssi-corks experiments");
}
