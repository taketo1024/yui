// KhI ssi on the SI-cork dataset; the ssi values established by the ssi-corks experiments serve as
// oracles. `ssi_invariant` is the streamlined cone vector pipeline; `ssi_invariant_v1` is the
// matrix cross-check.

use yui_core::poly::Poly;
use yui_core::num::FF2;
use yui_kh::khi::{ssi_invariant_v1, ssi_invariant};
use yui_kh::tng::builder::{CutOption, SymBuildConfig};

mod common;
use common::*;

type P = Poly<'H', FF2>;

#[test]
fn ssi_9_46() {
    let l = inv(k9_46());
    let c = P::variable();
    let ssi = ssi_invariant(&l, &c, false, SymBuildConfig::default());
    assert_eq!(ssi, (0, 2), "oracle from the ssi-corks experiments");
    assert_eq!(ssi, ssi_invariant_v1(&l, &c, false), "matrix-cone cross-check");
}

#[test]
#[ignore = "slow: 30-crossing interlock"]
fn ssi_interlock_9_46() {
    let (pd, cut) = interlock_9_46();
    let l = inv(pd);
    let config = SymBuildConfig { cut: CutOption::Manual(cut), ..Default::default() };
    assert_eq!(ssi_invariant(&l, &P::variable(), false, config), (0, 4), "oracle from the ssi-corks experiments");
}

#[test]
#[ignore = "heavy: 48-crossing interlock"]
fn ssi_interlock_17nh() {
    let (pd, cut) = interlock_17nh();
    let l = inv(pd);
    let config = SymBuildConfig { cut: CutOption::Manual(cut), ..Default::default() };
    assert_eq!(ssi_invariant(&l, &P::variable(), false, config), (0, 4), "oracle from the ssi-corks experiments");
}
