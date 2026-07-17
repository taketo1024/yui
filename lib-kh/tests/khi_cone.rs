// KhI ssi on the SI-cork dataset; the ssi values established by the ssi-corks experiments serve as
// oracles. `ssi_invariant` is the streamlined cone vector pipeline (V2); `SsiVersion::V1` is the
// matrix cross-check.

use yui_core::poly::Poly;
use yui_core::num::FF2;
use yui_kh::khi::{ssi_invariant, ssi_invariant_ver, SsiVersion};
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
    assert_eq!(ssi, ssi_invariant_ver(&l, &c, false, SymBuildConfig::default(), SsiVersion::V1), "matrix-cone cross-check");
}

#[test]
#[ignore = "slow: 30-crossing interlock"]
fn ssi_interlock_9_46() {
    let (pd, _) = interlock_9_46();
    let l = inv(pd);
    let config = SymBuildConfig { cut: CutOption::Auto(2), ..Default::default() };
    assert_eq!(ssi_invariant(&l, &P::variable(), false, config), (0, 4), "oracle from the ssi-corks experiments");
}

#[test]
#[ignore = "heavy: 48-crossing interlock"]
fn ssi_interlock_17nh() {
    let (pd, _) = interlock_17nh();
    let l = inv(pd);
    let config = SymBuildConfig { cut: CutOption::Auto(2), ..Default::default() };
    assert_eq!(ssi_invariant(&l, &P::variable(), false, config), (0, 4), "oracle from the ssi-corks experiments");
}
