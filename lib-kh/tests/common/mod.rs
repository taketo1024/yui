//! Helpers shared by the paper test-crates (`invkh1`, `invkh2`, `invkh3`).
//!
//! Knot data lives at the bottom of each test file, next to the claims that use it — only the
//! machinery is shared here.
#![allow(dead_code)] // not every test crate uses every helper

use std::ops::RangeInclusive;

use yui_kh::ss::{ssi_invariant_with, SsVersion};
use yui_kh::tng::builder::{CutOption, SymBuildConfig};
use yui_link::{InvLink, Link};

/// RUST_LOG-controlled logging to stdout (tests initialize no logger by default).
pub fn init_logger() {
    use env_logger::{Builder, Target};
    let _ = Builder::from_default_env().target(Target::Stdout).try_init();
}

/// The build window `ssi_invariant` uses by default: the complex starts at `-n_neg`, and the two
/// equivariant Lee classes sit at `h = 0` and `h = 1`. Without it the whole h-range gets built,
/// which is the difference between seconds and hours on the larger diagrams.
pub fn h_range(l: &InvLink) -> RangeInclusive<isize> {
    -(l.n_signed_crossings().1 as isize) ..= 1
}

/// `ssi` with the default window and no `expected` hint — the common case.
pub fn ssi(l: &InvLink) -> (i32, i32) {
    ssi_with(l, SymBuildConfig::default(), None)
}

/// `ssi` with an explicit build config and `expected` hint. A config that leaves `h_range` unset
/// gets [`h_range`]; passing one fixes the window instead.
pub fn ssi_with(l: &InvLink, config: SymBuildConfig, expected: Option<isize>) -> (i32, i32) {
    let config = SymBuildConfig {
        h_range: config.h_range.clone().or_else(|| Some(h_range(l))),
        ..config
    };
    ssi_invariant_with(l, false, config, expected, SsVersion::V2)
}

/// [`ssi_with`] on the positive-writhe representative (the cheaper build, since the window is
/// `-n_neg ..= 1`), converted back by `(s̲(K*), s̄(K*)) = (−s̄(K), −s̲(K))` (InvKh1, Prop 1.3).
/// `expected` is the guessed `s̄` on *that* side; a wrong guess costs time, not correctness.
pub fn ssi_positive(l: &InvLink, config: SymBuildConfig, expected: Option<isize>) -> (i32, i32) {
    if l.writhe() >= 0 {
        ssi_with(l, config, expected)
    } else {
        let (s_lo, s_hi) = ssi_with(&l.mirror(), config, expected);
        (-s_hi, -s_lo)
    }
}

/// Which representative to compute on.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Side {
    /// Compute `l` exactly as given — use where the un-mirrored computation is the point.
    AsGiven,
    /// Compute the positive-writhe representative and convert back by Prop 1.3. Usually much
    /// cheaper, since the build window is `-n_neg ..= 1`.
    Positive,
}

/// Assert the published pair. `expected` is the value for `l` itself, exactly as printed in the
/// paper; the solver hint is derived from it — `s̄(K)` when `K` is computed directly, `-s̲(K)`
/// when the mirror is.
pub fn assert_ssi(l: &InvLink, config: SymBuildConfig, side: Side, expected: (i32, i32)) {
    let got = match side {
        Side::AsGiven => ssi_with(l, config, Some(expected.1 as isize)),
        Side::Positive => {
            let hint = if l.writhe() >= 0 { expected.1 } else { -expected.0 };
            ssi_positive(l, config, Some(hint as isize))
        }
    };
    assert_eq!(got, expected);
}

/// A config that only sets the chunking.
pub fn cut_config(cut: CutOption) -> SymBuildConfig {
    SymBuildConfig { cut, ..Default::default() }
}

// `2K` built as an ordinary connected sum: `K # K` is already the flip diagram, and its axis is
// the band, so `conn_sum` leaves the base point exactly there and the involution comes back by
// traversal. The band's other side is the second on-axis edge, where `conn_sum` then splices.
pub fn flip_sum(l: &Link) -> InvLink {
    let sum = l.conn_sum(l);
    let base = sum.base_pt().expect("conn_sum bases K # K on the band");
    InvLink::from_symmetric_pd_code(sum.reindexed(base, 1).pd_code())
}

// `K̃ = 2K # (−K) # (−K)`, the construction of Proposition 1.5.
pub fn flip_construction(j: InvLink) -> InvLink {
    flip_sum(j.inner()).conn_sum(&j.mirror()).conn_sum(&j.mirror())
}
