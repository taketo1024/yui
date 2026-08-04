use itertools::Either;
use yui_core::lc::EitherKey;

use crate::kh::KhGen;

/// Generator type for the involutive Khovanov chain complex:
/// the disjoint union `B ⊔ Q` of two copies of `KhGen`. The `B` side maps to
/// the `Left` variant and the `Q` side to `Right`, matching the convention
/// used by `ChainMap::cone`.
pub type KhIGen = EitherKey<KhGen, KhGen>;

/// `KhIGen`-specific helpers not provided by [`EitherKey`].
pub trait KhIGenExt {
    fn rel_h_deg(&self) -> isize;
    fn rel_q_deg(&self) -> isize;
}

impl KhIGenExt for KhIGen {
    fn rel_h_deg(&self) -> isize {
        match self.inner() {
            Either::Left(x)  => x.rel_h_deg(),
            Either::Right(x) => x.rel_h_deg() + 1,
        }
    }

    fn rel_q_deg(&self) -> isize {
        match self.inner() {
            Either::Left(x) | Either::Right(x) => x.rel_q_deg(),
        }
    }
}
