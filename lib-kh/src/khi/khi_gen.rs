//! Generator type for the involutive Khovanov chain complex `CKhI = Cone(1 + τ)`
//! — the disjoint union of two copies of `KhGen` (`B`-side and `Q`-side in the
//! cone), realized as [`EitherKey<KhGen, KhGen>`].
//!
//! Reference:
//! - T. Sano, "Involutive Khovanov homology and equivariant knots",
//!   Algebr. Geom. Topol. 25 (2025), 5059–5111.
//!   <https://doi.org/10.2140/agt.2025.25.5059>, <https://arxiv.org/abs/2404.08568>

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

