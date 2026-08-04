//! Module over a [`Ring`]: an [`AddGrp`] equipped with scalar multiplication `r · m` by `R`.
//!
//! See: <https://en.wikipedia.org/wiki/Module_(mathematics)>

use std::ops::{Mul, MulAssign};
use crate::{AddGrp, AddGrpOps, Ring, RingOps};

/// Helper trait bundling [`AddGrpOps`] with scalar multiplication by `R`
/// (both `T * R` and `T * &R`) so [`RMod`] can require them via one HRTB.
pub trait RModOps<R, T>:
    AddGrpOps<T> +
    Mul<R, Output = T> +
    for<'a> Mul<&'a R, Output = T>
where
    R: Ring, for<'x> &'x R: RingOps<R>
{}

/// A (left) module over a ring `R = Self::R`: an [`AddGrp`] with a distributive,
/// associative scalar action `r · m` (for `r ∈ R`, `m ∈ Self`).
///
/// See: <https://en.wikipedia.org/wiki/Module_(mathematics)>
pub trait RMod:
    AddGrp +
    RModOps<Self::R, Self> +
    MulAssign<Self::R> +
    for<'a> MulAssign<&'a Self::R>
where
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    for<'a> &'a Self: RModOps<Self::R, Self>,
{
    /// The scalar ring.
    type R;
}

