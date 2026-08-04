//! Additive group: an [`AddMon`] in which every element has an additive inverse `-x`.
//!
//! See: <https://en.wikipedia.org/wiki/Group_(mathematics)>,
//! <https://en.wikipedia.org/wiki/Additive_group>

use std::ops::{Neg, Sub, SubAssign};
use crate::{AddMon, AddMonOps};

/// Helper trait extending [`AddMonOps`] with `Neg` and `Sub` reference variants,
/// so [`AddGrp`] can require them via one HRTB.
pub trait AddGrpOps<T = Self>:
    AddMonOps<T> +
    Neg<Output = T> +
    Sub<T, Output = T> +
    for<'a> Sub<&'a T, Output = T>
{}

/// An additive group: an [`AddMon`] with negation `-x` and subtraction `x - y`,
/// satisfying `x + (-x) = 0`.
///
/// See: <https://en.wikipedia.org/wiki/Group_(mathematics)>,
/// <https://en.wikipedia.org/wiki/Additive_group>
pub trait AddGrp:
    AddMon +
    AddGrpOps +
    SubAssign +
    for<'a> SubAssign<&'a Self>
where
    for<'a> &'a Self: AddGrpOps<Self>
{}