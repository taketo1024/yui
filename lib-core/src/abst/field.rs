//! Field: a commutative ring in which every nonzero element has a multiplicative inverse `x⁻¹`.
//!
//! See: <https://en.wikipedia.org/wiki/Field_(mathematics)>

use crate::{EucRing, EucRingOps};

/// Helper trait, currently identical to [`EucRingOps`], kept for symmetry
/// with the rest of the algebraic hierarchy.
pub trait FieldOps<T = Self>:
    EucRingOps<T>
{}

/// A field: a commutative ring where every `x ≠ 0` has a multiplicative inverse `x⁻¹`.
///
/// Every field is in particular a [`EucRing`] (with trivial remainder).
///
/// See: <https://en.wikipedia.org/wiki/Field_(mathematics)>
pub trait Field:
    EucRing +
    FieldOps
where
    for<'a> &'a Self: FieldOps<Self>
{}