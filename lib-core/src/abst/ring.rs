//! Ring: an additive group `+` together with an associative multiplication `·`
//! that distributes over addition.
//!
//! Multiplication is **not** assumed commutative.
//!
//! See: <https://en.wikipedia.org/wiki/Ring_(mathematics)>

use crate::num::Sign;
use crate::abst::{AddGrp, AddGrpOps, Mon, MonOps};

/// Helper trait bundling [`AddGrpOps`] and [`MonOps`] for [`Ring`].
pub trait RingOps<T = Self>:
    AddGrpOps<T> +
    MonOps<T>
{}

/// A ring: an [`AddGrp`] together with an associative multiplication `·`
/// distributing over addition `+`, with multiplicative identity [`One`](num_traits::One).
///
/// See: <https://en.wikipedia.org/wiki/Ring_(mathematics)>
pub trait Ring:
    AddGrp +
    Mon +
    RingOps
where
    for<'a> &'a Self: RingOps<Self>
{
    /// `+1` if `s` is positive, `-1` otherwise.
    fn from_sign(s: Sign) -> Self {
        if s.is_positive() {
            Self::one()
        } else {
            -Self::one()
        }
    }

    /// Multiplicative inverse, if this element is a unit; otherwise `None`.
    fn inv(&self) -> Option<Self>;

    /// `true` iff this element is a unit (has a multiplicative inverse).
    ///
    /// See: <https://en.wikipedia.org/wiki/Unit_(ring_theory)>
    fn is_unit(&self) -> bool;

    /// A unit `u` such that `self * u` is the chosen canonical associate of `self`.
    ///
    /// For example, in `ℤ` the normalizing unit of `-n` is `-1` (giving `n`);
    /// in a field every nonzero element is its own normalizer (returning `self.inv()`).
    ///
    /// See: <https://en.wikipedia.org/wiki/Associate_element>
    fn normalizing_unit(&self) -> Self;

    /// The canonical associate of `self` (i.e. `self * self.normalizing_unit()`).
    fn normalized(&self) -> Self {
        self.clone().into_normalized()
    }

    /// Like [`normalized`](Self::normalized), but consuming.
    fn into_normalized(self) -> Self {
        let u = self.normalizing_unit();
        if u.is_one() {
            self
        } else {
            self * u
        }
    }

    /// `true` iff `self` is `+1` or `-1`.
    fn is_pm_one(&self) -> bool {
        self.is_one() || (-self).is_one()
    }

    /// Heuristic cost of working with this element, used by chain-reduction
    /// pivot selection. Default: `0.0` for zero, `1.0` otherwise.
    fn c_weight(&self) -> f64 {
        if self.is_zero() {
            0.0
        } else {
            1.0
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::abst::Ring;
 
    #[test]
    fn is_pm_one() { 
        assert!(1.is_pm_one());
        assert!((-1).is_pm_one());
        assert!(!2.is_pm_one());
        assert!(!(-2).is_pm_one());
    }

    #[test]
    fn normalized() {
        assert_eq!(3.normalized(), 3);
        assert_eq!((-3).normalized(), 3);
    } 

}