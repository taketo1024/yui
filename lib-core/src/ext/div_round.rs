//! [`DivRound`]: integer division rounded to the nearest, halves away from zero.

/// Division rounded to the nearest integer.
pub trait DivRound {
    fn div_round(&self, rhs: &Self) -> Self;
}