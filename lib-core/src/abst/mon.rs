//! Multiplicative monoid: a set with associative multiplication `·` and identity `1`.
//!
//! Unlike [`AddMon`](crate::AddMon), commutativity is **not** assumed here.
//!
//! See: <https://en.wikipedia.org/wiki/Monoid>

use std::ops::{Mul, MulAssign};
use num_traits::One;
use crate::abst::MathType;

/// Helper trait bundling `Mul` impls so [`Mon`] can require all four
/// reference variants (`T * T`, `T * &T`, `&T * T`, `&T * &T`) via one HRTB.
pub trait MonOps<T = Self>:
    Sized +
    Mul<T, Output = T> +
    for<'a> Mul<&'a T, Output = T>
{}

/// A multiplicative monoid: a set with associative multiplication `·` and identity [`One`].
///
/// Not assumed commutative.
///
/// See: <https://en.wikipedia.org/wiki/Monoid>
pub trait Mon:
    MathType +
    One +
    MonOps +
    MulAssign +
    for<'a> MulAssign<&'a Self>
where
    for<'a> &'a Self: MonOps<Self>
{
    /// Multiply an iterator of factors into `Self`, folding with `*=`.
    ///
    /// `A` is any type for which `Self: MulAssign<A>` — typically `Self` or `&Self`.
    fn product<A, I>(itr: I) -> Self
    where
        Self: MulAssign<A>,
        I: IntoIterator<Item = A>
    {
        itr.into_iter().fold(Self::one(), |mut res, a| {
            res *= a;
            res
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn product() {
        let a = i64::product([4,5,6]);
        assert_eq!(a, 120);
    }
}