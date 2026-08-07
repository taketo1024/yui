//! Additive monoid: a set with associative, commutative addition `+` and identity `0`.
//!
//! See: <https://en.wikipedia.org/wiki/Monoid>,
//! <https://en.wikipedia.org/wiki/Commutative_monoid>

use std::ops::{Add, AddAssign};
use num_traits::Zero;
use crate::abst::MathType;

/// Helper trait bundling `Add` impls so [`AddMon`] can require all four
/// reference variants (`T + T`, `T + &T`, `&T + T`, `&T + &T`) via one HRTB.
pub trait AddMonOps<T = Self>:
    Sized +
    Add<T, Output = T> +              // S + T -> T
    for<'a> Add<&'a T, Output = T>    // S + &T -> T
{}

/// An additive monoid: a set with associative, commutative addition `+` and identity [`Zero`].
///
/// See: <https://en.wikipedia.org/wiki/Monoid>,
/// <https://en.wikipedia.org/wiki/Commutative_monoid>
pub trait AddMon:
    MathType +
    Zero +
    AddMonOps +                       // T + T -> T, T + &T -> T
    AddAssign +                       // T += T
    for<'a> AddAssign<&'a Self>       // T += &T
where
    for<'a> &'a Self: AddMonOps<Self> // &T + T -> T, &T + &T -> T
{
    /// Sum an iterator into `Self`, folding with `+=`.
    ///
    /// `A` is any type for which `Self: AddAssign<A>` — typically `Self` or `&Self`.
    fn sum<A, I>(itr: I) -> Self
    where
        Self: AddAssign<A>,
        I: IntoIterator<Item = A>
    {
        itr.into_iter().fold(Self::zero(), |mut res, a| {
            res += a;
            res
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sum() {
        let a = i64::sum([4,5,6]);
        assert_eq!(a, 15);
    }
}