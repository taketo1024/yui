//! Key types for [`Lc`](super::Lc): the trait [`LcKey`] and the two
//! constructions [`AsKey`] (wrap an arbitrary element as a key) and
//! [`EitherKey`] (disjoint union of two key sets).

use std::hash::Hash;
use derive_more::Display;
use itertools::Either;

use crate::lc::Lc;
use crate::{Elem, ElemBase, Ring, RingOps};

/// Marker trait for types usable as keys in [`Lc`](super::Lc) — i.e.
/// elements that are hashable and totally ordered.
pub trait LcKey: Elem + Hash + Ord {}

/// Wraps an arbitrary element `T` so it can be used as an [`LcKey`].
///
/// Used to build the [free module](crate::RMod) over any `T: ElemBase + Hash + Ord`.
#[derive(Debug, Display, Default, Hash, PartialEq, Eq, Clone, PartialOrd, Ord)]
#[display("<{}>", _0)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(transparent))]
pub struct AsKey<T>(pub T) where T: ElemBase;

impl<T> From<T> for AsKey<T> 
where T: ElemBase {
    fn from(value: T) -> Self {
        Self(value)
    }
}

impl<T> Elem for AsKey<T> 
where T: ElemBase { 
    fn math_symbol() -> String {
        let full_name = std::any::type_name::<T>();
        let name = full_name.split("::").last().unwrap_or(full_name);
        format!("Free<{}>", name)
    }
}

impl<T> LcKey for AsKey<T> 
where T: ElemBase + Hash + Ord {}

/// A disjoint union `X ⊔ Y` of two key sets, used to form direct sums of
/// linear combinations, e.g. the basis for `Lc<X, R> ⊕ Lc<Y, R> ≅ Lc<EitherKey<X, Y>, R>`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct EitherKey<X, Y>(Either<X, Y>) where X: LcKey, Y: LcKey;

impl<X, Y> EitherKey<X, Y> where X: LcKey, Y: LcKey {
    pub fn from_left(x: X) -> Self {
        Self(Either::Left(x))
    }

    pub fn from_right(y: Y) -> Self {
        Self(Either::Right(y))
    }

    pub fn entity(&self) -> Either<&X, &Y> {
        match &self.0 {
            Either::Left(x) => Either::Left(x),
            Either::Right(y) => Either::Right(y),
        }
    }

    pub fn is_left(&self) -> bool {
        matches!(self.0, Either::Left(_))
    }

    pub fn is_right(&self) -> bool {
        matches!(self.0, Either::Right(_))
    }

    pub fn inner(&self) -> &Either<X, Y> {
        &self.0
    }

    pub fn into_left(self) -> X {
        let Either::Left(x) = self.0 else {
            panic!();
        };
        x
    }

    pub fn into_right(self) -> Y {
        let Either::Right(y) = self.0 else {
            panic!();
        };
        y
    }
}

impl<X, Y> From<Either<X, Y>> for EitherKey<X, Y> where X: LcKey, Y: LcKey {
    fn from(e: Either<X, Y>) -> Self {
        Self(e)
    }
}

impl<X, Y> From<EitherKey<X, Y>> for Either<X, Y> where X: LcKey, Y: LcKey {
    fn from(e: EitherKey<X, Y>) -> Self {
        e.0
    }
}

impl<X, Y> Default for EitherKey<X, Y> where X: LcKey, Y: LcKey {
    fn default() -> Self {
        Self(Either::Left(X::default()))
    }
}

impl <X, Y> std::fmt::Display for EitherKey<X, Y> where X: LcKey, Y: LcKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.0 {
            Either::Left(x)  => std::fmt::Display::fmt(x, f),
            Either::Right(y) => std::fmt::Display::fmt(y, f),
        }
    }
}

impl<X, Y> Elem for EitherKey<X, Y> where X: LcKey, Y: LcKey {
    fn math_symbol() -> String {
        if X::math_symbol() == Y::math_symbol() {
            X::math_symbol()
        } else { 
            format!("E({},{})", X::math_symbol(), Y::math_symbol())
        }
    }
}

impl <X, Y> LcKey for EitherKey<X, Y> where X: LcKey, Y: LcKey {
}

pub fn split_lr<X, Y, R>(z: &Lc<EitherKey<X, Y>, R>) -> (Lc<X, R>, Lc<Y, R>)
where X: LcKey, Y: LcKey, R: Ring, for<'x> &'x R: RingOps<R>{
    let mut x = vec![];
    let mut y = vec![];
    for (e, r) in z.iter() { 
        if e.is_left() { 
            x.push((e.clone().into_left(), r.clone()));
        } else { 
            y.push((e.clone().into_right(), r.clone()));
        }
    } 
    (Lc::from_iter(x), Lc::from_iter(y))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lc::AsKey; // Assuming Free is defined in the crate

    #[test]
    fn test_either_key() {
        type T = EitherKey<AsKey<i32>, AsKey<String>>;
        let a = T::from_left(AsKey(42));
        let b = T::from_right(AsKey("hello".to_string()));

        assert_eq!(a.to_string(), "<42>");
        assert_eq!(b.to_string(), "<hello>");
    }
}