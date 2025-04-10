use std::hash::Hash;
use derive_more::Display;
use itertools::Either;

use crate::{Elem, ElemBase};

pub trait Gen: Elem + Hash + Ord {}

#[derive(Debug, Display, Default, Hash, PartialEq, Eq, Clone, PartialOrd, Ord)]
#[display("<{}>", _0)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(transparent))]
pub struct FreeGen<T>(pub T) where T: ElemBase;

impl<T> From<T> for FreeGen<T> 
where T: ElemBase {
    fn from(value: T) -> Self {
        Self(value)
    }
}

impl<T> Elem for FreeGen<T> 
where T: ElemBase { 
    fn math_symbol() -> String {
        let full_name = std::any::type_name::<T>();
        let name = full_name.split("::").last().unwrap_or(full_name);
        format!("Free<{}>", name)
    }
}

impl<T> Gen for FreeGen<T> 
where T: ElemBase + Hash + Ord {}

#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug, Display)]
pub struct EitherGen<X, Y>(Either<X, Y>) where X: Gen, Y: Gen;

impl<X, Y> EitherGen<X, Y> where X: Gen, Y: Gen {
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
}

impl<X, Y> From<Either<X, Y>> for EitherGen<X, Y> where X: Gen, Y: Gen {
    fn from(e: Either<X, Y>) -> Self {
        Self(e)
    }
}

impl<X, Y> From<EitherGen<X, Y>> for Either<X, Y> where X: Gen, Y: Gen {
    fn from(e: EitherGen<X, Y>) -> Self {
        e.0
    }
}

impl<X, Y> Default for EitherGen<X, Y> where X: Gen, Y: Gen {
    fn default() -> Self {
        Self(Either::Left(X::default()))
    }
}

impl<X, Y> Elem for EitherGen<X, Y> where X: Gen, Y: Gen {
    fn math_symbol() -> String {
        if X::math_symbol() == Y::math_symbol() {
            X::math_symbol()
        } else { 
            format!("E({},{})", X::math_symbol(), Y::math_symbol())
        }
    }
}

impl <X, Y> Gen for EitherGen<X, Y> where X: Gen, Y: Gen {
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lc::FreeGen; // Assuming Free is defined in the crate

    #[test]
    fn test_either_gen() {
        type T = EitherGen<FreeGen<i32>, FreeGen<String>>;
        let a = T::from_left(FreeGen(42));
        let b = T::from_right(FreeGen("hello".to_string()));

        assert_eq!(a.to_string(), "<42>");
        assert_eq!(b.to_string(), "<hello>");
    }
}