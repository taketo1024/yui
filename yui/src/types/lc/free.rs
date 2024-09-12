use std::hash::Hash;
use derive_more::Display;
use crate::Elem;
use crate::lc::Gen;

#[derive(Debug, Display, Default, Hash, PartialEq, Eq, Clone, PartialOrd, Ord)]
#[display("<{}>", _0)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(transparent))]
pub struct Free<T>(pub T) where T: Elem;

impl<T> From<T> for Free<T> 
where T: Elem {
    fn from(value: T) -> Self {
        Self(value)
    }
}

impl<T> Elem for Free<T> 
where T: Elem { 
    fn math_symbol() -> String {
        T::math_symbol()
    }
}

impl<T> Gen for Free<T> 
where T: Elem + Hash + Ord {}