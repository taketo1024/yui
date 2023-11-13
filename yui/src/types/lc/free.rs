use std::hash::Hash;
use derive_more::Display;
use crate::Elem;

use crate::Gen;

#[derive(Debug, Display, Default, Hash, PartialEq, Eq, Clone, PartialOrd, Ord)]
#[display(fmt = "<{}>", _0)]
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