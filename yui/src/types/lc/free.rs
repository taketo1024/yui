use std::hash::Hash;
use derive_more::Display;
use crate::{Elem, ElemBase};
use crate::lc::Gen;

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