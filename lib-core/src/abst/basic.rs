use std::fmt::{Debug, Display};
use std::hash::Hash;

/// A type that represents a mathematical object, identified by a symbol
/// returned from [`math_symbol`](MathType::math_symbol).
pub trait MathType:
    Default +
    PartialEq +
    Eq +
    Clone +
    Send +
    Sync +
    Display +
    Debug +
    'static
{
    fn math_symbol() -> String;
}

/// A type usable as an index or key — hashable, totally ordered, and carrying
/// the basic shape bounds. Blanket-implemented for any type satisfying them.
pub trait IndexType:
    Default +
    PartialEq +
    Eq +
    Hash +
    PartialOrd +
    Ord +
    Clone +
    Send +
    Sync +
    Display +
    Debug +
    'static
{}

impl<T> IndexType for T where T:
    Default +
    PartialEq +
    Eq +
    Hash +
    PartialOrd +
    Ord +
    Clone +
    Send +
    Sync +
    Display +
    Debug +
    'static
{}