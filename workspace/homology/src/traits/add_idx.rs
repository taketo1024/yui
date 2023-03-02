use std::fmt::Display;
use std::hash::Hash;
use std::ops::{Add, Neg, Sub};

pub trait AdditiveIndex:
    Clone
    + Copy
    + PartialEq
    + Eq
    + Hash
    + Display
    + Add<Output = Self>
    + Sub<Output = Self>
    + Neg<Output = Self>
{}

impl<T> AdditiveIndex for T where
    T: Clone
        + Copy
        + PartialEq
        + Eq
        + Hash
        + Display
        + Add<Output = Self>
        + Sub<Output = Self>
        + Neg<Output = Self>
{}

pub trait AdditiveIndexRange: Iterator + Clone
where
    Self::Item: AdditiveIndex,
{}

impl<T> AdditiveIndexRange for T
where
    T: Iterator + Clone,
    T::Item: AdditiveIndex,
{}
