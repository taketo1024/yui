use std::fmt::Display;
use std::hash::Hash;
use std::ops::{Add, Neg, Sub};

pub trait GridIdx:
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

impl<T> GridIdx for T where
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

pub trait GridItr: Iterator + Clone
where
    Self::Item: GridIdx,
{}

impl<T> GridItr for T
where
    T: Iterator + Clone,
    T::Item: GridIdx,
{}
