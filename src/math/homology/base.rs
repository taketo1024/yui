use std::fmt::Display;
use std::ops::{Add, Sub};
use std::hash::Hash;

pub trait AdditiveIndex: Clone + Copy + PartialEq + Eq + Hash + Display + Add<Output = Self> + Sub<Output = Self>{}
impl <T> AdditiveIndex for T
where T: Clone + Copy + PartialEq + Eq + Hash + Display + Add<Output = Self> + Sub<Output = Self>{}

pub trait Graded {
    type Index: AdditiveIndex;
    type IndexRange: Iterator<Item = Self::Index>;

    fn in_range(&self, k: Self::Index) -> bool;
    fn range(&self) -> Self::IndexRange;
}