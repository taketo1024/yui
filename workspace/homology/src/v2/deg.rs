use std::ops::{Add, Sub, Neg};
use std::hash::Hash;

use derive_more::{Add, Sub, Neg};
use num_traits::Zero;

pub trait Deg: Clone + Copy + PartialEq + Eq + Hash + Zero + Add + Sub + Neg {}

impl Deg for isize {}
impl Deg for isize2 {}
impl Deg for isize3 {}

#[allow(non_camel_case_types)]
#[derive(Clone, Copy, PartialEq, Eq, Hash, Add, Sub, Neg)]
pub struct isize2(isize, isize);

impl Zero for isize2 {
    fn zero() -> Self {
        Self(0, 0)
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero() && self.1.is_zero()
    }
}

#[allow(non_camel_case_types)]
#[derive(Clone, Copy, PartialEq, Eq, Hash, Add, Sub, Neg)]
pub struct isize3(isize, isize, isize);

impl Zero for isize3 {
    fn zero() -> Self {
        Self(0, 0, 0)
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero() && self.1.is_zero() && self.2.is_zero()
    }
}