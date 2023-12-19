use std::fmt::Display;
use std::hash::Hash;
use std::ops::{Add, Sub};

use derive_more::{Add, Display, Sub};
use num_traits::Zero;

pub trait GridDeg:
    Sized
    + Display
    + Clone
    + Copy
    + PartialEq
    + Eq
    + Hash
    + Zero
    + Add<Output = Self>
    + Sub<Output = Self>
    + Send
    + Sync
    + 'static
{}

impl GridDeg for isize {}
impl GridDeg for usize {}

macro_rules! make2 {
    ($name:ident, $t:ty) => {
        #[allow(non_camel_case_types)]
        #[derive(Display, Debug, Clone, Copy, Default, PartialEq, Eq, Hash, Add, Sub)]
        #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
        #[display(fmt = "({}, {})", _0, _1)]
        pub struct $name(pub $t, pub $t);

        impl From<($t, $t)> for $name {
            fn from(i: ($t, $t)) -> Self {
                Self(i.0, i.1)
            }
        }

        impl From<$name> for ($t, $t) {
            fn from(i: $name) -> Self {
                (i.0, i.1)
            }
        }

        impl Zero for $name {
            fn zero() -> Self {
                Self(0, 0)
            }

            fn is_zero(&self) -> bool {
                self.0.is_zero() && self.1.is_zero()
            }
        }
    };
}

make2!(isize2, isize);
make2!(usize2, usize);

impl GridDeg for isize2 {}
impl GridDeg for usize2 {}

macro_rules! make3 {
    ($name:ident, $t:ty) => {
        #[allow(non_camel_case_types)]
        #[derive(Display, Debug, Clone, Copy, Default, PartialEq, Eq, Hash, Add, Sub)]
        #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
        #[display(fmt = "({}, {}, {})", _0, _1, _2)]
        pub struct $name(pub $t, pub $t, pub $t);

        impl From<($t, $t, $t)> for $name {
            fn from(i: ($t, $t, $t)) -> Self {
                Self(i.0, i.1, i.2)
            }
        }

        impl From<$name> for ($t, $t, $t) {
            fn from(i: $name) -> Self {
                (i.0, i.1, i.2)
            }
        }

        impl Zero for $name {
            fn zero() -> Self {
                Self(0, 0, 0)
            }

            fn is_zero(&self) -> bool {
                self.0.is_zero() && self.1.is_zero() && self.2.is_zero()
            }
        }
    };
}

make3!(isize3, isize);
make3!(usize3, usize);

impl GridDeg for isize3 {}
impl GridDeg for usize3 {}