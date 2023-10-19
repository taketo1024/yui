use std::fmt::Display;
use std::ops::{Add, Sub};
use std::hash::Hash;

use derive_more::{Display, Add, Sub};
use num_traits::Zero;
pub trait Deg: Sized + Display + Clone + Copy + PartialEq + Eq + Hash + Zero + Add<Output = Self> + Sub<Output = Self> + 'static {}


impl Deg for isize  {}
impl Deg for usize  {}

macro_rules! make2 {
    ($name:ident, $t:ty) => {
        #[allow(non_camel_case_types)]
        #[derive(Display, Debug, Clone, Copy, PartialEq, Eq, Hash, Add, Sub)]
        #[display(fmt = "({}, {})", _0, _1)]
        pub struct $name(pub $t, pub $t);
        
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

impl Deg for isize2 {}
impl Deg for usize2 {}

macro_rules! make3 {
    ($name:ident, $t:ty) => {
        #[allow(non_camel_case_types)]
        #[derive(Display, Debug, Clone, Copy, PartialEq, Eq, Hash, Add, Sub)]
        #[display(fmt = "({}, {}, {})", _0, _1, _2)]
        pub struct $name(pub $t, pub $t, pub $t);
        
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

impl Deg for isize3 {}
impl Deg for usize3 {}