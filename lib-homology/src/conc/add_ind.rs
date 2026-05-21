use std::ops::{Add, Neg, Sub};
use num_traits::Zero;
use yui_core::IndexType;

pub trait AddInd:
    IndexType
    + Copy
    + Zero
    + Add<Output = Self>
    + Sub<Output = Self>
    + Neg<Output = Self>
{}

impl AddInd for isize {}

macro_rules! make {
    ($name:ident, $t:ty, $($idx:tt),+) => {
        #[allow(non_camel_case_types)]
        #[derive(
            Clone, Copy, Default,
            PartialEq, Eq, PartialOrd, Ord, Hash, Debug,
            derive_more::Add, derive_more::Sub, derive_more::Neg,
        )]
        pub struct $name($(pub make!(@unit $idx, $t)),+);

        impl Zero for $name {
            fn zero() -> Self {
                Self($(make!(@unit $idx, 0)),+)
            }

            fn is_zero(&self) -> bool {
                $(self.$idx.is_zero())&&+
            }
        }

        impl From<($(make!(@unit $idx, $t)),+)> for $name {
            fn from(i: ($(make!(@unit $idx, $t)),+)) -> Self {
                Self($(i.$idx),+)
            }
        }

        impl From<$name> for ($(make!(@unit $idx, $t)),+) {
            fn from(i: $name) -> Self {
                ($(i.$idx),+)
            }
        }

        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(f, "({})", [$(self.$idx.to_string()),+].join(", "))
            }
        }

        impl AddInd for $name {}
    };
    // Helper arm: drops the first token and returns the rest verbatim.
    // Lets `$(... $t ...),+` repeat using `$idx` as the binding while reusing `$t`.
    (@unit $_idx:tt, $($body:tt)*) => { $($body)* };
}

make!(isize2, isize, 0, 1);
make!(isize3, isize, 0, 1, 2);
