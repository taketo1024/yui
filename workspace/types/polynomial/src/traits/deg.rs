use std::ops::{Add, Sub};
use num_traits::Zero;

pub trait VarDeg: Add + Sub + Zero {
    fn is_negatable(&self) -> bool;
    fn neg_opt(&self) -> Option<Self>;
}

macro_rules! impl_deg_unsigned {
    ($t:ty) => {
        impl VarDeg for $t {
            fn is_negatable(&self) -> bool { 
                self.is_zero()
            }

            fn neg_opt(&self) -> Option<Self> {
                if self.is_zero() { 
                    Some(Self::zero())
                } else { 
                    None
                }
            }
        }
    };
}

macro_rules! impl_deg_signed {
    ($t:ty) => {
        impl VarDeg for $t {
            fn is_negatable(&self) -> bool { 
                true
            }

            fn neg_opt(&self) -> Option<Self> {
                Some(-self)
            }
        }
    };
}

impl_deg_unsigned!(usize);
impl_deg_signed!  (isize);

pub(crate) use {impl_deg_unsigned, impl_deg_signed};