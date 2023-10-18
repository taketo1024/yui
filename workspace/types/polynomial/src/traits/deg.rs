use std::ops::Add;
use num_traits::Zero;

pub trait MonoDeg: Add + Zero {
    fn is_add_unit(&self) -> bool;
    fn add_inv(&self) -> Option<Self>;
}

macro_rules! impl_deg_unsigned {
    ($t:ty) => {
        impl MonoDeg for $t {
            fn is_add_unit(&self) -> bool { 
                self.is_zero()
            }

            fn add_inv(&self) -> Option<Self> {
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
        impl MonoDeg for $t {
            fn is_add_unit(&self) -> bool { 
                true
            }

            fn add_inv(&self) -> Option<Self> {
                Some(-self)
            }
        }
    };
}

impl_deg_unsigned!(usize);
impl_deg_signed!  (isize);

pub(crate) use {impl_deg_unsigned, impl_deg_signed};