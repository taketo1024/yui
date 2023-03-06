use std::ops::Add;
use num_traits::Zero;

pub trait PolyDeg: Add + Zero {
    fn is_add_unit(&self) -> bool;
    fn add_inv(&self) -> Option<Self>;
}

macro_rules! impl_polydeg_unsigned {
    ($t:ty) => {
        impl PolyDeg for $t {
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

macro_rules! impl_polydeg_signed {
    ($t:ty) => {
        impl PolyDeg for $t {
            fn is_add_unit(&self) -> bool { 
                true
            }

            fn add_inv(&self) -> Option<Self> {
                Some(-self)
            }
        }
    };
}

impl_polydeg_unsigned!(usize);
impl_polydeg_signed!  (isize);

pub(crate) use {impl_polydeg_unsigned, impl_polydeg_signed};