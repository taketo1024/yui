use is_even::IsEven;
use num_traits::{One, Num};

pub trait Ring: Num + Clone + Default + Send + Sync + sprs::MulAcc + std::fmt::Debug {}
impl Ring for i32 {}

pub trait PowMod2<Rhs> {
    type Output;
    fn pow_mod2(self, rhs: Rhs) -> Self::Output;
}

// TODO use macro

impl PowMod2<i32> for i32 { 
    type Output = i32;
    fn pow_mod2(self, rhs: i32) -> i32 { 
        match rhs.is_even() {
            true  => Self::one(),
            false => self
        }
    }
}