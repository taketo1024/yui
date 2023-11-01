use std::ops::{Add, Sub};
use num_traits::Zero;

pub trait VarDeg: Add + Sub + Zero {
    fn is_negatable(&self) -> bool;
    fn neg_opt(&self) -> Option<Self>;

    fn strict_leq(&self, other: &Self) -> bool;
    fn strict_geq(&self, other: &Self) -> bool { 
        other.strict_leq(self)
    }
}

impl VarDeg for usize {
    fn is_negatable(&self) -> bool { 
        self.is_zero()
    }

    fn neg_opt(&self) -> Option<Self> {
        if self.is_zero() { 
            Some(0)
        } else { 
            None
        }
    }

    fn strict_leq(&self, other: &Self) -> bool {
        self <= other
    }
}

impl VarDeg for isize {
    fn is_negatable(&self) -> bool { 
        true
    }

    fn neg_opt(&self) -> Option<Self> {
        Some(-self)
    }

    fn strict_leq(&self, other: &Self) -> bool {
        self <= other
    }
}