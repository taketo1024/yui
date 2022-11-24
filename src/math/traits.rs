use std::iter::Sum;
use std::ops::{Rem, Neg};
use std::fmt::Debug;
use is_even::IsEven;
use num_traits::{One, Num};

// TODO split into RingOps
pub trait Ring: Clone + Default + Send + Sync + Debug + Num + Neg<Output = Self> + sprs::MulAcc + Sum {
    fn is_unit(&self) -> bool;
    fn normalizing_unit(&self) -> (Self, Self) { (Self::one(), Self::one()) }
}

// TODO use macro

impl Ring for i32 {
    fn is_unit(&self) -> bool {
        self == &1 || self == &-1
    }

    fn normalizing_unit(&self) -> (Self, Self) {
        if self >= &0 { (1, 1) } else { (-1, -1) }
    }
}

impl Ring for i64 {
    fn is_unit(&self) -> bool {
        self == &1 || self == &-1
    }

    fn normalizing_unit(&self) -> (Self, Self) {
        if self >= &0 { (1, 1) } else { (-1, -1) }
    }
}

pub trait EucRing: Ring + Rem {
    fn is_divisible(&self, y: Self) -> bool { 
        (self.clone() % y).is_zero()
    }

    fn gcd(mut x: Self, mut y: Self) -> Self { 
        while !y.is_zero() {
            let r = x.clone() % y.clone();
            (x, y) = (y, r);
        }

        let u = x.normalizing_unit().0;
        if !u.is_one() {
            x * u.clone()
        } else {
            x
        }
    }

    fn gcdx(mut x: Self, mut y: Self) -> (Self, Self, Self) { 
        let (mut s0, mut s1) = (Self::one(),  Self::zero());
        let (mut t0, mut t1) = (Self::zero(), Self::one() );

        while !y.is_zero() {
            let q = x.clone() / y.clone();
            let r = x.clone() % y.clone();
            (x, y) = (y, r);
            (s0, s1) = (s1.clone(), s0 - q.clone() * s1);
            (t0, t1) = (t1.clone(), t0 - q.clone() * t1);
        }

        (x, s0, t0)
    }
}

impl EucRing for i32 {}
impl EucRing for i64 {}

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

#[cfg(test)]
mod tests { 
    use super::EucRing;

    #[test]
    fn test_gcd_i32() {
        let (a, b) = (240, 46);
        let d = EucRing::gcd(a, b);
        assert_eq!(d, 2);
    }

    #[test]
    fn test_gcdx_i32() {
        let (a, b) = (240, 46);
        let (d, s, t) = EucRing::gcdx(a, b);
        assert_eq!(d, 2);
        assert_eq!(s * a + t * b, d);
    }
}