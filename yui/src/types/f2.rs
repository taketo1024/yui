use std::fmt::{Debug, Display};
use std::ops::{Add, Neg, Sub, Mul, Div, Rem, AddAssign, SubAssign, MulAssign, DivAssign, RemAssign};
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use auto_impl_ops::auto_ops;

use crate::{Elem, AddMonOps, AddGrpOps, MonOps, RingOps, FieldOps, EucRingOps, AddMon, AddGrp, Mon, Ring, EucRing, Field};

#[derive(Clone, Copy, PartialEq, Eq, Default)]
pub struct FF2(bool);

impl<I> From<I> for FF2
where I: ToPrimitive {
    fn from(a: I) -> Self {
        let b = a.to_i64().unwrap().is_odd();
        Self(b)
    }
}

impl Display for FF2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.0 { 
            write!(f, "1")
        } else { 
            write!(f, "0")
        }
    }
}

impl Debug for FF2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Display::fmt(self, f)
    }
}

impl Zero for FF2 {
    fn zero() -> Self {
        Self(false)
    }

    fn is_zero(&self) -> bool {
        !self.0
    }
}

impl One for FF2 {
    fn one() -> Self {
        Self(true)
    }

    fn is_one(&self) -> bool {
        self.0
    }
}

impl Neg for FF2 {
    type Output = Self;
    fn neg(self) -> Self::Output {
        self
    }
}

impl Neg for &FF2 {
    type Output = FF2;
    fn neg(self) -> Self::Output {
        self.clone()
    }
}

#[auto_ops]
impl<'a, 'b> Add<&'b FF2> for &'a FF2 {
    type Output = FF2;
    fn add(self, rhs: &'b FF2) -> Self::Output {
        FF2(self.0 != rhs.0)
    }
}

#[auto_ops]
impl<'a, 'b> Sub<&'b FF2> for &'a FF2 {
    type Output = FF2;
    fn sub(self, rhs: &'b FF2) -> Self::Output {
        Add::add(self, rhs)
    }
}

#[auto_ops]
impl<'a, 'b> Mul<&'b FF2> for &'a FF2 {
    type Output = FF2;
    fn mul(self, rhs: &'b FF2) -> Self::Output {
        FF2(self.0 && rhs.0)
    }
}

#[auto_ops]
impl<'a, 'b> Div<&'b FF2> for &'a FF2 {
    type Output = FF2;
    fn div(self, rhs: &'b FF2) -> Self::Output {
        assert!(!rhs.is_zero());
        self.clone()
    }
}

#[auto_ops]
impl<'a, 'b> Rem<&'b FF2> for &'a FF2 {
    type Output = FF2;
    fn rem(self, rhs: &'b FF2) -> Self::Output {
        assert!(!rhs.is_zero());
        FF2::zero()
    }
}

macro_rules! impl_alg_ops {
    ($trait:ident) => {
        impl $trait for FF2 {}
        impl<'a> $trait<FF2> for &'a FF2 {}
    };
}

impl_alg_ops!(AddMonOps);
impl_alg_ops!(AddGrpOps);
impl_alg_ops!(MonOps);
impl_alg_ops!(RingOps);
impl_alg_ops!(EucRingOps);
impl_alg_ops!(FieldOps);

impl Elem for FF2 {
    fn math_symbol() -> String {
        String::from("F₂")
    }
}

impl AddMon for FF2 {}
impl AddGrp for FF2 {}
impl Mon for FF2 {}

impl Ring for FF2 {
    fn inv(&self) -> Option<Self> {
        self.is_one().then_some(*self)
    }

    fn is_unit(&self) -> bool {
        !self.is_zero()
    }

    fn normalizing_unit(&self) -> Self {
        Self::one()
    }

    fn c_weight(&self) -> f64 {
        if self.is_zero() { 
            0f64
        } else { 
            1f64
        }
    }
}

impl EucRing for FF2 {}
impl Field for FF2 {}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn init() { 
        let a = FF2::from(0);
        assert_eq!(a.0, false);

        let a = FF2::from(1);
        assert_eq!(a.0, true);

        let a = FF2::from(2);
        assert_eq!(a.0, false);
    }

    #[test]
    fn display() { 
        let a = FF2::zero();
        assert_eq!(a.to_string(), "0");

        let a = FF2::one();
        assert_eq!(a.to_string(), "1");

        let a = FF2::from(2);
        assert_eq!(a.to_string(), "0");
    }

    #[test]
    fn add() { 
        let a = FF2::from(2);
        let b = FF2::from(4);

        assert_eq!(a + b, FF2::from(0));

        let a = FF2::from(3);
        let b = FF2::from(4);

        assert_eq!(a + b, FF2::from(1));

        let a = FF2::from(3);
        let b = FF2::from(5);

        assert_eq!(a + b, FF2::from(0));
    }

    #[test]
    fn add_assign() { 
        let mut a = FF2::from(3);
        a += FF2::from(4);
        assert_eq!(a, FF2::from(1));
    }

    #[test]
    fn neg() { 
        let a = FF2::from(3);
        assert_eq!(-a, FF2::from(1));
    }

    #[test]
    fn sub() { 
        let a = FF2::from(3);
        let b = FF2::from(5);

        assert_eq!(a - b, FF2::from(0));
    }

    #[test]
    fn sub_assign() { 
        let mut a = FF2::from(3);
        a -= FF2::from(4);
        assert_eq!(a, FF2::from(1));
    }

    #[test]
    fn mul() { 
        let a = FF2::from(3);
        let b = FF2::from(4);
        assert_eq!(a * b, FF2::from(0));

        let a = FF2::from(1);
        let b = FF2::from(5);
        assert_eq!(a * b, FF2::from(1));
    }

    #[test]
    fn mul_assign() { 
        let mut a = FF2::from(3);
        a *= FF2::from(4);
        assert_eq!(a, FF2::from(0));
    }

    #[test]
    fn div() { 
        let a = FF2::from(5);
        let b = FF2::from(3);
        assert_eq!(a / b, FF2::from(1));
    }

    #[test]
    fn div_assign() { 
        let mut a = FF2::from(4);
        a /= FF2::from(3);
        assert_eq!(a, FF2::from(0));
    }

    #[test]
    fn rem() { 
        let a = FF2::from(5);
        let b = FF2::from(3);
        assert_eq!(a % b, FF2::zero());
    }

    #[test]
    fn rem_assign() { 
        let mut a = FF2::from(5);
        a %= FF2::from(3);
        assert_eq!(a, FF2::zero());
    }
}