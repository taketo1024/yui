#![allow(non_upper_case_globals)]

use std::ops::{Add, Neg, Sub, Mul, Div, Rem, AddAssign, SubAssign, MulAssign, DivAssign, RemAssign};
use std::str::FromStr;
use derive_more::{Display, DebugCustom};
use num_traits::{Zero, One};
use auto_impl_ops::auto_ops;

use crate::{Elem, AddMonOps, AddGrpOps, MonOps, RingOps, FieldOps, EucRingOps, AddMon, AddGrp, Mon, Ring, EucRing, Field};

type I = i32;

#[derive(Clone, Copy, PartialEq, Eq, Default, Display, DebugCustom)]
#[display(fmt= "{}", _0)]
#[debug(fmt= "{}", _0)]
pub struct FF<const p: I>(I);

impl<const p: I> FF<p> { 
    pub fn new(a: I) -> Self { 
        assert!(p > 0);
        Self(a.rem_euclid(p))
    }

    pub fn rep(&self) -> &I { 
        &self.0
    }
}

impl<const p: I> From<I> for FF<p> {
    fn from(a: I) -> Self {
        Self::new(a)
    }
}

impl<const p: I> FromStr for FF<p> {
    type Err = <I as FromStr>::Err;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let a = s.parse::<I>()?;
        Ok(Self::from(a))
    }
}

impl<const p: I> Zero for FF<p> {
    fn zero() -> Self {
        Self(0)
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl<const p: I> One for FF<p> {
    fn one() -> Self {
        Self(1)
    }

    fn is_one(&self) -> bool {
        self.0.is_one()
    }
}

macro_rules! impl_unop {
    ($trait:ident, $method:ident) => {
        impl<const p: I> $trait for FF<p> {
            type Output = Self;
            fn $method(self) -> Self::Output {
                Self::new(self.0.$method())
            }
        }

        impl<'a, const p: I> $trait for &'a FF<p> {
            type Output = FF<p>;
            #[inline]
            fn $method(self) -> Self::Output {
                FF::new(self.0.$method())
            }
        }
    };
}

impl_unop!(Neg, neg);

macro_rules! impl_binop {
    ($trait:ident, $method:ident) => {
        #[auto_ops]
        impl<'a, 'b, const p: I> $trait<&'b FF<p>> for &'a FF<p> {
            type Output = FF<p>;
            fn $method(self, rhs: &'b FF<p>) -> Self::Output {
                FF::new(self.0.$method(&rhs.0))
            }
        }
    }
}

impl_binop!(Add, add);
impl_binop!(Sub, sub);
impl_binop!(Mul, mul);

#[auto_ops]
impl<'a, 'b, const p: I> Div<&'b FF<p>> for &'a FF<p> {
    type Output = FF<p>;
    fn div(self, rhs: &'b FF<p>) -> Self::Output {
        assert!(!rhs.is_zero());
        self * rhs.inv().unwrap()
    }
}

#[auto_ops]
impl<'a, 'b, const p: I> Rem<&'b FF<p>> for &'a FF<p> {
    type Output = FF<p>;
    fn rem(self, rhs: &'b FF<p>) -> Self::Output {
        assert!(!rhs.is_zero());
        FF::zero() // MEMO: FF<p> is a field. 
    }
}

macro_rules! impl_alg_ops {
    ($trait:ident) => {
        impl<const p: I> $trait for FF<p> {}
        impl<'a, const p: I> $trait<FF<p>> for &'a FF<p> {}
    };
}

impl_alg_ops!(AddMonOps);
impl_alg_ops!(AddGrpOps);
impl_alg_ops!(MonOps);
impl_alg_ops!(RingOps);
impl_alg_ops!(EucRingOps);
impl_alg_ops!(FieldOps);

impl<const p: I> Elem for FF<p> {
    fn math_symbol() -> String {
        use crate::util::format::subscript;
        format!("F{}", subscript(p as isize))
    }
}

impl<const p: I> AddMon for FF<p> {}
impl<const p: I> AddGrp for FF<p> {}
impl<const p: I> Mon for FF<p> {}

impl<const p: I> Ring for FF<p> {
    fn inv(&self) -> Option<Self> {
        if self.is_zero() { 
            None
        } else { 
            // 1 = ax + py  ->  ax = 1 mod p. 
            let (d, x, _y) = I::gcdx(&self.0, &p);
            
            assert!(d.is_one());
            
            let inv = Self::new(x);
            Some(inv)
        }
    }

    fn is_unit(&self) -> bool {
        !self.is_zero()
    }

    fn normalizing_unit(&self) -> Self {
        if self.is_zero() { 
            Self::one()
        } else { 
            self.inv().unwrap()
        }
    }
}

impl<const p: I> EucRing for FF<p> {}
impl<const p: I> Field for FF<p> {}

#[cfg(test)]
mod tests { 
    use super::*;

    type F3 = FF<3>;
    type F5 = FF<5>;

    #[test]
    fn init() { 
        let a = F3::new(-7);
        assert_eq!(a.0, 2);

        let a = F5::new(-7);
        assert_eq!(a.0, 3);
    }

    #[test]
    fn display() { 
        let a = F3::new(-7);
        assert_eq!(format!("{}", a), "2");

        let a = F5::new(-7);
        assert_eq!(format!("{}", a), "3");
    }

    #[test]
    fn debug() { 
        let a = F3::new(-7);
        assert_eq!(format!("{:?}", a), "2");

        let a = F5::new(-7);
        assert_eq!(format!("{:?}", a), "3");
    }

    #[test]
    fn add() { 
        let a = F5::new(3);
        let b = F5::new(4);

        assert_eq!(a + b, F5::new(2));
    }

    #[test]
    fn add_assign() { 
        let mut a = F5::new(3);
        a += F5::new(4);
        assert_eq!(a, F5::new(2));
    }

    #[test]
    fn neg() { 
        let a = F5::new(3);
        assert_eq!(-a, F5::new(2));
    }

    #[test]
    fn sub() { 
        let a = F5::new(3);
        let b = F5::new(4);

        assert_eq!(a - b, F5::new(4));
    }

    #[test]
    fn sub_assign() { 
        let mut a = F5::new(3);
        a -= F5::new(4);
        assert_eq!(a, F5::new(4));
    }

    #[test]
    fn mul() { 
        let a = F5::new(3);
        let b = F5::new(4);
        assert_eq!(a * b, F5::new(2));
    }

    #[test]
    fn mul_assign() { 
        let mut a = F5::new(3);
        a *= F5::new(4);
        assert_eq!(a, F5::new(2));
    }

    #[test]
    fn div() { 
        let a = F5::new(4);
        let b = F5::new(3);
        assert_eq!(a / b, F5::new(3));
    }

    #[test]
    fn div_assign() { 
        let mut a = F5::new(4);
        a /= F5::new(3);
        assert_eq!(a, F5::new(3));
    }


    #[test]
    fn rem() { 
        let a = F5::new(4);
        let b = F5::new(3);
        assert_eq!(a % b, F5::zero());
    }

    #[test]
    fn rem_assign() { 
        let mut a = F5::new(4);
        a %= F5::new(3);
        assert_eq!(a, F5::zero());
    }
}