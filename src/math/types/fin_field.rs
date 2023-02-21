#![allow(non_upper_case_globals)]

use std::fmt::Display;
use std::iter::{Sum, Product};
use std::ops::{Add, Neg, Sub, Mul, Div, Rem, AddAssign, SubAssign, MulAssign, DivAssign, RemAssign};
use std::str::FromStr;
use num_traits::{Zero, One};
use auto_impl_ops::auto_ops;

use crate::math::traits::{AlgBase, AddMonOps, AddGrpOps, MonOps, RingOps, FieldOps, EucRingOps, AddMon, AddGrp, Mon, Ring, EucRing, Field};

pub type I = i32;

#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
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

impl<const p: I> Display for FF<p> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
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

impl_unop!(Neg, neg);
impl_binop!(Add, add);
impl_binop!(Sub, sub);
impl_binop!(Mul, mul);

forward_accum!(Sum, sum, AddAssign, add_assign, zero);
forward_accum!(Product, product, MulAssign, mul_assign, one);

#[auto_ops]
impl<'a, const p: I> Div<&'a FF<p>> for FF<p> {
    type Output = Self;
    fn div(self, rhs: &'a Self) -> Self::Output {
        assert!(!rhs.is_zero());
        self * rhs.inv().unwrap()
    }
}

#[auto_ops]
impl<'a, const p: I> Rem<&'a FF<p>> for FF<p> {
    type Output = Self;
    fn rem(self, rhs: &'a Self) -> Self::Output {
        assert!(!rhs.is_zero());
        FF::zero() // MEMO: FF<p> is a field. 
    }
}

decl_alg_ops!(AddMonOps);
decl_alg_ops!(AddGrpOps);
decl_alg_ops!(MonOps);
decl_alg_ops!(RingOps);
decl_alg_ops!(EucRingOps);
decl_alg_ops!(FieldOps);

impl<const p: I> AlgBase for FF<p> {
    fn symbol() -> String {
        format!("F{p}")
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
            assert!(d.abs().is_one());
            
            let inv = if d.is_one() {
                Self::new(x)
            } else { 
                Self::new(-x)
            };
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

macro_rules! impl_binop {
    ($trait:ident, $method:ident) => {
        #[auto_ops]
        impl<'a, const p: I> $trait<&'a FF<p>> for FF<p> {
            type Output = Self;
            fn $method(self, rhs: &'a Self) -> Self::Output {
                Self::new(self.0.$method(&rhs.0))
            }
        }
    }
}

macro_rules! forward_accum {
    ($trait:ident, $method:ident, $accum_trait:ident, $accum_method:ident, $accum_init:ident) => {
        impl<const p: I> $trait for FF<p> {
            fn $method<Iter: Iterator<Item = Self>>(iter: Iter) -> Self {
                let mut res = Self::$accum_init();
                for r in iter { Self::$accum_method(&mut res, r) }
                return res;
            }
        }

        impl<'a, const p: I> $trait<&'a FF<p>> for FF<p> {
            fn $method<Iter: Iterator<Item = &'a FF<p>>>(iter: Iter) -> Self {
                let mut res = Self::$accum_init();
                for r in iter { Self::$accum_method(&mut res, r) }
                return res;
            }
        }
    }
}

macro_rules! decl_alg_ops {
    ($trait:ident) => {
        impl<const p: I> $trait for FF<p> {}
        impl<'a, const p: I> $trait<FF<p>> for &'a FF<p> {}
    };
}

pub(self) use {impl_unop, impl_binop, forward_accum, decl_alg_ops};

mod tests { 
    #![allow(unused)]
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