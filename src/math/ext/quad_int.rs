use std::fmt::Debug;
use std::iter::{Sum, Product};
use std::{fmt::Display};
use std::ops::{Add, Neg, Sub, Mul, AddAssign, SubAssign, MulAssign};
use num_traits::{Zero, One};
use super::int_ext::{Integer, IntOps};
use super::super::traits::*;

#[derive(Clone, Default, PartialEq, Eq, Debug)]
pub struct QuadInt<I, const D: i32>(I, I)
where I: Integer, for<'x> &'x I: IntOps<I>;

pub type GaussInt<I> = QuadInt<I, -1>;
pub type EisenInt<I> = QuadInt<I, -3>;

impl<I, const D: i32> QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    pub fn new(a: I, b: I) -> Self { 
        assert!(D % 4 != 0);
        Self(a, b)
    }

    pub fn is_rational(&self) -> bool { 
        self.1.is_zero()
    }

    pub fn left(&self) -> &I { 
        &self.0
    }

    pub fn right(&self) -> &I { 
        &self.1
    }

    pub fn pair_into(self) -> (I, I) { 
        (self.0, self.1)
    }

    pub fn pair(&self) -> (&I, &I) { 
        (&self.0, &self.1)
    }
}

impl<I, const D: i32> From<I> for QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    fn from(i: I) -> Self {
        Self::new(i, I::zero())
    }
}

impl<I, const D: i32> Symbol for QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    fn symbol() -> String {
        if D == -1 { 
            String::from("Z[i]")
        } else {
            String::from("Z[ω]")
        }
    }
}

impl<I, const D: i32> Display for QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        todo!()
    }
}

impl<I, const D: i32> Zero for QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    fn zero() -> Self {
        Self::new(I::zero(), I::zero())
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero() && self.1.is_zero()
    }
}

impl<I, const D: i32> One for QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    fn one() -> Self {
        Self::new(I::one(), I::zero())
    }

    fn is_one(&self) -> bool {
        self.0.is_one() && self.1.is_zero()
    }
}

impl<I, const D: i32> Mul<I> for QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    type Output = Self;

    fn mul(self, rhs: I) -> Self::Output {
        &self * &rhs
    }
}

impl<'a, I, const D: i32> Mul<&'a I> for &'a QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    type Output = QuadInt<I, D>;

    fn mul(self, rhs: &'a I) -> Self::Output {
        let (a, b) = self.pair();
        QuadInt(a * rhs, b * rhs)
    }
}

impl<I, const D: i32> Mul for QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        &self * &rhs
    }
}

impl<'a, I, const D: i32> Mul for &'a QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    type Output = QuadInt<I, D>;

    fn mul(self, rhs: Self) -> Self::Output {
        if rhs.is_rational() { 
            return self * &rhs.0
        }

        if self.is_rational() { 
            return rhs * &self.0
        } 

        let (a, b) = self.pair();
        let (c, d) = rhs.pair();

        match D.rem_euclid(4) { 
            1 => {
                let e = I::from((D - 1) / 4);
                let x = a * c + b * d * e;
                let y = a * d + b * c + b * d;
                QuadInt(x, y)
            },
            2 | 3 => {
                let e = I::from(D);
                let x = a * c + b * d * e;
                let y = a * d + b * c;
                QuadInt(x, y)
            },
            _ => panic!()
        }    
    }
}

impl<I, const D: i32> MulAssign for QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = &*self * &rhs;
    }
}

impl<'a, I, const D: i32> MulAssign<&'a QuadInt<I, D>> for QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    fn mul_assign(&mut self, rhs: &'a Self) {
        *self = &*self * rhs;
    }
}

macro_rules! impl_unop {
    ($trait:ident, $method:ident) => {
        impl<I, const D: i32> $trait for QuadInt<I, D>
        where I: Integer, for<'x> &'x I: IntOps<I> {
            type Output = Self;

            fn $method(self) -> Self::Output {
                let (a, b) = self.pair_into();
                Self(I::$method(a), I::$method(b))
            }
        }

        impl<'a, I, const D: i32> $trait for &'a QuadInt<I, D>
        where I: Integer, for<'x> &'x I: IntOps<I> {
            type Output = QuadInt<I, D>;

            fn $method(self) -> Self::Output {
                let (a, b) = self.pair();
                QuadInt(<&'a I>::$method(a), <&'a I>::$method(b))
            }
        }
    };
}

macro_rules! impl_binop {
    ($trait:ident, $method:ident) => {
        impl<I, const D: i32> $trait for QuadInt<I, D>
        where I: Integer, for<'x> &'x I: IntOps<I> {
            type Output = Self;

            fn $method(self, rhs: Self) -> Self::Output {
                let (a, b) = self.pair_into();
                let (c, d) =  rhs.pair_into();
                Self(I::$method(a, c), I::$method(b, d))
            }
        }

        impl<'a, I, const D: i32> $trait for &'a QuadInt<I, D>
        where I: Integer, for<'x> &'x I: IntOps<I> {
            type Output = QuadInt<I, D>;

            fn $method(self, rhs: Self) -> Self::Output {
                let (a, b) = self.pair();
                let (c, d) =  rhs.pair();
                QuadInt(<&'a I>::$method(a, c), <&'a I>::$method(b, d))
            }
        }
    };
}

macro_rules! impl_assop {
    ($trait:ident, $method:ident) => {
        impl<I, const D: i32> $trait for QuadInt<I, D>
        where I: Integer, for<'x> &'x I: IntOps<I> {
            fn $method(&mut self, rhs: Self) {
                let (a, b) = (&mut self.0, &mut self.1);
                let (c, d) = rhs.pair();
                a.$method(c);
                b.$method(d);
            }
        }

        impl<'a, I, const D: i32> $trait<&'a Self> for QuadInt<I, D>
        where I: Integer, for<'x> &'x I: IntOps<I> {
            fn $method(&mut self, rhs: &'a Self) {
                let (a, b) = (&mut self.0, &mut self.1);
                let (c, d) = rhs.pair();
                a.$method(c);
                b.$method(d);
            }
        }
    };
}

macro_rules! impl_accum {
    ($trait:ident, $method:ident, $accum_trait:ident, $accum_method:ident, $accum_init:ident) => {
        impl<I, const D: i32> $trait for QuadInt<I, D>
        where I: Integer, for<'x> &'x I: IntOps<I> {
            fn $method<Iter: Iterator<Item = Self>>(iter: Iter) -> Self {
                let mut res = Self::$accum_init();
                for r in iter { Self::$accum_method(&mut res, r) }
                return res;
            }
        }

        impl<'a, I, const D: i32> $trait<&'a QuadInt<I, D>> for QuadInt<I, D>
        where I: Integer, for<'x> &'x I: IntOps<I> {
            fn $method<Iter: Iterator<Item = &'a QuadInt<I, D>>>(iter: Iter) -> Self {
                let mut res = Self::$accum_init();
                for r in iter { Self::$accum_method(&mut res, r) }
                return res;
            }
        }
    }
}

impl_unop!(Neg, neg);
impl_binop!(Add, add);
impl_binop!(Sub, sub);
impl_assop!(AddAssign, add_assign);
impl_assop!(SubAssign, sub_assign);
impl_accum!(Sum, sum, AddAssign, add_assign, zero);
impl_accum!(Product, product, MulAssign, mul_assign, one);

impl<I, const D: i32> AlgBase for QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {}

impl<I, const D: i32> AddMonOps<Self> for QuadInt<I, D> 
where I: Integer, for<'x> &'x I: IntOps<I> {}

impl<'a, I, const D: i32> AddMonOps<QuadInt<I, D>> for &'a QuadInt<I, D> 
where I: Integer, for<'x> &'x I: IntOps<I> {}

impl<I, const D: i32> AddMon for QuadInt<I, D> 
where I: Integer, for<'x> &'x I: IntOps<I> {}

impl<I, const D: i32> AddGrpOps<Self> for QuadInt<I, D> 
where I: Integer, for<'x> &'x I: IntOps<I> {}

impl<'a, I, const D: i32> AddGrpOps<QuadInt<I, D>> for &'a QuadInt<I, D> 
where I: Integer, for<'x> &'x I: IntOps<I> {}

impl<I, const D: i32> AddGrp for QuadInt<I, D> 
where I: Integer, for<'x> &'x I: IntOps<I> {}

impl<I, const D: i32> MonOps<Self> for QuadInt<I, D> 
where I: Integer, for<'x> &'x I: IntOps<I> {}

impl<'a, I, const D: i32> MonOps<QuadInt<I, D>> for &'a QuadInt<I, D> 
where I: Integer, for<'x> &'x I: IntOps<I> {}

impl<I, const D: i32> Mon for QuadInt<I, D> 
where I: Integer, for<'x> &'x I: IntOps<I> {}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_traits::{Zero, One};
    use super::super::super::traits::*;
    use crate::math::ext::quad_int::QuadInt;

    #[test]
    fn check() { 
        fn check<T>() where T: AddMon, for<'a> &'a T: AddMonOps<T> {}

        type A = QuadInt<i32, -1>;
        type B = QuadInt<i64, -1>;
        type C = QuadInt<BigInt, -1>;

        check::<A>();
        check::<B>();
        check::<C>();
    }
 
    #[test]
    fn zero() { 
        type A = QuadInt<i32, -3>;
        let a = A::new(1, 3);
        let b = A::zero();
        let c = a + b;
        assert_eq!(c, A::new(1, 3));

        let a = A::new(0, 0);
        let b = A::new(0, 1);
        assert_eq!(a.is_zero(), true);
        assert_eq!(b.is_zero(), false);
    }

    #[test]
    fn one() { 
        type A = QuadInt<i32, -3>;

        let a = A::new(1, 3);
        let b = A::one();
        let c = a * b;
        assert_eq!(c, A::new(1, 3));

        let a = A::new(1, 0);
        let b = A::new(0, 1);
        let c = A::new(1, 1);

        assert_eq!(a.is_one(), true);
        assert_eq!(b.is_one(), false);
        assert_eq!(c.is_one(), false);
    }

    #[test]
    fn add() { 
        type A = QuadInt<i32, -3>;
        let a = A::new(1, 3);
        let b = A::new(-3, 2);
        let c = a + b;
        assert_eq!(c, A::new(-2, 5));
    }

    #[test]
    fn sub() { 
        type A = QuadInt<i32, -3>;
        let a = A::new(1, 3);
        let b = A::new(-3, 2);
        let c = a - b;
        assert_eq!(c, A::new(4, 1));
    }

    #[test]
    fn neg() { 
        type A = QuadInt<i32, -3>;
        let a = A::new(1, 3);
        assert_eq!(-a, A::new(-1, -3));
    }

    #[test]
    fn mul_eisen() { 
        type A = QuadInt<i32, -3>; // EisenInt
        let a = A::new(1, 3);
        let b = A::new(2, -1);
        let c = a * b;
        assert_eq!(c, A::new(5, 2));
    }

}