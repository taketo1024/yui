use std::fmt::Debug;
use std::iter::{Sum, Product};
use std::str::FromStr;
use std::{fmt::Display};
use std::ops::{Add, Neg, Sub, Mul, AddAssign, SubAssign, MulAssign, Rem, Div, RemAssign, DivAssign};
use num_traits::{Zero, One};
use auto_impl_ops::auto_ops;
use super::super::ext::int_ext::{Integer, IntOps};
use super::super::traits::*;

#[derive(Clone, Default, PartialEq, Eq, Debug)]
pub struct QuadInt<I, const D: i32>(I, I)
where I: Integer, for<'x> &'x I: IntOps<I>;

pub type GaussInt<I> = QuadInt<I, -1>;
pub type EisenInt<I> = QuadInt<I, -3>;

// The algebraic integers of Q(√D), 
// represented as Z[ω] where
//
//   ω =  { (1 + √D)/2 | D ≡ 1    (mod 4)
//        { √D         | D ≡ 2, 3 (mod 4).
//
// A general z ∈ Z[ω] is represented as
//
//   z = a + bω = { (a + b/2) + (b/2)√D | D ≡ 1
//                {        a  +    b √D | D ≡ 2, 3

impl<I, const D: i32> QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    pub fn new(a: I, b: I) -> Self { 
        assert!(D % 4 != 0);
        Self(a, b)
    }

    pub fn omega() -> Self { 
        Self::new(I::zero(), I::one())
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

    // When D ≡ 1,
    //
    //   bar(z) = (a + b/2) - (b/2)√D
    //          = (a + b) - b ω,
    //
    // when D ≡ 2, 3,
    //
    //   bar(z) = a - b √D
    //          = a - b ω
    //

    pub fn conj(&self) -> Self { 
        let (a, b) = self.pair();
        match D.rem_euclid(4) { 
            1     => QuadInt(a + b, -b),
            2 | 3 => QuadInt(a.clone(), -b),
            _     => panic!()
        }    
    }

    // When D ≡ 1,
    //
    //  N(z) = (a + b/2)^2 - (b/2)^2 D
    //       = a^2 + ab + b^2 (1 - D)/4,
    //
    // when D ≡ 2, 3,
    //
    //  bar(z) = a^2 - b^2 D.
    //
    
    pub fn norm(&self) -> I {
        let (a, b) = self.pair();
        match D.rem_euclid(4) { 
            1     => a * a + a * b + b * b * I::from( (1 - D) / 4),
            2 | 3 => a * a - b * b * I::from(D),
            _     => panic!()
        }    
    }
}

impl<I, const D: i32> From<I> for QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    fn from(i: I) -> Self {
        Self::new(i, I::zero())
    }
}

impl<I, const D: i32> FromStr for QuadInt<I, D>
where I: Integer + FromStr, for<'x> &'x I: IntOps<I> {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Ok(a) = s.parse::<I>() {
            Ok(Self::from(a))
        } else if let Ok((a, b)) = parse_tuple::<I, I>(s) { 
            Ok(Self::new(a, b))
        } else { 
            Err(())
        }
    }
}

fn parse_tuple<I, J>(s: &str) -> Result<(I, J), ()>
where I: FromStr, J: FromStr {
    let r = regex::Regex::new(r"\((.+)?,\s*(.+)?\)").unwrap();
    if let Some(c) = r.captures(&s) { 
        let (s1, s2) = (&c[1], &c[2]);
        if let (Ok(a), Ok(b)) = (s1.parse::<I>(), s2.parse::<J>()) {
            return Ok((a, b))
        }
    }
    Err(())
}
 

impl<I, const D: i32> Display for QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (a, b) = self.pair();
        let x = if D == -1 { "i" } else { "ω" };

        if b.is_zero() { 
            write!(f, "{a}")
        } else if a.is_zero() { 
            let b = 
                if b.is_one() { String::from("") } 
                else if (-b).is_one() { String::from("-") }
                else { b.to_string() };
            write!(f, "{b}{x}")
        } else {
            let sign = if !b.is_negative() { "+" } else { "-" };
            let b = 
                if b.is_unit() { String::from("") } 
                else if b.is_negative() { (-b).to_string() } 
                else { b.to_string() };
    
            write!(f, "{a} {sign} {b}{x}")
        }
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

impl_unop!(Neg, neg);
impl_add_op!(Add, add);
impl_add_op!(Sub, sub);

#[auto_ops]
impl<'a, I, const D: i32> Mul for &'a QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    type Output = QuadInt<I, D>;

    fn mul(self, rhs: Self) -> Self::Output {
        // When D ≡ 1,
        //
        //   ω^2 = (D + 1)/4 + √D/2 
        //       = (D - 1)/4 + ω,
        //
        // hence 
        // 
        //    (a + bω)(c + dω) 
        //  = (ac + bd(D - 1)/4) + (ad + bc + bd)ω.
        //
        // When D ≡ 2, 3,
        // 
        //   ω^2 = D,
        //
        // hence 
        // 
        //    (a + bω)(c + dω) 
        //  = (ac + bdD) + (ad + bc)ω.

        let (a, b) = self.pair();
        let (c, d) = rhs.pair();

        if b.is_zero() {
            return QuadInt(a * c, a * d)
        } else if d.is_zero() { 
            return QuadInt(a * c, b * c)
        }

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

impl_accum!(Sum, sum, AddAssign, add_assign, zero);
impl_accum!(Product, product, MulAssign, mul_assign, one);

// Div / Rem for GaussInt (D = -1).

impl<I> DivRound for GaussInt<I>
where I: Integer, for<'x> &'x I: IntOps<I> {
    fn div_round(&self, rhs: &Self) -> Self {
        let norm = rhs.norm();
        let w = self * &rhs.conj();
        let (x, y) = w.pair_into();
        QuadInt( 
            x.div_round(&norm),
            y.div_round(&norm)
        )
    }
}

#[auto_ops]
impl<'a, I> Div for &'a GaussInt<I>
where I: Integer, for<'x> &'x I: IntOps<I> {
    type Output = GaussInt<I>;

    fn div(self, rhs: Self) -> Self::Output {
        self.div_round(rhs)
    }
}

#[auto_ops]
impl<'a, I> Rem for &'a GaussInt<I>
where I: Integer, for<'x> &'x I: IntOps<I> {
    type Output = QuadInt<I, -1>;

    fn rem(self, rhs: Self) -> Self::Output {
        let q = self / rhs;
        self - &(rhs * &q)
    }
}

// Div / Rem for EisenInt (D = -3).

impl<I> DivRound for EisenInt<I>
where I: Integer, for<'x> &'x I: IntOps<I> {
    //  z / w = (x + y ω) / N(w) 
    //        = (x + y) / N(w) + y / N(w) (ω - 1).
    //
    // (m, n) = ([(x + y) / N(w)], [y / N(w) ]).
    //
    // [z / w] = m + n(ω - 1) = (m - n) + nω.

    fn div_round(&self, rhs: &Self) -> Self {
        let norm = rhs.norm();
        let w = self * &rhs.conj();
        let (x, y) = w.pair();
        let (m, n) = ( 
            (x + y).div_round(&norm),
            y.div_round(&norm)
        );
        QuadInt(&m - &n, n)
    }
}

#[auto_ops]
impl<'a, I> Div for &'a EisenInt<I>
where I: Integer, for<'x> &'x I: IntOps<I> {
    type Output = EisenInt<I>;

    fn div(self, rhs: Self) -> Self::Output {
        self.div_round(rhs)
    }
}

#[auto_ops]
impl<'a, I> Rem for &'a EisenInt<I>
where I: Integer, for<'x> &'x I: IntOps<I> {
    type Output = EisenInt<I>;

    fn rem(self, rhs: Self) -> Self::Output {
        let q = self / rhs;
        self - &(rhs * &q)
    }
}

impl_alg_op!(AddMonOps);
impl_alg_op!(AddGrpOps);
impl_alg_op!(MonOps);
impl_alg_op!(RingOps);
impl_alg_op_d!(EucRingOps, -1);
impl_alg_op_d!(EucRingOps, -3);

impl<I, const D: i32> AlgBase for QuadInt<I, D>
where I: Integer, for<'x> &'x I: IntOps<I> {
    fn symbol() -> String {
        if D == -1 { 
            String::from("Z[i]")
        } else {
            format!("Z[√{}]", D)
        }
    }
}

impl<I, const D: i32> AddMon for QuadInt<I, D> 
where I: Integer, for<'x> &'x I: IntOps<I> {}

impl<I, const D: i32> AddGrp for QuadInt<I, D> 
where I: Integer, for<'x> &'x I: IntOps<I> {}

impl<I, const D: i32> Mon for QuadInt<I, D> 
where I: Integer, for<'x> &'x I: IntOps<I> {}

impl<I, const D: i32> Ring for QuadInt<I, D> 
where I: Integer, for<'x> &'x I: IntOps<I> {
    // see: https://en.wikipedia.org/wiki/Quadratic_integer#Units
    fn is_unit(&self) -> bool {
        self.norm().is_unit()
    }

    fn inv(&self) -> Option<Self> {
        if let Some(u) = self.norm().inv() {
            let u = Self::from(u);
            Some(u * self.conj())
        } else { 
            None
        }
    }

    fn normalizing_unit(&self) -> Self {
        let (a, b) = self.pair();
        if D == -1 { 
            if a.is_positive() && !b.is_negative() {         // a > 0, b ≧ 0 -> 1
                Self::one()
            } else if !a.is_positive() && b.is_positive() {  // a ≦ 0, b > 0 -> -i
                -Self::omega()
            } else if a.is_negative() && !b.is_positive() {  // a < 0, b ≦ 0 -> -1
                -Self::one()
            } else if !a.is_negative() && b.is_negative() {  // a ≧ 0, b < 0 -> i
                Self::omega()
            } else {                                         // a = b = 0    -> 1
                Self::one()
            }
        } else if D == -3 { 
            let c = a + b;
            if a.is_positive() && !b.is_negative() {         // a > 0, b ≧ 0     -> 1
                Self::one()
            } else if !a.is_positive() && c.is_positive() {  // a ≦ 0, a + b > 0 -> 1/ω   = 1 - ω
                Self::new(I::one(), -I::one())
            } else if !c.is_positive() && b.is_positive() {  // a + b ≦ 0, b > 0 -> 1/ω^2 = -ω
                -Self::omega()
            } else if a.is_negative() && !b.is_positive() {  // a < 0, b ≦ 0     -> 1/ω^3 = -1
                -Self::one()
            } else if !a.is_negative() && c.is_negative() {  // a ≧ 0, a + b < 0 -> 1/ω^4 = ω - 1
                Self::new(-I::one(), I::one())
            } else if !c.is_negative() && b.is_negative() {  // a + b ≧ 0, b < 0 -> 1/ω^5 = ω
                Self::omega()
            } else { 
                Self::one()
            }
        } else { 
            if a.is_negative() { 
                -Self::one()
            } else { 
                Self::one()
            }
        }
    }
}

impl<I> EucRing for QuadInt<I, -1> 
where I: Integer, for<'x> &'x I: IntOps<I> {}

impl<I> EucRing for QuadInt<I, -3> 
where I: Integer, for<'x> &'x I: IntOps<I> {}

// -- macros -- //

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

macro_rules! impl_add_op {
    ($trait:ident, $method:ident) => {
        #[auto_ops]
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

macro_rules! impl_alg_op {
    ($trait:ident) => {
        impl<I, const D: i32> $trait<Self> for QuadInt<I, D> 
        where I: Integer, for<'x> &'x I: IntOps<I> {}

        impl<'a, I, const D: i32> $trait<QuadInt<I, D>> for &'a QuadInt<I, D> 
        where I: Integer, for<'x> &'x I: IntOps<I> {}
    };
}

macro_rules! impl_alg_op_d {
    ($trait:ident, $d:literal) => {
        impl<I> $trait<Self> for QuadInt<I, $d> 
        where I: Integer, for<'x> &'x I: IntOps<I> {}

        impl<'a, I> $trait<QuadInt<I, $d>> for &'a QuadInt<I, $d> 
        where I: Integer, for<'x> &'x I: IntOps<I> {}
    };
}

use {impl_unop, impl_add_op, impl_accum, impl_alg_op, impl_alg_op_d};

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_traits::{Zero, One};
    use super::super::super::traits::*;
    use crate::math::types::quad_int::QuadInt;

    #[test]
    fn check() { 
        fn check<T>() where T: Ring, for<'a> &'a T: RingOps<T> {}

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
    fn mul_gauss() { 
        type A = QuadInt<i32, -1>; // GaussInt
        let a = A::new(1, 3);
        let b = A::new(2, -1);
        let c = a * b;
        assert_eq!(c, A::new(5, 5));
    }

    #[test]
    fn mul_eisen() { 
        type A = QuadInt<i32, -3>; // EisenInt
        let a = A::new(1, 3);
        let b = A::new(2, -1);
        let c = a * b;
        assert_eq!(c, A::new(5, 2));
    }

    #[test]
    fn norm_gauss() { 
        type A = QuadInt<i32, -1>; // GaussInt
        let a = A::new(3, -2);
        assert_eq!(a.norm(), 13);
    }

    #[test]
    fn norm_eisen() { 
        type A = QuadInt<i32, -3>; // GaussInt
        let a = A::new(3, -2);
        assert_eq!(a.norm(), 7);
    }

    #[test]
    fn conj_gauss() { 
        type A = QuadInt<i32, -1>; // GaussInt
        let a = A::new(3, -2);
        assert_eq!(a.conj(), A::new(3, 2));
    }

    #[test]
    fn conj_eisen() { 
        type A = QuadInt<i32, -3>; // EisenInt
        let a = A::new(3, -2);
        assert_eq!(a.conj(), A::new(1, 2));
    }

    #[test]
    fn unit_gauss() { 
        type A = QuadInt<i32, -1>; // GaussInt
        assert_eq!(A::new(1,  0).is_unit(), true);
        assert_eq!(A::new(0,  1).is_unit(), true);
        assert_eq!(A::new(-1, 0).is_unit(), true);
        assert_eq!(A::new(0, -1).is_unit(), true);
        assert_eq!(A::new(1,  1).is_unit(), false);
    }

    #[test]
    fn unit_eisen() { 
        type A = QuadInt<i32, -3>; // EisenInt
        assert_eq!(A::new(1,  0).is_unit(), true);
        assert_eq!(A::new(0,  1).is_unit(), true);
        assert_eq!(A::new(-1, 1).is_unit(), true);
        assert_eq!(A::new(-1, 0).is_unit(), true);
        assert_eq!(A::new(0, -1).is_unit(), true);
        assert_eq!(A::new(1, -1).is_unit(), true);
        assert_eq!(A::new(1,  1).is_unit(), false);
    }

    #[test]
    fn inv_gauss() { 
        type A = QuadInt<i32, -1>; // GaussInt
        assert_eq!(A::new(1,  0).inv(), Some(A::new(1,  0)));
        assert_eq!(A::new(-1, 0).inv(), Some(A::new(-1, 0)));
        assert_eq!(A::new(0,  1).inv(), Some(A::new(0, -1)));
        assert_eq!(A::new(0, -1).inv(), Some(A::new(0,  1)));
        assert_eq!(A::new(1,  1).inv(), None);
    }

    #[test]
    fn inv_eisen() { 
        type A = QuadInt<i32, -3>; // EisenInt
        assert_eq!(A::new(1,  0).inv(), Some(A::new(1,  0)));
        assert_eq!(A::new(0,  1).inv(), Some(A::new(1, -1)));
        assert_eq!(A::new(-1, 1).inv(), Some(A::new(0, -1)));
        assert_eq!(A::new(-1, 0).inv(), Some(A::new(-1, 0)));
        assert_eq!(A::new(0, -1).inv(), Some(A::new(-1, 1)));
        assert_eq!(A::new(1, -1).inv(), Some(A::new(0,  1)));
        assert_eq!(A::new(1,  1).inv(), None);
    }

    #[test]
    fn normalizing_unit_gauss() { 
        type A = QuadInt<i32, -1>; // GaussInt
        assert_eq!(A::new(1,  0).normalizing_unit(), A::new(1,  0));
        assert_eq!(A::new(-1, 0).normalizing_unit(), A::new(-1, 0));
        assert_eq!(A::new(2,  0).normalizing_unit(), A::new(1,  0));
        assert_eq!(A::new(0,  1).normalizing_unit(), A::new(0, -1));
        assert_eq!(A::new(0, -1).normalizing_unit(), A::new(0,  1));
        assert_eq!(A::new(0,  2).normalizing_unit(), A::new(0, -1));
        
        assert_eq!(A::new(1,  1).normalizing_unit(), A::new(1,  0));
        assert_eq!(A::new(-1, 1).normalizing_unit(), A::new(0, -1));
        assert_eq!(A::new(-1,-1).normalizing_unit(), A::new(-1, 0));
        assert_eq!(A::new(1, -1).normalizing_unit(), A::new(0,  1));
    }

    #[test]
    fn normalizing_unit_eisen() { 
        type A = QuadInt<i32, -3>; // EisenInt
        assert_eq!(A::new(1,  0).normalizing_unit(), A::new(1,  0));
        assert_eq!(A::new(0,  1).normalizing_unit(), A::new(1, -1));
        assert_eq!(A::new(-1, 1).normalizing_unit(), A::new(0, -1));
        assert_eq!(A::new(-1, 0).normalizing_unit(), A::new(-1, 0));
        assert_eq!(A::new(0, -1).normalizing_unit(), A::new(-1, 1));
        assert_eq!(A::new(1, -1).normalizing_unit(), A::new(0,  1));
    }
    
    #[test]
    fn rem_gauss() {
        type A = QuadInt<i32, -1>; // GaussInt
        let a = A::new(49, -58);
        let b = A::new(7, 9);
        let q = &a / &b;
        let r = &a % &b;

        assert!(!r.is_zero());
        assert!(r.norm() < a.norm());
        assert_eq!(a, b * q + r);
    }

    #[test]
    fn rem_eisen() {
        type A = QuadInt<i32, -3>; // EisenInt
        let a = A::new(49, -58);
        let b = A::new(7, 9);
        let q = &a / &b;
        let r = &a % &b;

        assert!(!r.is_zero());
        assert!(r.norm() < a.norm());
        assert_eq!(a, b * q + r);
    }
}