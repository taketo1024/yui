use std::fmt::Display;
use std::str::FromStr;
use std::cmp;
use std::iter::{Sum, Product};
use std::ops::{Mul, Add, Sub, Neg, AddAssign, SubAssign, MulAssign, Div, DivAssign, Rem, RemAssign};
use num_traits::{Zero, One};
use auto_impl_ops::auto_ops;
use crate::math::ext::int_ext::{Integer, IntOps};
use crate::math::traits::{EucRing, EucRingOps, Elem, Mon, AddMon, AddGrp, AddMonOps, AddGrpOps, MonOps, RingOps, Ring, FieldOps, Field};

#[derive(Copy, Clone, Debug)]
pub struct Ratio<T> {
    numer: T,
    denom: T,
}

impl<T> Ratio<T> {
    #[inline]
    const fn new_raw(numer: T, denom: T) -> Ratio<T> {
        Ratio { numer, denom }
    }

    #[inline]
    pub const fn numer(&self) -> &T {
        &self.numer
    }

    #[inline]
    pub const fn denom(&self) -> &T {
        &self.denom
    }

    #[inline]
    fn inv_raw(self) -> Self {
        Self::new_raw(self.denom, self.numer)
    }
}

impl<T> From<T> for Ratio<T>
where T: One {
    fn from(a: T) -> Self {
        Self::new_raw(a, T::one())
    }
}

impl<T> Default for Ratio<T>
where T: Default + One {
    fn default() -> Self {
        Self::from(T::default())
    }
}

impl<T> Display for Ratio<T>
where T: Display {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fn par(s: String) -> String {
            if s.contains(' ') { 
                format!("({s})")
            } else { 
                s
            }
        }
        
        let p = par(self.numer.to_string());
        let q = par(self.denom.to_string());

        if &q == "1" { 
            write!(f, "{}", p)
        } else { 
            write!(f, "{}/{}", p, q)
        }
    }
}

impl<T> Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    #[inline]
    pub fn new(numer: T, denom: T) -> Ratio<T> {
        let mut ret = Ratio::new_raw(numer, denom);
        ret.reduce();
        ret
    }

    pub fn reduce(&mut self) {
        if self.denom.is_zero() {
            panic!("denominator == 0");
        }

        if self.numer.is_zero() {
            self.denom.set_one();
            return;
        }

        if self.numer == self.denom {
            self.set_one();
            return;
        }

        let g: T = EucRing::gcd(&self.numer, &self.denom);

        #[inline]
        fn replace_with<T: Zero>(x: &mut T, f: impl FnOnce(T) -> T) {
            let y = core::mem::replace(x, T::zero());
            *x = f(y);
        }

        // self.numer /= g;
        replace_with(&mut self.numer, |x| x / g.clone());

        // self.denom /= g;
        replace_with(&mut self.denom, |x| x / g);

        // normalize denom
        let u = self.denom.normalizing_unit();
        if !u.is_one() {
            replace_with(&mut self.numer, |x| x * u.clone());
            replace_with(&mut self.denom, |x| x * u);
        }
    }
}

impl<T> From<(T, T)> for Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    fn from(pair: (T, T)) -> Self {
        let (p, q) = pair;
        Self::new(p, q)
    }
}

impl<T> FromStr for Ratio<T>
where T: EucRing + FromStr, for<'x> &'x T: EucRingOps<T> {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Ok(a) = s.parse::<T>() {
            return Ok(Self::from(a))
        } 
        
        let r = regex::Regex::new(r"(.+)/(.+)").unwrap();
        if let Some(c) = r.captures(&s) { 
            let (s1, s2) = (&c[1], &c[2]);
            if let (Ok(a), Ok(b)) = (s1.parse::<T>(), s2.parse::<T>()) {
                return Ok(Self::new(a, b))
            }
        }
        Err(())
    }
}

impl<T> PartialEq for Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    // MEMO: it is assumed that p/q is uniquely reduced.
    fn eq(&self, other: &Self) -> bool {
        self.numer == other.numer && 
        self.denom == other.denom
    }
}

impl<T> Eq for Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {}

impl<T> Zero for Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    fn zero() -> Self {
        Self::new(T::zero(), T::one())
    }

    fn is_zero(&self) -> bool {
        self.numer.is_zero()
    }
}

impl<T> One for Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    fn one() -> Self {
        Self::new_raw(T::one(), T::one())
    }

    fn is_one(&self) -> bool {
        self.numer == self.denom
    }
}

impl_unop!(Neg, neg);
impl_add_op!(Add, add);
impl_add_op!(Sub, sub);

#[auto_ops]
impl<'a, 'b, T> Mul<&'b Ratio<T>> for &'a Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    type Output = Ratio<T>;
    fn mul(self, rhs: &'b Ratio<T>) -> Self::Output {
        let gcd_ad = EucRing::gcd(&self.numer, &rhs.denom);
        let gcd_bc = EucRing::gcd(&self.denom, &rhs.numer);

        Ratio::new(
            &self.numer / &gcd_ad * (&rhs.numer / &gcd_bc),
            &self.denom / &gcd_bc * (&rhs.denom / &gcd_ad),
        )
    }
}

#[auto_ops]
impl<'a, 'b, T> Div<&'b Ratio<T>> for &'a Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    type Output = Ratio<T>;
    fn div(self, rhs: &'b Ratio<T>) -> Self::Output {
        assert!(!rhs.is_zero());
        self * rhs.clone().inv_raw()
    }
}

#[auto_ops]
impl<'a, 'b, T> Rem<&'b Ratio<T>> for &'a Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    type Output = Ratio<T>;
    fn rem(self, rhs: &'b Ratio<T>) -> Self::Output {
        assert!(!rhs.is_zero());
        Ratio::zero() // MEMO Frac<T> is a field. 
    }
}

impl_accum!(Sum, sum, AddAssign, add_assign, zero);
impl_accum!(Product, product, MulAssign, mul_assign, one);

decl_alg_ops!(AddMonOps);
decl_alg_ops!(AddGrpOps);
decl_alg_ops!(MonOps);
decl_alg_ops!(RingOps);
decl_alg_ops!(EucRingOps);
decl_alg_ops!(FieldOps);

impl<T> Elem for Ratio<T> 
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    fn set_symbol() -> String {
        let t = T::set_symbol();
        if &t == "Z" { 
            String::from("Q")
        } else { 
            format!("Frac({})", T::set_symbol())
        }
    }
}

impl<T> Mon for Ratio<T> 
where T: EucRing, for<'x> &'x T: EucRingOps<T> {}

impl<T> AddMon for Ratio<T> 
where T: EucRing, for<'x> &'x T: EucRingOps<T> {}

impl<T> AddGrp for Ratio<T> 
where T: EucRing, for<'x> &'x T: EucRingOps<T> {}

impl<T> Ring for Ratio<T> 
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    fn inv(&self) -> Option<Self> {
        if self.is_zero() { 
            None
        } else { 
            let mut inv = self.clone().inv_raw();
            inv.reduce();
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

impl<T> EucRing for Ratio<T> 
where T: EucRing, for<'x> &'x T: EucRingOps<T> {}

impl<T> Field for Ratio<T> 
where T: EucRing, for<'x> &'x T: EucRingOps<T> {}

impl<T> Ratio<T>
where T: Integer, for<'x> &'x T: IntOps<T> {
    pub fn abs(&self) -> Self {
        if self.numer.is_negative() { 
            -self
        } else { 
            self.clone()
        }
    }

    pub fn to_f64(&self) -> f64 { 
        let p = self.numer.to_f64().unwrap();
        let q = self.denom.to_f64().unwrap();
        p / q
    }
}

impl<T> Ord for Ratio<T>
where T: Integer, for<'x> &'x T: IntOps<T> {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        let l = self.to_f64();
        let r = other.to_f64();
        l.total_cmp(&r)
    }
}

impl<T> PartialOrd for Ratio<T> 
where T: Integer, for<'x> &'x T: IntOps<T> {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

macro_rules! impl_unop {
    ($trait:ident, $method:ident) => {
        impl<T> $trait for Ratio<T>
        where T: EucRing, for<'x> &'x T: EucRingOps<T> {
            type Output = Self;
            fn $method(self) -> Self::Output {
                Ratio::new(self.numer.$method(), self.denom)
            }
        }

        impl<'a, T> $trait for &'a Ratio<T>
        where T: EucRing, for<'x> &'x T: EucRingOps<T> {
            type Output = Ratio<T>;
            #[inline]
            fn $method(self) -> Self::Output {
                Ratio::new((&self.numer).$method(), self.denom.clone())
            }
        }
    };
}

macro_rules! impl_add_op {
    ($trait:ident, $method:ident) => {
        #[auto_ops]
        impl<'a, 'b, T> $trait<&'b Ratio<T>> for &'a Ratio<T>
        where T: EucRing, for<'x> &'x T: EucRingOps<T> {
            type Output = Ratio<T>;
            fn $method(self, rhs: &'b Ratio<T>) -> Self::Output {
                if self.denom == rhs.denom {
                    return Ratio::new( 
                        (&self.numer).$method(&rhs.numer), 
                        rhs.denom.clone()
                    );
                }
                let lcm = EucRing::lcm(&self.denom, &rhs.denom);
                let lhs_numer = &self.numer * (&lcm / &self.denom);
                let rhs_numer = &rhs.numer * (&lcm / &rhs.denom);
                Ratio::new(lhs_numer.$method(rhs_numer), lcm)
            }
        }
    };
}

macro_rules! impl_accum {
    ($trait:ident, $method:ident, $accum_trait:ident, $accum_method:ident, $accum_init:ident) => {
        impl<T> $trait for Ratio<T>
        where T: EucRing, for<'x> &'x T: EucRingOps<T> {
            fn $method<Iter: Iterator<Item = Self>>(iter: Iter) -> Self {
                let mut res = Self::$accum_init();
                for r in iter { Self::$accum_method(&mut res, r) }
                return res;
            }
        }

        impl<'a, T> $trait<&'a Ratio<T>> for Ratio<T>
        where T: EucRing, for<'x> &'x T: EucRingOps<T> {
            fn $method<Iter: Iterator<Item = &'a Ratio<T>>>(iter: Iter) -> Self {
                let mut res = Self::$accum_init();
                for r in iter { Self::$accum_method(&mut res, r) }
                return res;
            }
        }
    }
}

macro_rules! decl_alg_ops {
    ($trait:ident) => {
        impl<T> $trait for Ratio<T>
        where T: EucRing, for<'x> &'x T: EucRingOps<T> {}

        impl<T> $trait<Ratio<T>> for &Ratio<T>
        where T: EucRing, for<'x> &'x T: EucRingOps<T> {}
    };
}

pub(self) use {impl_unop, impl_add_op, impl_accum, decl_alg_ops};

#[cfg(test)]
mod tests { 
    use crate::math::types::quad_int::GaussInt;
    use super::*;

    #[test]
    fn set_symbol() {
        assert_eq!(Ratio::<i32>::set_symbol(), "Q");
        assert_eq!(Ratio::<GaussInt<i32>>::set_symbol(), "Frac(Z[i])");
    }

    #[test]
    fn constants() {
        assert_eq!(Ratio::zero(), Ratio::new_raw(0, 1));
        assert_eq!(Ratio::one(),  Ratio::new_raw(1, 1));
    }

    #[test]
    fn reduce() {
        let a = Ratio::new(6, -8);
        assert_eq!(a.numer, -3);
        assert_eq!(a.denom, 4);

        let a = Ratio::new(0, i32::MIN);
        assert_eq!(a.numer, 0);
        assert_eq!(a.denom, 1);
    }

    #[test]
    fn display() {
        assert_eq!(Ratio::new(-3, 1).to_string(), "-3");
        assert_eq!(Ratio::new(-3, 4).to_string(), "-3/4");

        use GaussInt as T;
        assert_eq!(Ratio::new(T::new(1, 1), T::new(2, 3)).to_string(), "(1 + i)/(2 + 3i)");
    }

    #[test]
    fn add() { 
        let a = Ratio::new(1, 2);
        let b = Ratio::new(3, 5);
        assert_eq!(a + b, Ratio::new(11, 10));
    }

    #[test]
    fn add_assign() { 
        let mut a = Ratio::new(1, 2);
        a += Ratio::new(3, 5);

        assert_eq!(a, Ratio::new(11, 10));
    }

    #[test]
    fn neg() { 
        let a = Ratio::new(1, 2);
        assert_eq!(-a, Ratio::new(-1, 2));
    }

    #[test]
    fn sub() { 
        let a = Ratio::new(1, 2);
        let b = Ratio::new(3, 5);
        assert_eq!(a - b, Ratio::new(-1, 10));
    }

    #[test]
    fn sub_assign() { 
        let mut a = Ratio::new(1, 2);
        a -= Ratio::new(3, 5);
        assert_eq!(a, Ratio::new(-1, 10));
    }

    #[test]
    fn mul() { 
        let a = Ratio::new(3, 10);
        let b = Ratio::new(2, 7);
        assert_eq!(a * b, Ratio::new(3, 35));
    }

    #[test]
    fn mul_assign() { 
        let mut a = Ratio::new(3, 10);
        a *= Ratio::new(2, 7);
        assert_eq!(a, Ratio::new(3, 35));
    }

    #[test]
    fn div() { 
        let a = Ratio::new(3, 10);
        let b = Ratio::new(2, 7);
        assert_eq!(a / b, Ratio::new(21, 20));
    }

    #[test]
    fn div_assign() { 
        let mut a = Ratio::new(3, 10);
        a /= Ratio::new(2, 7);
        assert_eq!(a, Ratio::new(21, 20));
    }

    #[test]
    fn rem() { 
        let a = Ratio::new(3, 10);
        let b = Ratio::new(2, 7);
        assert_eq!(a % b, Ratio::zero());
    }

    #[test]
    fn rem_assign() { 
        let mut a = Ratio::new(3, 10);
        a %= Ratio::new(2, 7);
        assert_eq!(a, Ratio::zero());
    }

    #[test]
    fn inv() { 
        let a = Ratio::new(-3, 10);
        assert_eq!(a.inv(), Some(Ratio::new(-10, 3)));

        let a = Ratio::<i32>::zero();
        assert_eq!(a.inv(), None);
    }

    #[test]
    fn is_unit() { 
        let a = Ratio::new(-3, 10);
        assert!(a.is_unit());

        let a = Ratio::<i32>::zero();
        assert_eq!(a.is_unit(), false);
    }

    #[test]
    fn normalizing_unit() { 
        let a = Ratio::new(-3, 10);
        assert_eq!(a.normalizing_unit(), Ratio::new(-10, 3));

        let a = Ratio::<i32>::zero();
        assert_eq!(a.normalizing_unit(), Ratio::one());
    }

    #[test]
    fn cmp() { 
        let a = Ratio::new(3, 5);
        let b = Ratio::new(4, 7);
        assert!(a > b);
    }
}