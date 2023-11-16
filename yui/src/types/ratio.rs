use std::fmt::{Display, Debug};
use std::str::FromStr;
use std::cmp;
use std::iter::{Sum, Product};
use std::ops::{Mul, Add, Sub, Neg, AddAssign, SubAssign, MulAssign, Div, DivAssign, Rem, RemAssign};
use num_traits::{Zero, One};
use auto_impl_ops::auto_ops;
use crate::{EucRing, EucRingOps, Elem, Mon, AddMon, AddGrp, AddMonOps, AddGrpOps, MonOps, RingOps, Ring, FieldOps, Field, Integer, IntOps};

#[derive(Copy, Clone, PartialEq, Eq)]
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
}

impl<T> Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    #[inline]
    pub fn new(numer: T, denom: T) -> Ratio<T> {
        assert!(!denom.is_zero());

        let mut ret = Ratio::new_raw(numer, denom);
        ret.reduce();
        ret
    }

    fn reduce(&mut self) {
        if self.numer.is_zero() {
            if !self.denom.is_one() { 
                self.denom.set_one();
            }
            return;
        }

        let u = self.denom.normalizing_unit();

        if !u.is_one() { 
            self.numer *= &u;
            self.denom *= &u;
        }

        if self.denom.is_one() || self.numer.is_unit() { 
            return
        }

        let g = EucRing::gcd(&self.numer, &self.denom); // normalized

        if !g.is_one() {
            self.numer /= &g;
            self.denom /= &g;
        }
    }
}

impl<T> Ratio<T>
where T: One {
    pub fn from_numer(a: T) -> Self {
        Self::new_raw(a, T::one())
    }
}

impl<T> Ratio<T>
where T: One + PartialEq {
    pub fn is_numer(&self) -> bool { 
        self.denom.is_one()
    }
}

impl<T> From<i32> for Ratio<T>
where T: One + From<i32> {
    fn from(i: i32) -> Self {
        Self::from_numer(T::from(i))
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
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Ok(a) = s.parse::<T>() {
            return Ok(Self::from_numer(a))
        } 
        
        let r = regex::Regex::new(r"(.+)/(.+)").unwrap();
        if let Some(c) = r.captures(s) { 
            let (s1, s2) = (&c[1], &c[2]);
            if let (Ok(a), Ok(b)) = (s1.parse::<T>(), s2.parse::<T>()) {
                return Ok(Self::new(a, b))
            }
        }

        Err(format!("cannot parse string: '{s}'"))
    }
}

impl<T> Default for Ratio<T>
where T: Default + One {
    fn default() -> Self {
        Self::from_numer(T::default())
    }
}

impl<T> Display for Ratio<T>
where T: Display {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fmt(f, self.numer.to_string(), self.denom.to_string())
    }
}

impl<T> Debug for Ratio<T>
where T: Debug {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fmt(f, format!("{:?}", self.numer), format!("{:?}", self.denom))
    }
}

fn fmt(f: &mut std::fmt::Formatter<'_>, numer: String, denom: String) -> std::fmt::Result { 
    fn par(s: String) -> String {
        if s.contains(' ') { 
            format!("({s})")
        } else { 
            s
        }
    }

    let p = par(numer);
    let q = par(denom);

    if &q == "1" { 
        write!(f, "{}", p)
    } else { 
        write!(f, "{}/{}", p, q)
    }
}

impl<T> Zero for Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    fn zero() -> Self {
        Self::from_numer(T::zero())
    }

    fn is_zero(&self) -> bool {
        self.numer.is_zero()
    }
}

impl<T> One for Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    fn one() -> Self {
        Self::from_numer(T::one())
    }

    fn is_one(&self) -> bool {
        self.numer == self.denom
    }
}

macro_rules! impl_add_assign_op {
    ($trait:ident, $method:ident) => {
        #[auto_ops]
        impl<T> $trait<&Ratio<T>> for Ratio<T>
        where T: EucRing, for<'x> &'x T: EucRingOps<T> {
            fn $method(&mut self, rhs: &Ratio<T>) {
                let (_, b) = (&self.numer, &self.denom);
                let (c, d) = ( &rhs.numer,  &rhs.denom);
                
                if rhs.is_zero() { 
                    // do nothing
                } else if self.is_zero() { 
                    self.numer.$method(c);  // 0 -> 0 ± c
                    self.denom = d.clone(); // 1 -> d
                } else if b == d { 
                    self.numer.$method(c);  // a -> a ± c
                    self.reduce()
                } else { 
                    let l = EucRing::lcm(b, d); // l = xb = yd
                    self.numer *= (&l / b);     // a -> xa ± yc
                    self.numer.$method((&l / d) * c);
                    self.denom = l;             // b -> l
                    self.reduce()
                }
            }
        }
    };
}

impl_add_assign_op!(AddAssign, add_assign);
impl_add_assign_op!(SubAssign, sub_assign);

impl<T> Neg for Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Ratio::new(-&self.numer, self.denom)
    }
}

impl<T> Neg for &Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    type Output = Ratio<T>;
    fn neg(self) -> Self::Output {
        Ratio::new(-&self.numer, self.denom.clone())
    }
}

#[auto_ops]
impl<T> MulAssign<&Ratio<T>> for Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    fn mul_assign(&mut self, rhs: &Ratio<T>) {
        let (a, b) = (&self.numer, &self.denom);
        let (c, d) = ( &rhs.numer,  &rhs.denom);

        if self.is_zero() || rhs.is_one() { 
            // do nothing
        } else if rhs.is_zero() { 
            self.set_zero();             // a -> 0, b -> 1
        } else if rhs.is_numer() { 
            let k = EucRing::gcd(b, c);  // b = kb', c = kc'
            self.numer *= c / &k;        // a -> a * c'
            self.denom /= &k;            // b -> b'
        } else if self.is_numer() { 
            let k = EucRing::gcd(a, d);  // a = ka', d = kd'
            self.numer /= &k;            // a -> a' * c
            self.numer *= c;             // 
            self.denom = d / &k;         // 1 ->      d'
        } else {
            let k = EucRing::gcd(a, d);  // a = ka', d = kd'
            let l = EucRing::gcd(b, c);  // b = lb', c = lc'
            self.numer /= &k;            // a -> a' * c'
            self.numer *= c / &l;        //      
            self.denom /= &l;            // b -> b' * d'
            self.denom *= d / &k;        //      
        }
    }
}

#[auto_ops]
impl<T> DivAssign<&Ratio<T>> for Ratio<T>
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    fn div_assign(&mut self, rhs: &Ratio<T>) {
        assert!(!rhs.is_zero());
        *self *= rhs.inv().unwrap()
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

macro_rules! impl_accum {
    ($trait:ident, $method:ident, $accum_trait:ident, $accum_method:ident, $accum_init:ident) => {
        impl<T> $trait for Ratio<T>
        where T: EucRing, for<'x> &'x T: EucRingOps<T> {
            fn $method<Iter: Iterator<Item = Self>>(iter: Iter) -> Self {
                iter.fold(Self::$accum_init(), |mut res, r| { 
                    Self::$accum_method(&mut res, r);
                    res
                })
            }
        }

        impl<'a, T> $trait<&'a Ratio<T>> for Ratio<T>
        where T: EucRing, for<'x> &'x T: EucRingOps<T> {
            fn $method<Iter: Iterator<Item = &'a Ratio<T>>>(iter: Iter) -> Self {
                iter.fold(Self::$accum_init(), |mut res, r| { 
                    Self::$accum_method(&mut res, r);
                    res
                })
            }
        }
    }
}

impl_accum!(Sum, sum, AddAssign, add_assign, zero);
impl_accum!(Product, product, MulAssign, mul_assign, one);

macro_rules! decl_alg_ops {
    ($trait:ident) => {
        impl<T> $trait for Ratio<T>
        where T: EucRing, for<'x> &'x T: EucRingOps<T> {}

        impl<T> $trait<Ratio<T>> for &Ratio<T>
        where T: EucRing, for<'x> &'x T: EucRingOps<T> {}
    };
}

decl_alg_ops!(AddMonOps);
decl_alg_ops!(AddGrpOps);
decl_alg_ops!(MonOps);
decl_alg_ops!(RingOps);
decl_alg_ops!(EucRingOps);
decl_alg_ops!(FieldOps);

impl<T> Elem for Ratio<T> 
where T: EucRing, for<'x> &'x T: EucRingOps<T> {
    fn math_symbol() -> String {
        let t = T::math_symbol();
        if &t == "Z" { 
            String::from("Q")
        } else { 
            format!("Q({})", T::math_symbol())
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
            let inv = Self::new(self.denom.clone(), self.numer.clone());
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

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn math_symbol() {
        assert_eq!(Ratio::<i32>::math_symbol(), "Q");
    }

    #[test]
    fn constants() {
        assert_eq!(Ratio::zero(), Ratio::new_raw(0, 1));
        assert_eq!(Ratio::one(),  Ratio::new_raw(1, 1));
    }

    #[test]
    fn reduce() {
        let a = Ratio::new(0, -4);
        assert_eq!(a.numer, 0);
        assert_eq!(a.denom, 1);

        let a = Ratio::new(-3, 1);
        assert_eq!(a.numer, -3);
        assert_eq!(a.denom, 1);

        let a = Ratio::new(1, -3);
        assert_eq!(a.numer, -1);
        assert_eq!(a.denom, 3);

        let a = Ratio::new(6, -8);
        assert_eq!(a.numer, -3);
        assert_eq!(a.denom, 4);
    }

    #[test]
    fn display() {
        assert_eq!(format!("{}", Ratio::new(-3, 1)), "-3");
        assert_eq!(format!("{}", Ratio::new(-3, 4)), "-3/4");
    }

    #[test]
    fn debug() {
        assert_eq!(format!("{:?}", Ratio::new(-3, 1)), "-3");
        assert_eq!(format!("{:?}", Ratio::new(-3, 4)), "-3/4");
    }

    #[test]
    fn add() { 
        let a = Ratio::new(1, 2);
        let b = Ratio::new(3, 5);
        assert_eq!(a + b, Ratio::new(11, 10));

        let a = Ratio::new(1, 2);
        let o = Ratio::zero();
        assert_eq!(&a + &o, a);
        assert_eq!(&o + &a, a);

        let a = Ratio::new(1, 3);
        let b = Ratio::new(2, 3);
        assert_eq!(a + b, Ratio::new(1, 1));

        let a = Ratio::new(1, 6);
        let b = Ratio::new(1, 3);
        assert_eq!(a + b, Ratio::new(1, 2));
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

        let a = Ratio::new(1, 2);
        let o = Ratio::zero();
        assert_eq!(&a - &o, a);
        assert_eq!(&o - &a, -a);
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
        let b = Ratio::new(-2, 7);
        assert_eq!(a * b, Ratio::new(-3, 35));

        let a = Ratio::new(3, 4);
        let e = Ratio::one();
        assert_eq!(&a * &e, a);
        assert_eq!(&e * &a, a);

        let a = Ratio::new(3, 4);
        let e = -Ratio::one();
        assert_eq!(&a * &e, -a);
        assert_eq!(&e * &a, -a);

        let a = Ratio::new(3, 4);
        let o = Ratio::zero();
        assert_eq!(&a * &o, Ratio::zero());
        assert_eq!(&o * &a, Ratio::zero());
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
        assert!(!a.is_unit());
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