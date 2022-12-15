use std::fmt::{Debug, Display};
use std::ops::{Add, Neg, Sub, Mul, Rem, Div, AddAssign, SubAssign, MulAssign, RemAssign, DivAssign};
use std::hash::Hash;
use std::iter::{Sum, Product};
use is_even::IsEven;
use num_traits::{Zero, One};
use super::sign::Sign;

// Set Elements
pub trait MathElem: 
    Debug + Display + Clone + PartialEq + Eq + Hash
{
    fn symbol() -> String;
}

macro_rules! impl_math_elem {
    ($type:ident, $symbol:literal) => {
        impl MathElem for $type {
            fn symbol() -> String {
                String::from($symbol)
            }
        }
    };
}

// Additive Monoids 

pub trait AddMonOps<T>: 
    Sized + Add<Output = T>
{}

pub trait AddMon: 
    MathElem + AddMonOps<Self> + AddAssign + Sum<Self> + Zero
where 
    for<'a> &'a Self: AddMonOps<Self>,
    for<'a> Self: Sum<&'a Self>
{}

macro_rules! impl_add_mon {
    ($type:ident) => {
        impl AddMonOps<$type> for $type {}
        impl<'a> AddMonOps<$type> for &'a $type {}
        impl AddMon for $type {}
    };
}

// Additive Groups 

pub trait AddGrpOps<T>: 
    AddMonOps<T> + Neg<Output = T> + Sub<Output = T>
{}

pub trait AddGrp: 
    AddMon + AddGrpOps<Self> + SubAssign + From<Sign>
where 
    for<'a> &'a Self: AddGrpOps<Self>
{}

macro_rules! impl_add_grp {
    ($type:ident) => {
        impl_add_mon!($type);
        impl AddGrpOps<$type> for $type {}
        impl<'a> AddGrpOps<$type> for &'a $type {}
        impl AddGrp for $type {}
        impl From<Sign> for $type { 
            fn from(e: Sign) -> Self {
                if e.is_positive() { Self::one() } else { -Self::one() }
            }
        }
    };
}

// Monoids (multiplicative)

pub trait MonOps<T>: 
    Sized + Mul<Output = T> 
{}

pub trait Mon: 
    MathElem + MonOps<Self> + MulAssign + Product<Self> + One
where
    for<'a> &'a Self: MonOps<Self>,
    for<'a> Self: Product<&'a Self>
{}

macro_rules! impl_mon {
    ($type:ident) => {
        impl MonOps<$type> for $type {}
        impl<'a> MonOps<$type> for &'a $type {}
        impl Mon for $type {}
    };
}

// Rings 
pub trait RingOps<T>: 
    AddGrpOps<T> + MonOps<T>
{}

pub trait RingMethods: 
    Sized
{
    fn inv(&self) -> Option<Self>;
    fn is_unit(&self) -> bool;
    fn normalizing_unit(&self) -> Self;
}

pub trait Ring: 
    AddGrp + Mon + RingOps<Self> + RingMethods + One 
where
    for<'a> &'a Self: RingOps<Self>
{}

macro_rules! impl_ring {
    ($type:ident) => {
        impl_add_grp!($type);
        impl_mon!($type);
        impl RingOps<$type> for $type {}
        impl<'a> RingOps<$type> for &'a $type {}
        impl Ring for $type {}
    };
}

macro_rules! impl_ring_methods_integer {
    ($type:ident) => {
        impl RingMethods for $type {
            fn inv(&self) -> Option<Self> {
                match self.is_unit() {
                    true => Some(self.clone()),
                    false => None
                }
            }

            fn is_unit(&self) -> bool {
                self == &1 || self == &-1
            }

            fn normalizing_unit(&self) -> Self {
                match self >= &0 { 
                    true  => 1,
                    false => -1
                }
            }
        }
    };
}

// Euclidean Rings
pub trait EucRingOps<T>: 
    RingOps<T> + Rem<Output = T> + Div<Output = T>
{}

pub trait EucRing: 
    Ring + EucRingOps<Self> + RemAssign + DivAssign
where 
    for<'a> &'a Self: EucRingOps<Self>
{
    fn divides(&self, y: &Self) -> bool { 
        !self.is_zero() && (y % self).is_zero()
    }

    fn gcd(x: &Self, y: &Self) -> Self {
        if x.is_zero() && y.is_zero() { return Self::zero() }
        if x.divides(y) { return x.clone() }
        if y.divides(x) { return y.clone() }

        let (mut x, mut y) = (x.clone(), y.clone());

        while !y.is_zero() {
            let r = &x % &y;
            (x, y) = (y, r);
        }

        let u = x.normalizing_unit();

        match u.is_one() { 
            true  => x,
            false => x * u
        }
    }

    fn gcdx(x: &Self, y: &Self) -> (Self, Self, Self) {
        if x.is_zero() && y.is_zero() { return (Self::zero(), Self::zero(), Self::zero()) }
        if x.divides(y) { return (x.clone(), Self::one(), Self::zero()) }
        if y.divides(x) { return (y.clone(), Self::zero(), Self::one()) }

        let (mut x,  mut y)  = (x.clone(), y.clone());
        let (mut s0, mut s1) = (Self::one(),  Self::zero());
        let (mut t0, mut t1) = (Self::zero(), Self::one() );

        while !y.is_zero() {
            let q = &x / &y;
            let r = &x % &y;

            (x, y) = (y, r);
            (s1, s0) = (s0 - &q * &s1, s1);
            (t1, t0) = (t0 - &q * &t1, t1);
        }

        (x, s0, t0)
    }
}

macro_rules! impl_euc_ring {
    ($type:ident) => {
        impl_ring!($type);
        impl EucRingOps<$type> for $type {}
        impl<'a> EucRingOps<$type> for &'a $type {}
        impl EucRing for $type {}
    };
}

macro_rules! impl_euc_ring_integer {
    ($type:ident) => {
        impl_math_elem!($type, "Z");
        impl_ring_methods_integer!($type);
        impl_euc_ring!($type);
    }
}

impl_euc_ring_integer!(i32);
impl_euc_ring_integer!(i64);

// PowMod2

pub trait PowMod2<Rhs> 
where Self: One, Rhs: IsEven {
    fn pow_mod2(self, rhs: Rhs) -> Self { 
        match rhs.is_even() {
            true  => Self::one(),
            false => self
        }
    }
}

impl<T, Rhs> PowMod2<Rhs> for T
where T: One, Rhs: IsEven {}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn check_add_mon() {
        fn check<T>() where T: AddMon, for<'a> &'a T: AddMonOps<T> {}
        check::<i32>();
        check::<i64>();
    }
    
    #[test]
    fn check_add_grp() {
        fn check<T>() where T: AddGrp, for<'a> &'a T: AddGrpOps<T> {}
        check::<i32>();
        check::<i64>();
    }
    
    #[test]
    fn check_mon() {
        fn check<T>() where T: Mon, for<'a> &'a T: MonOps<T> {}
        check::<i32>();
        check::<i64>();
    }

    #[test]
    fn check_ring() {
        fn check<T>() where T: Ring, for<'a> &'a T: RingOps<T> {}
        check::<i32>();
        check::<i64>();
    }

    #[test]
    fn check_eucring() {
        fn check<T>() where T: EucRing, for<'a> &'a T: EucRingOps<T> {}
        check::<i32>();
        check::<i64>();
    }

    #[test]
    fn int_is_unit() { 
        assert!(1.is_unit());
        assert!((-1).is_unit());
        assert_eq!(2.is_unit(), false);
    }

    #[test]
    fn int_inv() { 
        assert_eq!(1.inv(), Some(1));
        assert_eq!((-1).inv(), Some(-1));
        assert_eq!(2.inv(), None);
    }

    #[test]
    fn int_normalizing_unit() { 
        assert_eq!(1.normalizing_unit(), 1);
        assert_eq!((-1).normalizing_unit(), -1);
        assert_eq!(2.normalizing_unit(), 1);
    }

    #[test]
    fn int_divides() {
        assert!(2.divides(&4));
        assert_eq!(3.divides(&4), false);
        assert_eq!(0.divides(&1), false);
    }

    #[test]
    fn gcd_i32() {
        let (a, b) = (240, 46);
        let d = i32::gcd(&a, &b);
        assert_eq!(d, 2);

        let (a, b) = (24, 0);
        let d = i32::gcd(&a, &b);
        assert_eq!(d, 24);

        let (a, b) = (0, 0);
        let d = i32::gcd(&a, &b);
        assert_eq!(d, 0);
    }

    #[test]
    fn gcdx_i32() {
        let (a, b) = (240, 46);
        let (d, s, t) = i32::gcdx(&a, &b);
        assert_eq!(d, 2);
        assert_eq!(&s * &a + &t * &b, d);

        let (a, b) = (24, 0);
        let (d, s, t) = i32::gcdx(&a, &b);
        assert_eq!(d, 24);
        assert_eq!(s, 1);
        assert_eq!(t, 0);

        let (a, b) = (0, 0);
        let (d, s, t) = i32::gcdx(&a, &b);
        assert_eq!(d, 0);
        assert_eq!(s, 0);
        assert_eq!(t, 0);
    }

    #[test]
    fn pow_mod_2() {
        assert_eq!((-1).pow_mod2(10), 1);
        assert_eq!((-1).pow_mod2(-11), -1);
    }
}