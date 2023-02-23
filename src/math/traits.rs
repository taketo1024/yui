use std::fmt::{Debug, Display};
use std::ops::{Add, Neg, Sub, Mul, Rem, Div, AddAssign, SubAssign, MulAssign, RemAssign, DivAssign};
use std::iter::{Sum, Product};
use is_even::IsEven;
use num_traits::{Zero, One};
use super::types::sign::Sign;

pub trait AlgBase: 
    Default + 
    PartialEq + 
    Eq + 
    Clone + 
    Send + 
    Sync + 
    Display + 
    Debug + 
    'static
{
    fn symbol() -> String; // TODO rename to `set_symbol`
}

// Additive Monoids 

pub trait AddMonOps<T = Self>: 
    Sized + 
    Add<Output = T>
{}

pub trait AddMon: 
    AlgBase + 
    AddMonOps + 
    AddAssign + 
    for<'a> AddAssign<&'a Self> + 
    Sum<Self> + 
    for<'a> Sum<&'a Self> +
    Zero
where 
    for<'a> &'a Self: AddMonOps<Self>
{}

// Additive Groups 

pub trait AddGrpOps<T = Self>: 
    AddMonOps<T> + 
    Neg<Output = T> + 
    Sub<Output = T>
{}

pub trait AddGrp: 
    AddMon + 
    AddGrpOps + 
    SubAssign + 
    for<'a> SubAssign<&'a Self>
where 
    for<'a> &'a Self: AddGrpOps<Self>
{}

// Monoids (multiplicative)

pub trait MonOps<T = Self>: 
    Sized + 
    Mul<Output = T> 
{}

pub trait Mon: 
    AlgBase + 
    MonOps + 
    MulAssign + 
    for<'a> MulAssign<&'a Self> + 
    Product<Self> + 
    for<'a> Product<&'a Self> +
    One
where
    for<'a> &'a Self: MonOps<Self>
{}

// Rings 

pub trait RingOps<T = Self>: 
    AddGrpOps<T> + 
    MonOps<T>
{}

pub trait Ring: 
    AddGrp + 
    Mon + 
    RingOps + 
    sprs::MulAcc
where
    for<'a> &'a Self: RingOps<Self>
{
    fn inv(&self) -> Option<Self>;
    fn is_unit(&self) -> bool;
    fn normalizing_unit(&self) -> Self;
    fn from_sign(e: Sign) -> Self {
        if e.is_positive() { Self::one() } else { -Self::one() }
    }
}

// Euclidean Rings

pub trait EucRingOps<T = Self>: 
    RingOps<T> + 
    Div<Output = T> +
    Rem<Output = T> 
{}

pub trait EucRing: 
    Ring + 
    EucRingOps + 
    DivAssign +
    for<'a> DivAssign<&'a Self> +
    RemAssign + 
    for<'a> RemAssign<&'a Self>
where 
    for<'a> &'a Self: EucRingOps<Self>,
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

    fn lcm(x: &Self, y: &Self) -> Self { 
        let g = Self::gcd(x, y);
        x * &(y / &g)
    }
}

// Fields

pub trait FieldOps<T = Self>: 
    EucRingOps<T>
{}

pub trait Field: 
    EucRing + 
    FieldOps
where 
    for<'a> &'a Self: FieldOps<Self> 
{}

// R-Modules

pub trait RModOps<R, S, T>: AddGrpOps<T> + Mul<S, Output = T>
where 
    R: Ring, for<'x> &'x R: RingOps<R>, 
    S: RingOps<R> // S = R or &'a R
{}

pub trait RMod:
    AddGrp + 
    RModOps<Self::R, Self::R, Self> + 
    MulAssign<Self::R> +
    for<'a> MulAssign<&'a Self::R>
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>, 
    for<'a> &'a Self: RModOps<Self::R, &'a Self::R, Self>,
{
    type R;
}

// PowMod2

pub trait PowMod2<Rhs> 
where Self: One, Rhs: IsEven {
    fn pow_mod2(self, rhs: Rhs) -> Self { 
        if rhs.is_even() { Self::one() } else { self }
    }
}

impl<T, Rhs> PowMod2<Rhs> for T
where T: One, Rhs: IsEven {}

// DivRound

pub trait DivRound { 
    fn div_round(&self, rhs: &Self) -> Self;
}