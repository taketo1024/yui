use std::fmt::Debug;
use std::ops::{Add, Neg, Sub, Mul, Rem};
use num_traits::{Zero, One};
use std::iter::{Sum, Product};
use is_even::IsEven;

// Set Elements
pub trait MathElem: 
    Clone + PartialEq + Eq + Debug
{}

impl<T> MathElem for T where 
    T: Clone + PartialEq + Eq + Debug
{}

// Additive Monoids 

pub trait AddMonOps<T>: 
    Sized + Add<Output = T>
{}

pub trait AddMon: 
    MathElem + AddMonOps<Self> + Sum<Self> + Zero
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
    AddMon + AddGrpOps<Self>
where 
    for<'a> &'a Self: AddGrpOps<Self>
{}

macro_rules! impl_add_grp {
    ($type:ident) => {
        impl_add_mon!($type);
        impl AddGrpOps<$type> for $type {}
        impl<'a> AddGrpOps<$type> for &'a $type {}
        impl AddGrp for $type {}
    };
}

// Monoids (multiplicative)

pub trait MonOps<T>: 
    Sized + Mul<Output = T> 
{}

pub trait Mon: 
    MathElem + MonOps<Self> + Product<Self> + One
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

macro_rules! impl_ring_integer {
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
        impl_ring!($type);
    };
}

impl_ring_integer!(i8);
impl_ring_integer!(i16);
impl_ring_integer!(i32);
impl_ring_integer!(i64);

pub trait PowMod2<Rhs> {
    type Output;
    fn pow_mod2(self, rhs: Rhs) -> Self::Output;
}

macro_rules! impl_powmod2_integer {
    ($type:ident) => {
        impl PowMod2<$type> for $type { 
            type Output = $type;
            fn pow_mod2(self, rhs: $type) -> $type { 
                match rhs.is_even() {
                    true  => Self::one(),
                    false => self
                }
            }
        }
    };
}

impl_powmod2_integer!(i32);
impl_powmod2_integer!(i64);

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn check_add_mon() {
        fn check<T>() where T: AddMon, for<'a> &'a T: AddMonOps<T> {}
        check::<i8>();
        check::<i16>();
        check::<i32>();
        check::<i64>();
    }
    
    #[test]
    fn check_add_grp() {
        fn check<T>() where T: AddGrp, for<'a> &'a T: AddGrpOps<T> {}
        check::<i8>();
        check::<i16>();
        check::<i32>();
        check::<i64>();
    }
    
    #[test]
    fn check_mon() {
        fn check<T>() where T: Mon, for<'a> &'a T: MonOps<T> {}
        check::<i8>();
        check::<i16>();
        check::<i32>();
        check::<i64>();
    }

    #[test]
    fn check_ring() {
        fn check<T>() where T: Ring, for<'a> &'a T: RingOps<T> {}
        check::<i8>();
        check::<i16>();
        check::<i32>();
        check::<i64>();
    }
}