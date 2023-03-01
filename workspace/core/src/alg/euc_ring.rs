// Euclidean Rings

use std::ops::{Div, DivAssign, Rem, RemAssign};
use crate::{Ring, RingOps};

pub trait EucRingOps<T = Self>: 
    RingOps<T> + 
    Div<T, Output = T> +
    for<'a> Div<&'a T, Output = T> +
    Rem<T, Output = T> +
    for<'a> Rem<&'a T, Output = T> +
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
        x * (y / g)
    }
}