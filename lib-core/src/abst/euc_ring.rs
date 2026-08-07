//! Euclidean ring (a.k.a. Euclidean domain): a [`Ring`] equipped with division `/`
//! and remainder `%`, supporting Euclidean algorithms (gcd, extended gcd, lcm).
//!
//! The name is "EucRing" here rather than the more standard "Euclidean domain"
//! to make the trait hierarchy `EucRing: Ring` immediately visible.
//!
//! See: <https://en.wikipedia.org/wiki/Euclidean_domain>

use std::ops::{Div, DivAssign, Rem, RemAssign};
use crate::abst::{Ring, RingOps};

/// Helper trait extending [`RingOps`] with `Div` and `Rem` reference variants
/// so [`EucRing`] can require them via one HRTB.
pub trait EucRingOps<T = Self>:
    RingOps<T> +
    Div<T, Output = T> +
    for<'a> Div<&'a T, Output = T> +
    Rem<T, Output = T> +
    for<'a> Rem<&'a T, Output = T> +
{}

/// A Euclidean ring: a [`Ring`] with division `/` and remainder `%`
/// satisfying the Euclidean property — i.e. for any `x, y` with `y ≠ 0`,
/// `x = (x / y) · y + (x % y)` with `x % y` strictly "smaller" than `y`.
///
/// See: <https://en.wikipedia.org/wiki/Euclidean_domain>
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
    /// `true` iff `self` divides `y` (and `self ≠ 0`).
    fn divides(&self, y: &Self) -> bool {
        !self.is_zero() && (y % self).is_zero()
    }

    /// Greatest common divisor, returned in normalized form.
    ///
    /// `gcd(0, 0) = 0`.
    ///
    /// See: <https://en.wikipedia.org/wiki/Euclidean_algorithm>
    fn gcd(x: &Self, y: &Self) -> Self {
        if x.is_zero() && y.is_zero() { return Self::zero() }
        if x.divides(y) { return x.normalized() }
        if y.divides(x) { return y.normalized() }

        let (mut x, mut y) = (x.clone(), y.clone());

        while !y.is_zero() {
            let r = &x % &y;
            (x, y) = (y, r);
        }

        x.into_normalized()
    }

    /// Extended gcd: returns `(d, s, t)` such that `d = gcd(x, y) = s·x + t·y`,
    /// with `d` normalized.
    ///
    /// See: <https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm>
    fn gcdx(x: &Self, y: &Self) -> (Self, Self, Self) {
        if x.is_zero() && y.is_zero() { return (Self::zero(), Self::zero(), Self::zero()) }

        // `d` is normalized, so the cofactor is the normalizing unit rather than `1`.
        if x.divides(y) {
            let u = x.normalizing_unit();
            return (x * &u, u, Self::zero())
        }
        if y.divides(x) {
            let u = y.normalizing_unit();
            return (y * &u, Self::zero(), u)
        }

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

        let (d, s, t) = (x, s0, t0);

        let u = d.normalizing_unit();
        match u.is_one() {
            true  => (d, s, t),
            false => (d * &u, s * &u, t * &u)
        }
    }

    /// Least common multiple, returned in normalized form.
    fn lcm(x: &Self, y: &Self) -> Self {
        if x.is_zero() || y.is_zero() { return Self::zero() }

        let g = Self::gcd(x, y);
        let m = x * (y / g);
        m.into_normalized()
    }
}