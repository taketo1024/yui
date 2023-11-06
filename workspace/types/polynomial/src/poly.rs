use std::fmt::{Display, Debug};
use std::iter::{Sum, Product};
use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Neg, DivAssign, RemAssign, Div, Rem};
use std::str::FromStr;
use num_traits::{Zero, One, Pow};
use auto_impl_ops::auto_ops;

use yui_core::{Elem, AddMon, AddMonOps, AddGrp, AddGrpOps, Mon, MonOps, Ring, RingOps, EucRing, EucRingOps, Field, FieldOps};
use yui_lin_comb::{LinComb, Gen};

use crate::{MultiDeg, Univar, BiVar, MultiVar};

// A polynomial is a linear combination of monomials over R.

// Univar-type (ordinary, Laurent)
pub type Poly  <const X: char, R> = PolyBase<Univar<X, usize>, R>;          
pub type LPoly <const X: char, R> = PolyBase<Univar<X, isize>, R>;

// Bivar-type (ordinary, Laurent)
pub type Poly2 <const X: char, const Y: char, R> = PolyBase<BiVar<X, Y, usize>, R>;
pub type LPoly2<const X: char, const Y: char, R> = PolyBase<BiVar<X, Y, isize>, R>;

// Multivar-type (ordinary, Laurent)
pub type PolyN <const X: char, R> = PolyBase<MultiVar<X, usize>, R>;
pub type LPolyN<const X: char, R> = PolyBase<MultiVar<X, isize>, R>;

pub trait Mono: 
    Mul<Output = Self> + 
    Div<Output = Self> + 
    One + 
    PartialOrd + 
    Ord + 
    From<Self::Deg> +
    Gen
{
    type Deg;

    fn deg(&self) -> Self::Deg;
    fn is_unit(&self) -> bool;
    fn inv(&self) -> Option<Self>;
    fn divides(&self, other: &Self) -> bool;
}

#[derive(Clone, PartialEq, Eq, Default)]
pub struct PolyBase<X, R>
where 
    X: Mono, 
    R: Ring, for<'x> &'x R: RingOps<R>
{
    data: LinComb<X, R>,
    zero: (X, R)
}

impl<X, R> PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn new(data: LinComb<X, R>) -> Self { 
        Self { data, zero: (X::one(), R::zero()) }
    }

    pub fn from_const(r: R) -> Self {
        Self::from((X::one(), r))
    }

    pub fn inner(&self) -> &LinComb<X, R> { 
        &self.data
    }

    pub fn len(&self) -> usize { 
        self.data.len()
    }

    pub fn coeff(&self, x: &X) -> &R {
        self.data.coeff(x)
    }

    pub fn coeff_for(&self, i: X::Deg) -> &R {
        self.data.coeff(&X::from(i))
    }

    pub fn iter(&self) -> impl Iterator<Item = (&X, &R)> {
        self.data.iter()
    }

    pub fn into_iter(self) -> impl Iterator<Item = (X, R)> {
        self.data.into_iter()
    }

    pub fn reduce(&mut self) { 
        self.data.reduce()
    }

    pub fn reduced(&self) -> Self { 
        Self::from(self.data.reduced())
    }

    pub fn is_const(&self) -> bool { 
        self.iter().all(|(i, r)| 
            i.is_one() || !i.is_one() && r.is_zero()
        )
    }

    pub fn is_mono(&self) -> bool { 
        self.iter().filter(|(_, r)| 
            !r.is_zero()
        ).count() <= 1
    }

    pub fn const_term(&self) -> &R { 
        self.coeff(&X::one())
    }

    pub fn lead_term(&self) -> (&X, &R) { 
        self.iter()
            .filter(|(_, r)| !r.is_zero())
            .max_by_key(|(i, _)| *i)
            .unwrap_or((&self.zero.0, &self.zero.1))
    }

    pub fn lead_coeff(&self) -> &R { 
        self.lead_term().1
    }

    pub fn lead_deg(&self) -> X::Deg { 
        self.lead_term().0.deg()
    }

    pub fn map_coeffs<R2, F>(&self, f: F) -> PolyBase<X, R2>
    where 
        R2: Ring, for<'x> &'x R2: RingOps<R2>, 
        F: Fn(&R) -> R2
    {
        PolyBase::<X, R2>::from( self.data.map_coeffs(f) )
    }
}

macro_rules! impl_univar {
    ($I:ty) => {
        impl<const X: char, R> PolyBase<Univar<X, $I>, R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            pub fn variable() -> Self { 
                Self::from(Self::mono(1))
            }

            pub fn mono(i: $I) -> Univar<X, $I> {
                Univar::from(i)
            }
        }
    };
}

impl_univar!(usize);
impl_univar!(isize);

macro_rules! impl_bivar {
    ($I:ty) => {
        impl<const X: char, const Y: char, R> PolyBase<BiVar<X, Y, $I>, R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            pub fn variable(i: usize) -> Self { 
                assert!(i < 2);
                let d = if i == 0 { (1, 0) } else { (0, 1) };
                Self::from(BiVar::from(d))
            }

            pub fn mono(i: $I, j: $I) -> BiVar<X, Y, $I> {
                BiVar::from((i, j))
            }
        }
    };
}

impl_bivar!(usize);
impl_bivar!(isize);

macro_rules! impl_mvar {
    ($I:ty) => {
        impl<const X: char, R> PolyBase<MultiVar<X, $I>, R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            pub fn variable(i: usize) -> Self { 
                let d = MultiDeg::from((i, 1));
                Self::from(Univar::from(d)) // x^1
            }

            pub fn mono<const N: usize>(degs: [$I; N]) -> MultiVar<X, $I> {
                MultiVar::from(degs)
            }
        }
    };
}

impl_mvar!(usize);
impl_mvar!(isize);

macro_rules! impl_polyn_funcs {
    ($I:ty) => {
        impl<const X: char, R> PolyBase<MultiVar<X, $I>, R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            pub fn lead_term_for(&self, k: usize) -> (&MultiVar<X, $I>, &R) { 
                self.iter()
                    .filter(|(x, r)| !r.is_zero() && x.deg_for(k) > 0)
                    .max_by(|(x, _), (y, _)|
                        Ord::cmp(
                            &x.deg_for(k), &y.deg_for(k)
                        ).then(Ord::cmp(
                            &x, &y
                        ))
                    )
                    .unwrap_or(self.lead_term())
            }

            pub fn lead_coeff_for(&self, k: usize) -> &R { 
                self.lead_term_for(k).1
            }

            pub fn lead_deg_for(&self, k: usize) -> MultiDeg<$I> { 
                self.lead_term_for(k).0.deg()
            }
        }
    };
}

impl_polyn_funcs!(usize);
impl_polyn_funcs!(isize);

impl<X, R> From<X> for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(x: X) -> Self {
        Self::from((x, R::one()))
    }
}

impl<X, R> From<(X, R)> for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(pair: (X, R)) -> Self {
        let t = LinComb::from(pair);
        Self::from(t)
    }
}

impl<X, R> FromIterator<(X, R)> for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from_iter<T: IntoIterator<Item = (X, R)>>(iter: T) -> Self {
        Self::from(LinComb::from_iter(iter))
    }
}

impl<X, R> From<i32> for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(i: i32) -> Self {
        Self::from_const(R::from(i))
    }
}

impl<X, R> From<LinComb<X, R>> for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(data: LinComb<X, R>) -> Self {
        Self::new(data)
    }
}

impl<X, R> From<PolyBase<X, R>> for LinComb<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(poly: PolyBase<X, R>) -> Self {
        poly.data
    }
}

impl<X, R> FromStr for PolyBase<X, R>
where X: Mono + FromStr, R: Ring + FromStr, for<'x> &'x R: RingOps<R> {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Ok(r) = R::from_str(s) { 
            Ok(Self::from_const(r))
        } else if let Ok(x) = X::from_str(s) { 
            Ok(Self::from(x))
        } else {
            // TODO support more complex format.
            Err(())
        }
    }
}

impl<X, R> Display for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data.fmt(f, false) // descending order
    }
}

impl<X, R> Debug for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data.fmt(f, false) // descending order
    }
}

impl<X, R> Zero for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn zero() -> Self {
        Self::from(LinComb::zero())
    }

    fn is_zero(&self) -> bool {
        self.is_const() && self.const_term().is_zero()
    }
}

impl<X, R> One for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn one() -> Self {
        Self::from((X::one(), R::one()))
    }

    fn is_one(&self) -> bool {
        self.is_const() && self.const_term().is_one()
    }
}

impl<X, R> Neg for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self::from(-self.data)
    }
}

impl<X, R> Neg for &PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = PolyBase<X, R>;
    fn neg(self) -> Self::Output {
        PolyBase::from(-&self.data)
    }
}

macro_rules! impl_assop {
    ($trait:ident, $method:ident) => {
        #[auto_ops]
        impl<X, R> $trait<&PolyBase<X, R>> for PolyBase<X, R>
        where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
            fn $method(&mut self, rhs: &PolyBase<X, R>) {
                self.data.$method(&rhs.data)
            }
        }
    };
}

impl_assop!(AddAssign, add_assign);
impl_assop!(SubAssign, sub_assign);

#[auto_ops]
impl<X, R> MulAssign<&R> for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn mul_assign(&mut self, rhs: &R) {
        self.data *= rhs
    }
}

#[auto_ops]
impl<X, R> MulAssign<&PolyBase<X, R>> for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn mul_assign(&mut self, rhs: &PolyBase<X, R>) {
        if rhs.is_one() {
            // do nothing
        } else if rhs.is_const() { 
            *self *= rhs.const_term()
        } else if self.is_const() { 
            *self = rhs * self.const_term()
        } else { 
            self.data *= &rhs.data
        }
    }
}

macro_rules! impl_accum {
    ($trait:ident, $method:ident, $accum_trait:ident, $accum_method:ident, $accum_init:ident) => {
        impl<X, R> $trait for PolyBase<X, R>
        where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
            fn $method<Iter: Iterator<Item = Self>>(iter: Iter) -> Self {
                let mut res = Self::$accum_init();
                for r in iter { Self::$accum_method(&mut res, r) }
                return res;
            }
        }

        impl<'a, X, R> $trait<&'a Self> for PolyBase<X, R>
        where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
            fn $method<Iter: Iterator<Item = &'a Self>>(iter: Iter) -> Self {
                let mut res = Self::$accum_init();
                for r in iter { Self::$accum_method(&mut res, r) }
                return res;
            }
        }
    };
}

impl_accum!(Sum, sum, AddAssign, add_assign, zero);
impl_accum!(Product, product, MulAssign, mul_assign, one);

macro_rules! impl_pow_unsigned {
    ($t:ty) => {
        impl<X, R> Pow<$t> for &PolyBase<X, R>
        where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
            type Output = PolyBase<X, R>;
            fn pow(self, n: $t) -> Self::Output {
                let mut res = PolyBase::one();
                for _ in 0..n { 
                    res *= self
                }
                res
            }
        }
    };
}

impl_pow_unsigned!(u32);
impl_pow_unsigned!(u64);
impl_pow_unsigned!(usize);

macro_rules! impl_pow_signed {
    ($t:ty) => {
        impl<X, R> Pow<$t> for &PolyBase<X, R>
        where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
            type Output = PolyBase<X, R>;
            fn pow(self, n: $t) -> Self::Output {
                if n >= 0 { 
                    self.pow(n as usize)
                } else {
                    let inv = self.inv().unwrap();
                    (&inv).pow(-n as usize)
                }
            }
        }
    }
}

impl_pow_signed!(i32);
impl_pow_signed!(i64);
impl_pow_signed!(isize);

// TODO support eval for other cases

macro_rules! impl_eval_bivar {
    ($I:ty) => {
        impl<const X: char, const Y: char, R> PolyBase<BiVar<X, Y, $I>, R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            pub fn eval(&self, x: &R, y: &R) -> R
            where for<'x, 'y> &'x R: Pow<&'y $I, Output = R> { 
                self.iter().map(|(i, r)| { 
                    r * i.eval(x, y)
                }).sum()
            }
        }
    };
}

impl_eval_bivar!(usize);
impl_eval_bivar!(isize);

macro_rules! impl_alg_op {
    ($trait:ident) => {
        impl<X, R> $trait<Self> for PolyBase<X, R>
        where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {}

        impl<X, R> $trait<PolyBase<X, R>> for &PolyBase<X, R>
        where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {}
    };
}

impl_alg_op!(AddMonOps);
impl_alg_op!(AddGrpOps);
impl_alg_op!(MonOps);
impl_alg_op!(RingOps);

impl<X, R> Elem for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn math_symbol() -> String {
        format!("{}[{}]", R::math_symbol(), X::math_symbol())
    }
}

impl<X, R> AddMon for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {}

impl<X, R> AddGrp for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {}

impl<X, R> Mon for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {}

impl<X, R> Ring for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn inv(&self) -> Option<Self> {
        if !self.is_mono() { 
            return None
        }

        let (x, a) = self.lead_term(); // (a x^i)^{-1} = a^{-1} x^{-i}
        if let (Some(xinv), Some(ainv)) = (x.inv(), a.inv()) { 
            Some(Self::from((xinv, ainv)))
        } else { 
            None
        }
    }

    fn is_unit(&self) -> bool {
        if self.is_mono() { 
            let (x, a) = self.lead_term();
            x.is_unit() && a.is_unit()
        } else { 
            false
        }
    }

    fn normalizing_unit(&self) -> Self {
        let u = self.lead_coeff().normalizing_unit();
        Self::from_const(u)
    }
}

// UPoly over R: Field

impl<const X: char, R> Poly<X, R>
where R: Field, for<'x> &'x R: FieldOps<R> {
    pub fn div_rem(&self, rhs: &Self) -> (Self, Self) { 
        let iter = |f: Self, g: &Self| -> (Self, Self) { 
            if f.lead_deg() < g.lead_deg() { 
                return (Self::zero(), f)
            }

            let (i, a) = f.lead_term(); // ax^i
            let (j, b) = g.lead_term(); // bx^j
            
            let k = i.deg() - j.deg(); // >= 0
            let c = a / b;
            let x = Univar::from(k);
            let q = Poly::from((x, c));   // cx^k = (a/b) x^{i-j}.
            let r = f - &q * g;
            
            (q, r)
        };

        let mut q = Self::zero();
        let mut r = self.clone();

        let i = self.lead_deg();
        let j =  rhs.lead_deg();
        
        for _ in j ..= i { // passes if j > i. 
            let (q1, r1) = iter(r, rhs);
            q += q1;
            r = r1;
        }

        r.reduce();

        (q, r)
    }
}

#[auto_ops]
impl<const X: char, R> Div<&Poly<X, R>> for Poly<X, R>
where R: Field, for<'x> &'x R: FieldOps<R> {
    type Output = Self;

    fn div(self, rhs: &Poly<X, R>) -> Self {
        self.div_rem(rhs).0
    }
}

#[auto_ops]
impl<const X: char, R> Rem<&Poly<X, R>> for Poly<X, R>
where R: Field, for<'x> &'x R: FieldOps<R> {
    type Output = Self;

    fn rem(self, rhs: &Poly<X, R>) -> Self::Output {
        self.div_rem(rhs).1
    }
}

impl<const X: char, R> EucRingOps<Poly<X, R>> for Poly<X, R>
where R: Field, for<'x> &'x R: FieldOps<R> {}

impl<const X: char, R> EucRingOps<Poly<X, R>> for &Poly<X, R>
where R: Field, for<'x> &'x R: FieldOps<R> {}

impl<const X: char, R> EucRing for Poly<X, R>
where R: Field, for<'x> &'x R: FieldOps<R> {}

#[cfg(test)]
mod tests {
    use yui_utils::map;
    use super::*;

    #[test]
    fn init() { 
        type P = Poly::<'x', i32>; 

        let x = P::mono;
        let f = P::from_iter([(x(0), 1), (x(1), 2), (x(2), -3)]);
        assert_eq!(f.data.coeff(&x(0)), &1);
        assert_eq!(f.data.coeff(&x(1)), &2);
        assert_eq!(f.data.coeff(&x(2)), &-3);
        assert_eq!(f.data.coeff(&x(3)), &0);
    }
 
    #[test]
    fn display_poly() { 
        type P = Poly::<'x', i32>; 

        let x = P::mono;
        let f = P::from_iter([(x(0), 1), (x(1), 2), (x(2), -3)]);
        assert_eq!(&f.to_string(), "-3x² + 2x + 1");
    }
 
    #[test]
    fn display_lpoly() { 
        type P = LPoly::<'x', i32>; 

        let x = P::mono;
        let f = P::from_iter([(x(-1), 4), (x(0), 2), (x(2), 3)]);
        assert_eq!(&f.to_string(), "3x² + 2 + 4x⁻¹");
    }

    #[test]
    fn display_poly2() { 
        type P = Poly2::<'x', 'y', i32>; 

        let xy = |i, j| P::mono(i, j);
        let f = P::from_iter([(xy(0, 0), 3), (xy(1, 0), 2), (xy(2, 3), 3)]);
        assert_eq!(&f.to_string(), "3x²y³ + 2x + 3");
    }
 
    #[test]
    fn display_mpoly() { 
        type P = PolyN::<'x', i32>; 

        let xn = |i| P::mono(i);
        let f = P::from_iter([
            (xn([0,0,0]),  3),
            (xn([1,0,0]), -1),
            (xn([3,0,2]),  2),
        ]);
        assert_eq!(&f.to_string(), "2x₀³x₂² - x₀ + 3");
    }

    #[test]
    fn display_mlpoly() { 
        type P = LPolyN::<'x', i32>; 

        let xn = |i| P::mono(i);
        let f = P::from_iter([
            (xn([ 0,0,0]), 3),
            (xn([ 1,0,0]), 1),
            (xn([-3,1,3]), 2),
        ]);
        assert_eq!(&f.to_string(), "x₀ + 2x₀⁻³x₁x₂³ + 3");
    }

    #[test]
    fn reduce() { 
        type P = Poly::<'x', i32>;

        let x = P::mono;
        let data = map!{
            x(0) => 0, 
            x(1) => 1, 
            x(2) => 0 
        };
        let f = P::from(LinComb::new(data));
        assert_eq!(f.len(), 3);

        let f = f.reduced();
        assert_eq!(f.len(), 1);
    }

    #[test]
    fn zero() {
        type P = Poly::<'x', i32>;

        let zero = P::zero();
        assert_eq!(zero, P::from_iter([]));

        let x = P::mono;
        let data = map!{
            x(0) => 0 
        };
        let f = P::from(LinComb::new(data));

        assert!(f.is_zero());
    }

    #[test]
    fn one() {
        type P = Poly::<'x', i32>;

        let x = P::mono;
        let p = P::one();
        assert_eq!(p, P::from_iter([(x(0), 1)]));
    }

    #[test]
    fn variable() {
        type P = Poly::<'x', i32>;

        let x = P::mono;
        let p = P::variable();
        assert_eq!(p, P::from_iter([(x(1), 1)]));
    }

    #[test]
    fn variable_bivar() {
        type P = Poly2::<'x', 'y', i32>;

        let xy = |i, j| P::mono(i, j);
        let p = P::variable(0);
        let q = P::variable(1);
        assert_eq!(p, P::from_iter([(xy(1, 0), 1)]));
        assert_eq!(q, P::from_iter([(xy(0, 1), 1)]));
    }

    #[test]
    fn variable_mvar() {
        type P = PolyN::<'x', i32>;

        let xn = |i| P::mono(i);
        let p = P::variable(0);
        let q = P::variable(1);
        assert_eq!(p, P::from_iter([(xn([1, 0]), 1)]));
        assert_eq!(q, P::from_iter([(xn([0, 1]), 1)]));
    }

    #[test]
    fn coeff() { 
        type P = Poly::<'x', i32>;

        let x = P::mono;
        let f = P::from_iter([(x(0), 2), (x(1), 3), (x(2), -4)]);
        
        assert_eq!(f.coeff_for(0), &2);
        assert_eq!(f.coeff_for(1), &3);
        assert_eq!(f.coeff_for(2), &-4);
        assert_eq!(f.coeff_for(3), &0);
    }

    #[test]
    fn const_term() { 
        type P = Poly::<'x', i32>;
        let x = P::mono;
        let f = P::from_iter([(x(0), 2), (x(1), 3), (x(2), -4)]);
        assert_eq!(f.const_term(), &2);

        let f = P::from_iter([(x(1), 3), (x(2), -4)]);
        assert_eq!(f.const_term(), &0);
    }

    #[test]
    fn lead_term() { 
        type P = Poly::<'x', i32>;

        let x = P::mono;
        let f = P::from_iter([(x(0), 2), (x(1), 3), (x(2), -4)]);
        let (m, a) = f.lead_term();

        assert_eq!(m.deg(), 2);
        assert_eq!(a, &-4);

        let f = P::zero();
        let (m, a) = f.lead_term();
        assert_eq!(m.deg(), 0);
        assert_eq!(a, &0);
    }

    #[test]
    fn add() { 
        type P = Poly::<'x', i32>;

        let x = P::mono;
        let f = P::from_iter([(x(0), 2), (x(1), 3), (x(2), -4)]);
        let g = P::from_iter([(x(0), -3), (x(1), -3), (x(3), 5)]);

        assert_eq!(f + g, P::from_iter([(x(0), -1), (x(1), 0), (x(2), -4), (x(3), 5)]));
    }

    #[test]
    fn neg() { 
        type P = Poly::<'x', i32>;

        let x = P::mono;
        let f = P::from_iter([(x(0), 2), (x(1), 3), (x(2), -4)]);

        assert_eq!(-f, P::from_iter([(x(0), -2), (x(1), -3), (x(2), 4)]));
    }

    #[test]
    fn sub() { 
        type P = Poly::<'x', i32>;
        let x = P::mono;
        let f = P::from_iter([(x(0), 2), (x(1), 3), (x(2), -4)]);
        let g = P::from_iter([(x(0), -3), (x(1), -3), (x(3), 5)]);
        assert_eq!(f - g, P::from_iter([(x(0), 5), (x(1), 6), (x(2), -4), (x(3), -5)]));
    }

    #[test]
    fn mul() { 
        type P = Poly::<'x', i32>;

        let x = P::mono;
        let f = P::from_iter([(x(0), 2), (x(1), 3), (x(2), -4)]);
        let g = P::from_iter([(x(0), -3), (x(1), -3), (x(3), 5)]);

        assert_eq!(f * g, P::from_iter([(x(0), -6), (x(1), -15), (x(2), 3), (x(3), 22), (x(4), 15), (x(5), -20)]));
    }

    #[test]
    fn mul_const() { 
        type P = Poly::<'x', i32>;

        let x = P::mono;
        let f = P::from_iter([(x(0), 2), (x(1), 3), (x(2), -4)]);
        let g = P::from_const(3);

        assert_eq!(&f * &g, P::from_iter([(x(0), 6), (x(1), 9), (x(2), -12)]));
        assert_eq!(&g * &f, P::from_iter([(x(0), 6), (x(1), 9), (x(2), -12)]));
    }

    #[test]
    fn pow() { 
        type P = Poly::<'x', i32>;

        let x = P::mono;
        let f = P::from_iter([(x(0), 3), (x(1), 2)]);

        assert_eq!(f.pow(0), P::one());
        assert_eq!(f.pow(1), f);
        assert_eq!(f.pow(0), P::one());
        assert_eq!(f.pow(2), P::from_iter([(x(0), 9), (x(1), 12), (x(2), 4)]));
    }

    #[test]
    fn pow_laurent() { 
        type P = LPoly::<'x', i32>;

        let x = P::mono;
        let f = P::variable();

        assert_eq!(f.pow(0), P::one());
        assert_eq!(f.pow(-1), P::from((x(-1), 1)));
        assert_eq!(f.pow(0), P::one());
        assert_eq!(f.pow(-2), P::from((x(-2), 1)));
    }

    #[test]
    fn inv() { 
        type P = Poly::<'x', i32>;

        let x = P::mono;
        let f = P::from_const(1);
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from_const(1)));

        let f = P::from_const(0);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::from_const(2);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::variable();
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::from_iter([(x(0), 1), (x(1), 1)]);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn inv_rat() { 
        use yui_ratio::Ratio;
        type R = Ratio<i32>;
        type P = Poly::<'x', R>;

        let x = P::mono;
        let f = P::from_const(R::from_numer(1));

        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from_const(R::from_numer(1))));

        let f = P::from_const(R::zero());
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::from_const(R::from_numer(2));
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from_const(R::new(1, 2))));

        let f = P::variable();
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::from_iter([(x(0), R::one()), (x(1), R::one())]);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn inv_laurent() { 
        type P = LPoly::<'x', i32>;

        let x = P::mono;
        let f = P::from_const(1);
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from_const(1)));

        let f = P::from_const(0);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::from_const(2);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::variable();
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from((x(-1), 1))));

        let f = P::from((x(1), 2));
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::from_iter([(x(0), 1), (x(1), 1)]);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn inv_laurent_rat() { 
        use yui_ratio::Ratio;
        type R = Ratio<i32>;
        type P = LPoly::<'x', R>;

        let x = P::mono;
        let f = P::from_const(R::from_numer(1));
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from_const(R::from_numer(1))));

        let f = P::from_const(R::zero());
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::from_const(R::from_numer(2));
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from_const(R::new(1, 2))));

        let f = P::variable();
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from((x(-1), R::one()))));

        let f = P::from((x(1), R::from_numer(2)));
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from((x(-1), R::new(1, 2)))));

        let f = P::from_iter([(x(0), R::one()), (x(1), R::one())]);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn div_rem() { 
        use yui_ratio::Ratio;
        type R = Ratio<i32>;
        type P = Poly::<'x', R>;

        let x = P::mono;
        let f = P::from_iter([(x(0), R::from_numer(1)), (x(1), R::from_numer(2)), (x(2), R::from_numer(1))]);
        let g = P::from_iter([(x(0), R::from_numer(3)), (x(1), R::from_numer(2))]);
        let (q, r) = f.div_rem(&g);

        assert_eq!(q, P::from_iter([(x(0), R::new(1, 4)), (x(1), R::new(1, 2))]));
        assert_eq!(r, P::from_const( R::new(1, 4)) );
        assert_eq!(f, q * &g + r);

        let (q, r) = g.div_rem(&f);
        assert_eq!(q, P::zero());
        assert_eq!(r, g);
    }

    #[test]
    fn from_str() { 
        type P = Poly::<'x', i32>;

        assert_eq!(P::from_str("-3"), Ok(P::from_const(-3)));
        assert_eq!(P::from_str("x"), Ok(P::variable()));
        assert_eq!(P::from_str("y"), Err(()));

        // TODO support more complex types
    }

    #[test]
    fn eval_bivar() { 
        type P = Poly2::<'x', 'y', i32>;

        let xy = |i, j| P::mono(i, j);
        let p = P::from_iter([(xy(0,0),3), (xy(1,0),2), (xy(0,1),-1), (xy(1,1),4)]);
        let v = p.eval(&2, &3); // 3 + 2(2) - 1(3) + 4(2*3)
        assert_eq!(v, 28);
    }

    #[test]
    fn lead_term_for() {
        type P = PolyN::<'x', i32>;

        let xn = |i| P::mono(i);
        let f = P::from_iter([
            (xn([1,2,3]), 1),
            (xn([2,1,3]), 2),
            (xn([0,2,4]), 3),
            (xn([5,5,5]), 0), // should be ignored
        ]);
        assert_eq!(f.lead_term(),      (&xn([2,1,3]), &2));
        assert_eq!(f.lead_term_for(0), (&xn([2,1,3]), &2));
        assert_eq!(f.lead_term_for(1), (&xn([1,2,3]), &1));
        assert_eq!(f.lead_term_for(2), (&xn([0,2,4]), &3));
        assert_eq!(f.lead_term_for(3), f.lead_term());
    }
}