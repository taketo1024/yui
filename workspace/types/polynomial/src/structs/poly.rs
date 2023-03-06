use std::fmt::Display;
use std::iter::{Sum, Product};
use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Neg, DivAssign, RemAssign, Div, Rem};
use std::str::FromStr;
use itertools::Itertools;
use num_traits::{Zero, One, Pow};
use auto_impl_ops::auto_ops;

use yui_core::{Elem, AddMon, AddMonOps, AddGrp, AddGrpOps, Mon, MonOps, Ring, RingOps, EucRing, EucRingOps, Field, FieldOps};
use yui_lin_comb::LinComb;

use crate::{PolyGen, MDegree, Poly};

#[derive(Clone, PartialEq, Eq, Default, Debug)]
pub struct PolyBase<I, R>
where 
    I: PolyGen, 
    R: Ring, for<'x> &'x R: RingOps<R>
{
    data: LinComb<I, R>,
    zero: (I, R)
}

impl<I, R> PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new_raw(data: LinComb<I, R>) -> Self { 
        Self { data, zero: (I::one(), R::zero()) }
    }

    pub fn new(data: LinComb<I, R>) -> Self { 
        let mut new = Self::new_raw(data);
        new.reduce();
        new
    }

    pub fn from_deg(data: Vec<(I::Degree, R)>) -> Self {
        let data = data.into_iter().map(|(i, r)| 
            (I::from(i), r)
        ).collect_vec();
        Self::from(data)
    }

    pub fn mono(i: I::Degree, coeff: R) -> Self { 
        Self::from((I::from(i), coeff))
    }

    pub fn as_lincomb(&self) -> &LinComb<I, R> { 
        &self.data
    }

    pub fn len(&self) -> usize { 
        self.data.len()
    }

    pub fn coeff(&self, i: &I) -> &R {
        self.data.coeff(i)
    }

    pub fn coeff_for(&self, i: I::Degree) -> &R {
        self.data.coeff(&I::from(i))
    }

    pub fn iter(&self) -> impl Iterator<Item = (&I, &R)> {
        self.data.iter()
    }

    pub fn into_iter(self) -> impl Iterator<Item = (I, R)> {
        self.data.into_iter()
    }

    pub fn reduce(&mut self) { 
        self.data.reduce()
    }

    pub fn reduced(&self) -> Self { 
        Self::new_raw(self.data.reduced())
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
        self.coeff(&I::one())
    }

    pub fn lead_term(&self) -> (&I, &R) { 
        self.iter()
            .filter(|(_, r)| !r.is_zero())
            .max_by_key(|(i, _)| *i)
            .unwrap_or((&self.zero.0, &self.zero.1))
    }

    pub fn lead_coeff(&self) -> &R { 
        self.lead_term().1
    }

    pub fn lead_deg(&self) -> I::Degree { 
        self.lead_term().0.degree()
    }
}

impl<I, R> PolyBase<I, R>
where I: PolyGen, I::Degree: One, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn variable() -> Self { 
        Self::mono(I::Degree::one(), R::one()) // x^1
    }
}

impl<I, R> From<R> for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(a: R) -> Self {
        Self::new(LinComb::from((I::one(), a)))
    }
}

impl<I, R> From<(I, R)> for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(data: (I, R)) -> Self {
        Self::new(LinComb::from(data))
    }
}

impl<I, R> From<Vec<(I, R)>> for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(data: Vec<(I, R)>) -> Self {
        Self::new(LinComb::from(data))
    }
}

// for bivar-poly
impl<I, J, R> PolyBase<I, R>
where I: PolyGen<Degree = MDegree<J>>, J: Zero, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_deg2(data: Vec<((J, J), R)>) -> Self {
        let data = data.into_iter().map(|(d, r)| {
            let (d0, d1) = d;
            (vec![d0, d1], r)
        }).collect_vec();
        Self::from_mdeg(data)
    }
}

// for multivar-poly
impl<I, J, R> PolyBase<I, R>
where I: PolyGen<Degree = MDegree<J>>, J: Zero, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_mdeg(data: Vec<(Vec<J>, R)>) -> Self {
        let data = data.into_iter().map(|(v, r)| {
            let mdeg = MDegree::from_vec(v);
            (I::from(mdeg), r)
        }).collect_vec();
        Self::from(data)
    }
}

impl<I, R> FromStr for PolyBase<I, R>
where I: PolyGen + FromStr, R: Ring + FromStr, for<'x> &'x R: RingOps<R> {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Ok(r) = R::from_str(s) { 
            Ok(Self::from(r))
        } else if let Ok(i) = I::from_str(s) { 
            Ok(Self::from((i, R::one())))
        } else {
            // TODO support more complex format.
            Err(())
        }
    }
}

impl<I, R> Display for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data.fmt(f)
    }
}

impl<I, R> Zero for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn zero() -> Self {
        Self::new_raw(LinComb::zero())
    }

    fn is_zero(&self) -> bool {
        self.is_const() && self.const_term().is_zero()
    }
}

impl<I, R> One for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn one() -> Self {
        Self::new_raw(LinComb::from((I::one(), R::one())))
    }

    fn is_one(&self) -> bool {
        self.is_const() && self.const_term().is_one()
    }
}

impl<I, R> Neg for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self::new(-self.data)
    }
}

impl<I, R> Neg for &PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = PolyBase<I, R>;
    fn neg(self) -> Self::Output {
        PolyBase::new(-&self.data)
    }
}

macro_rules! impl_assop {
    ($trait:ident, $method:ident) => {
        #[auto_ops]
        impl<I, R> $trait<&PolyBase<I, R>> for PolyBase<I, R>
        where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
            fn $method(&mut self, rhs: &PolyBase<I, R>) {
                self.data.$method(&rhs.data)
            }
        }
    };
}

impl_assop!(AddAssign, add_assign);
impl_assop!(SubAssign, sub_assign);

#[auto_ops]
impl<I, R> MulAssign<&R> for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn mul_assign(&mut self, rhs: &R) {
        self.data *= rhs
    }
}

#[auto_ops]
impl<I, R> MulAssign<&PolyBase<I, R>> for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn mul_assign(&mut self, rhs: &PolyBase<I, R>) {
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
        impl<I, R> $trait for PolyBase<I, R>
        where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
            fn $method<Iter: Iterator<Item = Self>>(iter: Iter) -> Self {
                let mut res = Self::$accum_init();
                for r in iter { Self::$accum_method(&mut res, r) }
                return res;
            }
        }

        impl<'a, I, R> $trait<&'a Self> for PolyBase<I, R>
        where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
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
        impl<I, R> Pow<$t> for &PolyBase<I, R>
        where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
            type Output = PolyBase<I, R>;
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
        impl<I, R> Pow<$t> for &PolyBase<I, R>
        where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
            type Output = PolyBase<I, R>;
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

macro_rules! impl_alg_op {
    ($trait:ident) => {
        impl<I, R> $trait<Self> for PolyBase<I, R>
        where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {}

        impl<I, R> $trait<PolyBase<I, R>> for &PolyBase<I, R>
        where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {}
    };
}

impl_alg_op!(AddMonOps);
impl_alg_op!(AddGrpOps);
impl_alg_op!(MonOps);
impl_alg_op!(RingOps);

impl<I, R> Elem for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn set_symbol() -> String {
        format!("{}[{}]", R::set_symbol(), I::set_symbol())
    }
}

impl<I, R> AddMon for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {}

impl<I, R> AddGrp for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {}

impl<I, R> Mon for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {}

impl<I, R> Ring for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
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
        Self::from(self.lead_coeff().normalizing_unit())
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
            
            let k = i.degree() - j.degree(); // >= 0
            let c = a / b;
            let q = Poly::mono(k, c);   // cx^k = (a/b) x^{i-j}.
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
    use crate::{Mono, Poly2, PolyN, LPoly, LPolyN};
 
    #[test]
    fn init_poly() { 
        type P = Poly::<'x', i32>; 
        let f = P::from_deg(vec![(0, 3), (1, 2), (2, 3)]);
        assert_eq!(&f.to_string(), "3 + 2x + 3x²");
    }
 
    #[test]
    fn init_lpoly() { 
        type P = LPoly::<'x', i32>; 
        let f = P::from_deg(vec![(-1, 4), (0, 2), (2, 3)]);
        assert_eq!(&f.to_string(), "4x⁻¹ + 2 + 3x²");
    }

    #[test]
    fn init_poly2() { 
        type P = Poly2::<'x', 'y', i32>; 
        let f = P::from_deg2(vec![((0, 0), 3), ((1, 0), 2), ((2, 3), 3)]);
        assert_eq!(&f.to_string(), "3 + 2x + 3x²y³");
    }
 
    #[test]
    fn init_mpoly() { 
        type P = PolyN::<'x', i32>; 
        let f = P::from_mdeg(vec![
            (vec![], 3),
            (vec![1], -1),
            (vec![3, 0, 2], 2),
        ]);
        assert_eq!(&f.to_string(), "3 - x₀ + 2x₀³x₂²");
    }

    #[test]
    fn init_mlpoly() { 
        type P = LPolyN::<'x', i32>; 
        let f = P::from_mdeg(vec![
            (vec![], 3),
            (vec![1], 1),
            (vec![-3, 1, 3], 2),
        ]);
        assert_eq!(&f.to_string(), "3 + 2x₀⁻³x₁x₂³ + x₀");
    }

    #[test]
    fn reduce() { 
        type P = Poly::<'x', i32>;
        type X = Mono<'x', usize>;

        let data = map!{
            X::from(0) => 0, 
            X::from(1) => 1, 
            X::from(2) => 0 
        };
        let mut f = P::new_raw(LinComb::new_raw(data));
        assert_eq!(f.len(), 3);

        f.reduce();
        assert_eq!(f.len(), 1);
    }

    #[test]
    fn zero() {
        type P = Poly::<'x', i32>;
        let x = P::zero();
        assert_eq!(x, P::from_deg(vec![]));
    }

    #[test]
    fn one() {
        type P = Poly::<'x', i32>;
        let x = P::one();
        assert_eq!(x, P::from_deg(vec![(0, 1)]));
    }

    #[test]
    fn variable() {
        type P = Poly::<'x', i32>;
        let x = P::variable();
        assert_eq!(x, P::from_deg(vec![(1, 1)]));
    }

    #[test]
    fn coeff() { 
        type P = Poly::<'x', i32>;
        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        
        assert_eq!(f.coeff_for(0), &2);
        assert_eq!(f.coeff_for(1), &3);
        assert_eq!(f.coeff_for(2), &-4);
        assert_eq!(f.coeff_for(3), &0);
    }

    #[test]
    fn const_term() { 
        type P = Poly::<'x', i32>;
        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        assert_eq!(f.const_term(), &2);

        let f = P::from_deg(vec![(1, 3), (2, -4)]);
        assert_eq!(f.const_term(), &0);
    }

    #[test]
    fn lead_term() { 
        type P = Poly::<'x', i32>;

        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        let (x, a) = f.lead_term();
        assert_eq!(x.degree(), 2);
        assert_eq!(a, &-4);

        let f = P::zero();
        let (x, a) = f.lead_term();
        assert_eq!(x.degree(), 0);
        assert_eq!(a, &0);
    }

    #[test]
    fn add() { 
        type P = Poly::<'x', i32>;

        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        let g = P::from_deg(vec![(0, -3), (1, -3), (3, 5)]);
        assert_eq!(f + g, P::from_deg(vec![(0, -1), (1, 0), (2, -4), (3, 5)]));
    }

    #[test]
    fn neg() { 
        type P = Poly::<'x', i32>;

        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        assert_eq!(-f, P::from_deg(vec![(0, -2), (1, -3), (2, 4)]));
    }

    #[test]
    fn sub() { 
        type P = Poly::<'x', i32>;

        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        let g = P::from_deg(vec![(0, -3), (1, -3), (3, 5)]);
        assert_eq!(f - g, P::from_deg(vec![(0, 5), (1, 6), (2, -4), (3, -5)]));
    }

    #[test]
    fn mul() { 
        type P = Poly::<'x', i32>;

        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        let g = P::from_deg(vec![(0, -3), (1, -3), (3, 5)]);
        assert_eq!(f * g, P::from_deg(vec![(0, -6), (1, -15), (2, 3), (3, 22), (4, 15), (5, -20)]));
    }

    #[test]
    fn mul_const() { 
        type P = Poly::<'x', i32>;

        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        let g = P::from(3);

        assert_eq!(&f * &g, P::from_deg(vec![(0, 6), (1, 9), (2, -12)]));
        assert_eq!(&g * &f, P::from_deg(vec![(0, 6), (1, 9), (2, -12)]));
    }

    #[test]
    fn pow() { 
        type P = Poly::<'x', i32>;

        let f = P::from_deg(vec![(0, 3), (1, 2)]);
        assert_eq!(f.pow(0), P::one());
        assert_eq!(f.pow(1), f);       assert_eq!(f.pow(0), P::one());
        assert_eq!(f.pow(2), P::from_deg(vec![(0, 9), (1, 12), (2, 4)]));
    }

    #[test]
    fn pow_laurent() { 
        type P = LPoly::<'x', i32>;

        let f = P::variable();
        assert_eq!(f.pow(0), P::one());
        assert_eq!(f.pow(-1), P::mono(-1, 1));       assert_eq!(f.pow(0), P::one());
        assert_eq!(f.pow(-2), P::mono(-2, 1));
    }

    #[test]
    fn inv() { 
        type P = Poly::<'x', i32>;

        let f = P::from(1);
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from(1)));

        let f = P::from(0);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::from(2);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::variable();
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::from_deg(vec![(0, 1), (1, 1)]);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn inv_rat() { 
        use yui_ratio::Ratio;
        type R = Ratio<i32>;
        type P = Poly::<'x', R>;

        let f = P::from(R::from(1));
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from(R::from(1))));

        let f = P::from(R::zero());
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::from(R::from(2));
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from(R::new(1, 2))));

        let f = P::variable();
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::from_deg(vec![(0, R::one()), (1, R::one())]);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn inv_laurent() { 
        type P = LPoly::<'x', i32>;

        let f = P::from(1);
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from(1)));

        let f = P::from(0);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::from(2);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::variable();
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::mono(-1, 1)));

        let f = P::mono(1, 2);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::from_deg(vec![(0, 1), (1, 1)]);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn inv_laurent_rat() { 
        use yui_ratio::Ratio;
        type R = Ratio<i32>;
        type P = LPoly::<'x', R>;

        let f = P::from(R::from(1));
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from(R::from(1))));

        let f = P::from(R::zero());
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);

        let f = P::from(R::from(2));
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::from(R::new(1, 2))));

        let f = P::variable();
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::mono(-1, R::one())));

        let f = P::mono(1, R::from(2));
        assert_eq!(f.is_unit(), true);
        assert_eq!(f.inv(), Some(P::mono(-1, R::new(1, 2))));

        let f = P::from_deg(vec![(0, R::one()), (1, R::one())]);
        assert_eq!(f.is_unit(), false);
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn div_rem() { 
        use yui_ratio::Ratio;
        type R = Ratio<i32>;
        type P = Poly::<'x', R>;

        let f = P::from_deg(vec![(0, R::from(1)), (1, R::from(2)), (2, R::from(1))]);
        let g = P::from_deg(vec![(0, R::from(3)), (1, R::from(2))]);
        let (q, r) = f.div_rem(&g);

        assert_eq!(q, P::from_deg( vec![(0, R::new(1, 4)), (1, R::new(1, 2))]) );
        assert_eq!(r, P::from( R::new(1, 4)) );
        assert_eq!(f, q * &g + r);

        let (q, r) = g.div_rem(&f);
        assert_eq!(q, P::zero());
        assert_eq!(r, g);
    }

    #[test]
    fn from_str() { 
        type P = Poly::<'x', i32>;

        assert_eq!(P::from_str("-3"), Ok(P::from(-3)));
        assert_eq!(P::from_str("x"), Ok(P::variable()));
        assert_eq!(P::from_str("y"), Err(()));

        // TODO support more complex types
    }
}