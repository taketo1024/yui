use std::fmt::{Display, Debug};
use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Neg, DivAssign, RemAssign, Div, Rem};
use std::str::FromStr;
use delegate::delegate;
use num_traits::{Zero, One, Pow};
use auto_impl_ops::auto_ops;

use crate::{Elem, AddMon, AddMonOps, AddGrp, AddGrpOps, Mon, MonOps, Ring, RingOps, EucRing, EucRingOps, Field, FieldOps};
use crate::lc::Lc;
use super::{MultiDeg, Var, Var2, Var3,MultiVar, Mono, MonoOrd};

// A polynomial is a linear combination of monomials over R.

// Var-type (ordinary, Laurent)
pub type Poly  <const X: char, R> = PolyBase<Var<X, usize>, R>;          
pub type LPoly <const X: char, R> = PolyBase<Var<X, isize>, R>;

// Bivar-type (ordinary, Laurent)
pub type Poly2 <const X: char, const Y: char, R> = PolyBase<Var2<X, Y, usize>, R>;
pub type LPoly2<const X: char, const Y: char, R> = PolyBase<Var2<X, Y, isize>, R>;

// Trivar-type (ordinary, Laurent)
pub type Poly3 <const X: char, const Y: char, const Z: char, R> = PolyBase<Var3<X, Y, Z, usize>, R>;
pub type LPoly3<const X: char, const Y: char, const Z: char, R> = PolyBase<Var3<X, Y, Z, isize>, R>;

// Multivar-type (ordinary, Laurent)
pub type PolyN <const X: char, R> = PolyBase<MultiVar<X, usize>, R>;
pub type LPolyN<const X: char, R> = PolyBase<MultiVar<X, isize>, R>;

#[derive(Clone, PartialEq, Eq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(transparent))]
pub struct PolyBase<X, R>
where 
    X: Mono, 
    R: Ring, for<'x> &'x R: RingOps<R>
{
    data: Lc<X, R>,
    #[cfg_attr(feature = "serde", serde(skip))]
    zero: (X, R)
}

impl<X, R> PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn new(data: Lc<X, R>) -> Self { 
        Self { data, zero: (X::one(), R::zero()) }
    }

    pub fn from_const(r: R) -> Self {
        Self::from((X::one(), r))
    }

    pub fn inner(&self) -> &Lc<X, R> { 
        &self.data
    }

    delegate! { 
        to self.data {
            pub fn nterms(&self) -> usize;
            pub fn any_term(&self) -> Option<(&X, &R)>;
            #[call(is_gen)] pub fn is_mono(&self) -> bool;
            #[call(as_gen)] pub fn as_mono(&self) -> Option<X>;
            pub fn coeff(&self, x: &X) -> &R;
            pub fn iter(&self) -> impl Iterator<Item = (&X, &R)>;
        }
    }

    pub fn coeff_for(&self, i: X::Deg) -> &R {
        self.data.coeff(&X::from(i))
    }

    pub fn is_const(&self) -> bool { 
        self.iter().all(|(x, _)| x.is_one())
    }

    pub fn const_term(&self) -> &R { 
        self.coeff(&X::one())
    }

    pub fn lead_term(&self) -> (&X, &R) { 
        self.iter().max_by(|t1, t2| MonoOrd::cmp_grlex(t1.0, t2.0))
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

    pub fn sort_terms_by<F>(&self, cmp: F) -> impl Iterator<Item = (&X, &R)>
    where F: Fn(&X, &X) -> std::cmp::Ordering { 
        self.data.sort_terms_by(cmp)
    }

    pub fn to_string_by<F>(&self, cmp: F, descending: bool) -> String
    where F: Fn(&X, &X) -> std::cmp::Ordering {
        self.data.to_string_by(cmp, descending)
    }
}

macro_rules! impl_var_specific {
    ($I:ty) => {
        // Univar
        impl<const X: char, R> PolyBase<Var<X, $I>, R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            pub fn variable() -> Self { 
                Self::from(Self::mono(1))
            }

            pub fn mono(i: $I) -> Var<X, $I> {
                Var::from(i)
            }

            pub fn eval(&self, x: &R) -> R
            where for<'x, 'y> &'x R: Pow<&'y $I, Output = R> { 
                R::sum(self.iter().map(|(i, r)| { 
                    r * i.eval(x)
                }))
            }
        }

        // Bivar
        impl<const X: char, const Y: char, R> PolyBase<Var2<X, Y, $I>, R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            pub fn variable(i: usize) -> Self { 
                assert!(i < 2);
                let d = if i == 0 { (1, 0) } else { (0, 1) };
                Self::from(Var2::from(d))
            }

            pub fn mono(i: $I, j: $I) -> Var2<X, Y, $I> {
                Var2::from((i, j))
            }

            pub fn eval(&self, x: &R, y: &R) -> R
            where for<'x, 'y> &'x R: Pow<&'y $I, Output = R> { 
                R::sum(self.iter().map(|(i, r)| { 
                    r * i.eval(x, y)
                }))
            }
        }

        // Trivar
        impl<const X: char, const Y: char, const Z: char, R> PolyBase<Var3<X, Y, Z, $I>, R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            pub fn variable(i: usize) -> Self { 
                assert!(i < 3);
                let d = match i { 
                    0 => (1, 0, 0),
                    1 => (0, 1, 0),
                    2 => (0, 0, 1),
                    _ => panic!()
                };
                Self::from(Var3::from(d))
            }

            pub fn mono(i: $I, j: $I, k: $I) -> Var3<X, Y, Z, $I> {
                Var3::from((i, j, k))
            }

            pub fn eval(&self, x: &R, y: &R, z: &R) -> R
            where for<'x, 'y> &'x R: Pow<&'y $I, Output = R> { 
                R::sum(self.iter().map(|(i, r)| { 
                    r * i.eval(x, y, z)
                }))
            }
        }

        // MultiVar
        impl<const X: char, R> PolyBase<MultiVar<X, $I>, R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            pub fn variable(i: usize) -> Self { 
                let d = MultiDeg::from((i, 1));
                Self::from(MultiVar::from(d)) // x^1
            }

            pub fn mono<const N: usize>(degs: [$I; N]) -> MultiVar<X, $I> {
                MultiVar::from(degs)
            }

            pub fn lead_term_for(&self, k: usize) -> Option<(&MultiVar<X, $I>, &R)> { 
                self.iter()
                    .filter(|(x, _)| x.deg_for(k) > 0)
                    .max_by(|(x, _), (y, _)|
                        Ord::cmp( &x.deg_for(k), &y.deg_for(k))
                        .then_with(|| 
                            MultiVar::cmp_grlex(&x, &y)
                        )
                    )
            }
        }
    };
}

impl_var_specific!(usize);
impl_var_specific!(isize);

impl<X, R> From<X> for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(x: X) -> Self {
        Self::from((x, R::one()))
    }
}

impl<X, R> From<(X, R)> for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(pair: (X, R)) -> Self {
        let t = Lc::from(pair);
        Self::from(t)
    }
}

impl<X, R> FromIterator<(X, R)> for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from_iter<T: IntoIterator<Item = (X, R)>>(iter: T) -> Self {
        Self::from(Lc::from_iter(iter))
    }
}

impl<X, R> From<Lc<X, R>> for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(data: Lc<X, R>) -> Self {
        Self::new(data)
    }
}

impl<X, R> From<PolyBase<X, R>> for Lc<X, R>
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

impl<X, R> IntoIterator for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    type Item = (X, R);
    type IntoIter = std::collections::hash_map::IntoIter<X, R>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

impl<X, R> Display for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.to_string_by(X::cmp_for_display, true))
    }
}

impl<X, R> Debug for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Display::fmt(self, f)
    }
}

impl<X, R> Zero for PolyBase<X, R>
where X: Mono, R: Ring, for<'x> &'x R: RingOps<R> {
    fn zero() -> Self {
        Self::from(Lc::zero())
    }

    fn is_zero(&self) -> bool {
        self.data.is_zero()
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
        if self.nterms() != 1 { 
            return None
        }

        let (x, a) = self.any_term()?; // (a x^i)^{-1} = a^{-1} x^{-i}
        let (xinv, ainv) = (x.inv()?, a.inv()?);
        let inv = Self::from((xinv, ainv));
        Some(inv)
    }

    fn is_unit(&self) -> bool {
        if self.nterms() == 1 { 
            let (x, a) = self.any_term().unwrap();
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
            let x = Var::from(k);
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

cfg_if::cfg_if! { 
    if #[cfg(feature = "tex")] {
        use crate::TeX;

        impl<X, R> TeX for PolyBase<X, R>
        where X: Mono + TeX, R: Ring + TeX, for<'x> &'x R: RingOps<R> {
            fn to_tex_string(&self) -> String {
                use crate::util::format::lc;
                let terms = self.sort_terms_by(|x, y| X::cmp_lex(x, y).reverse()).map(|(x, r)|
                    (x.to_tex_string(), r.to_tex_string())
                );
                lc(terms)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::Ratio;
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

        let xy = P::mono;
        let f = P::from_iter([(xy(0, 0), 3), (xy(1, 0), 2), (xy(2, 3), 3)]);
        assert_eq!(&f.to_string(), "3x²y³ + 2x + 3");
    }
 
    #[test]
    fn display_mpoly() { 
        type P = PolyN::<'x', i32>; 

        let xn = P::mono;
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

        let xn = P::mono;
        let f = P::from_iter([
            (xn([ 0,0,0]), 3),
            (xn([ 1,0,0]), 1),
            (xn([-3,1,3]), 2),
        ]);
        assert_eq!(&f.to_string(), "x₀ + 3 + 2x₀⁻³x₁x₂³");
    }

    #[test]
    fn zero() {
        type P = Poly::<'x', i32>;

        let z = P::zero();
        assert_eq!(z, P::from_iter([]));
        assert!(z.is_zero());
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

        let xy = P::mono;
        let p = P::variable(0);
        let q = P::variable(1);
        assert_eq!(p, P::from_iter([(xy(1, 0), 1)]));
        assert_eq!(q, P::from_iter([(xy(0, 1), 1)]));
    }

    #[test]
    fn variable_mvar() {
        type P = PolyN::<'x', i32>;

        let xn = P::mono;
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
        assert!(f.is_unit());
        assert_eq!(f.inv(), Some(P::from_const(1)));

        let f = P::from_const(0);
        assert!(!f.is_unit());
        assert_eq!(f.inv(), None);

        let f = P::from_const(2);
        assert!(!f.is_unit());
        assert_eq!(f.inv(), None);

        let f = P::variable();
        assert!(!f.is_unit());
        assert_eq!(f.inv(), None);

        let f = P::from_iter([(x(0), 1), (x(1), 1)]);
        assert!(!f.is_unit());
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn inv_rat() { 
        type R = Ratio<i32>;
        type P = Poly::<'x', R>;

        let x = P::mono;
        let f = P::from_const(R::from(1));

        assert!(f.is_unit());
        assert_eq!(f.inv(), Some(P::from_const(R::from(1))));

        let f = P::from_const(R::zero());
        assert!(!f.is_unit());
        assert_eq!(f.inv(), None);

        let f = P::from_const(R::from(2));
        assert!(f.is_unit());
        assert_eq!(f.inv(), Some(P::from_const(R::new(1, 2))));

        let f = P::variable();
        assert!(!f.is_unit());
        assert_eq!(f.inv(), None);

        let f = P::from_iter([(x(0), R::one()), (x(1), R::one())]);
        assert!(!f.is_unit());
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn inv_laurent() { 
        type P = LPoly::<'x', i32>;

        let x = P::mono;
        let f = P::from_const(1);
        assert!(f.is_unit());
        assert_eq!(f.inv(), Some(P::from_const(1)));

        let f = P::from_const(0);
        assert!(!f.is_unit());
        assert_eq!(f.inv(), None);

        let f = P::from_const(2);
        assert!(!f.is_unit());
        assert_eq!(f.inv(), None);

        let f = P::variable();
        assert!(f.is_unit());
        assert_eq!(f.inv(), Some(P::from((x(-1), 1))));

        let f = P::from((x(1), 2));
        assert!(!f.is_unit());
        assert_eq!(f.inv(), None);

        let f = P::from_iter([(x(0), 1), (x(1), 1)]);
        assert!(!f.is_unit());
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn inv_laurent_rat() { 
        type R = Ratio<i32>;
        type P = LPoly::<'x', R>;

        let x = P::mono;
        let f = P::from_const(R::from(1));
        assert!(f.is_unit());
        assert_eq!(f.inv(), Some(P::from_const(R::from(1))));

        let f = P::from_const(R::zero());
        assert!(!f.is_unit());
        assert_eq!(f.inv(), None);

        let f = P::from_const(R::from(2));
        assert!(f.is_unit());
        assert_eq!(f.inv(), Some(P::from_const(R::new(1, 2))));

        let f = P::variable();
        assert!(f.is_unit());
        assert_eq!(f.inv(), Some(P::from((x(-1), R::one()))));

        let f = P::from((x(1), R::from(2)));
        assert!(f.is_unit());
        assert_eq!(f.inv(), Some(P::from((x(-1), R::new(1, 2)))));

        let f = P::from_iter([(x(0), R::one()), (x(1), R::one())]);
        assert!(!f.is_unit());
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn div_rem() { 
        type R = Ratio<i32>;
        type P = Poly::<'x', R>;

        let x = P::mono;
        let f = P::from_iter([(x(0), R::from(1)), (x(1), R::from(2)), (x(2), R::from(1))]);
        let g = P::from_iter([(x(0), R::from(3)), (x(1), R::from(2))]);
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

        let xy = P::mono;
        let p = P::from_iter([(xy(0,0),3), (xy(1,0),2), (xy(0,1),-1), (xy(1,1),4)]);
        let v = p.eval(&2, &3); // 3 + 2(2) - 1(3) + 4(2*3)
        assert_eq!(v, 28);
    }

    #[test]
    fn lead_term_for() {
        type P = PolyN::<'x', i32>;

        let xn = P::mono;
        let f = P::from_iter([
            (xn([1,2,3]), 1),
            (xn([2,1,3]), 2),
            (xn([0,2,4]), 3),
            (xn([5,5,5]), 0), // should be ignored
        ]);
        assert_eq!(f.lead_term_for(0), Some((&xn([2,1,3]), &2)));
        assert_eq!(f.lead_term_for(1), Some((&xn([1,2,3]), &1)));
        assert_eq!(f.lead_term_for(2), Some((&xn([0,2,4]), &3)));
        assert_eq!(f.lead_term_for(3), None);
    }

    #[test]
    #[cfg(feature = "serde")]
    fn serialize_univar() { 
        type P = Poly::<'x', i32>; 

        let x = P::mono;
        let f = P::from_iter([(x(0), 1), (x(1), 2), (x(2), -3)]);

        let ser = serde_json::to_string(&f).unwrap();
        let des = serde_json::from_str(&ser).unwrap();
        assert_eq!(f, des);
    }

    #[test]
    #[cfg(feature = "serde")]
    fn serialize_bivar() { 
        type P = LPoly2::<'x', 'y', i32>; 

        let xy = P::mono;
        let f = P::from_iter([
            (xy(0, 0), 3), 
            (xy(1, 0), 1), 
            (xy(-2, 13), -3)
        ]);

        let ser = serde_json::to_string(&f).unwrap();
        dbg!(&ser);
        let des = serde_json::from_str(&ser).unwrap();
        assert_eq!(f, des);
    }

    #[test]
    #[cfg(feature = "serde")]
    fn serialize_mvar() { 
        type P = LPolyN::<'x', i32>; 

        let xn = P::mono;
        let f = P::from_iter([
            (xn([0,0,0]),  3),
            (xn([1,0,0]), -1),
            (xn([12,0,-2]),  12),
        ]);

        let ser = serde_json::to_string(&f).unwrap();
        let des = serde_json::from_str(&ser).unwrap();
        assert_eq!(f, des);
    }

    #[test]
    #[cfg(feature = "tex")]
    fn tex() { 
        type P = LPolyN::<'x', i32>; 

        let xn = P::mono;
        let f = P::from_iter([
            (xn([ 0,0,0]),   3),
            (xn([ 3,0,1]),  -1),
            (xn([-2,1,3]), -2),
        ]);

        assert_eq!(f.to_tex_string(), "-x_0^3x_2 + 3 - 2x_0^{-2}x_1x_2^3");
    }
}