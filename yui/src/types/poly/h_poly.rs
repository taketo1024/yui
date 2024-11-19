use std::fmt::Display;
use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Neg, DivAssign, RemAssign, Div, Rem};
use std::str::FromStr;
use num_traits::{Zero, One};
use auto_impl_ops::auto_ops;

use crate::poly::var::fmt_mono;
use crate::poly::Poly;
use crate::{AddGrp, AddGrpOps, AddMon, AddMonOps, Elem, EucRing, EucRingOps, Field, FieldOps, Mon, MonOps, Ring, RingOps};

use super::Mono;

// Homogeneous polynomial
#[derive(Clone, Copy, Debug, Default)]
pub struct HPoly<const X: char, R> { 
    deg: usize,
    coeff: R
}

impl<const X: char, R> HPoly<X, R> { 
    pub fn new(deg: usize, coeff: R) -> Self { 
        Self { deg, coeff }
    }

    pub fn coeff(&self) -> &R { 
        &self.coeff
    }

    pub fn deg(&self) -> usize { 
        self.deg
    }

    pub fn from_const(r: R) -> Self { 
        Self::new(0, r)
    }

    pub fn variable() -> Self
    where R: One { 
        Self::new(1, R::one())
    }
}

impl<const X: char, R> Display for HPoly<X, R>
where R: Display {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use crate::util::format::lc;
        let x = fmt_mono(&X.to_string(), &self.deg, true);
        let t = lc([(x, &self.coeff)].into_iter());
        t.fmt(f)
    }
}

impl<const X: char, R> FromStr for HPoly<X, R>
where R: Ring + FromStr, for<'x> &'x R: RingOps<R> {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let p = Poly::<X, R>::from_str(s)?;
        if p.is_zero() { 
            Ok(Self::zero())
        } else if p.nterms() == 1 { 
            let (x, a) = p.any_term().unwrap();
            let p = Self::new(x.deg(), a.clone());
            Ok(p)
        } else {
            Err(())
        }
    }
}

impl<const X: char, R> Zero for HPoly<X, R>
where R: AddMon, for<'x> &'x R: AddMonOps<R> {
    fn zero() -> Self {
        Self::new(0, R::zero())
    }

    fn is_zero(&self) -> bool {
        self.coeff.is_zero()
    }
}

impl<const X: char, R> One for HPoly<X, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn one() -> Self {
        Self::new(0, R::one())
    }

    fn is_one(&self) -> bool {
        self.deg == 0 && self.coeff.is_one()
    }
}

impl<const X: char, R> PartialEq for HPoly<X, R>
where R: PartialEq + Zero {
    fn eq(&self, other: &Self) -> bool {
        if self.coeff.is_zero() && other.coeff.is_zero() {
            true
        } else {
            self.deg == other.deg && self.coeff == other.coeff
        }
    }
}

impl<const X: char, R> Eq for HPoly<X, R>
where R: Eq + Zero {}

#[auto_ops]
impl<const X: char, R> AddAssign<&HPoly<X, R>> for HPoly<X, R>
where R: AddMon, for<'x> &'x R: AddMonOps<R> {
    fn add_assign(&mut self, rhs: &HPoly<X, R>) {
        if self.is_zero() { 
            *self = rhs.clone()
        } else if rhs.is_zero() { 
            return
        } else { 
            assert_eq!(self.deg, rhs.deg, "{self} + {rhs} is not homogeneous.");
            self.coeff.add_assign(&rhs.coeff)
        }
    }
}

#[auto_ops]
impl<const X: char, R> SubAssign<&HPoly<X, R>> for HPoly<X, R>
where R: AddGrp, for<'x> &'x R: AddGrpOps<R> {
    fn sub_assign(&mut self, rhs: &HPoly<X, R>) {
        if self.is_zero() { 
            *self = -rhs
        } else if rhs.is_zero() { 
            return
        } else { 
            assert_eq!(self.deg, rhs.deg, "{self} - {rhs} is not homogeneous.");
            self.coeff.sub_assign(&rhs.coeff)
        }
    }
}

impl<const X: char, R> Neg for HPoly<X, R>
where R: AddGrp, for<'x> &'x R: AddGrpOps<R> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self::new(self.deg, -self.coeff)
    }
}

impl<const X: char, R> Neg for &HPoly<X, R>
where R: AddGrp, for<'x> &'x R: AddGrpOps<R> {
    type Output = HPoly<X, R>;
    fn neg(self) -> Self::Output {
        HPoly::new(self.deg, -&self.coeff)
    }
}

#[auto_ops]
impl<const X: char, R> MulAssign<&R> for HPoly<X, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn mul_assign(&mut self, rhs: &R) {
        if rhs.is_one() { 
            return
        }
        self.coeff *= rhs
    }
}

#[auto_ops]
impl<const X: char, R> MulAssign<&HPoly<X, R>> for HPoly<X, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn mul_assign(&mut self, rhs: &HPoly<X, R>) {
        if rhs.is_one() { 
            return
        }
        self.deg += rhs.deg;
        self.coeff *= &rhs.coeff;
    }
}

macro_rules! impl_alg_op {
    ($trait:ident) => {
        impl<const X: char, R> $trait<Self> for HPoly<X, R>
        where R: Ring, for<'x> &'x R: RingOps<R> {}

        impl<const X: char, R> $trait<HPoly<X, R>> for &HPoly<X, R>
        where R: Ring, for<'x> &'x R: RingOps<R> {}
    };
}

impl_alg_op!(AddMonOps);
impl_alg_op!(AddGrpOps);
impl_alg_op!(MonOps);
impl_alg_op!(RingOps);

impl<const X: char, R> Elem for HPoly<X, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn math_symbol() -> String {
        format!("{}[{}]", R::math_symbol(), X)
    }
}

impl<const X: char, R> AddMon for HPoly<X, R>
where R: Ring, for<'x> &'x R: RingOps<R> {}

impl<const X: char, R> AddGrp for HPoly<X, R>
where R: Ring, for<'x> &'x R: RingOps<R> {}

impl<const X: char, R> Mon for HPoly<X, R>
where R: Ring, for<'x> &'x R: RingOps<R> {}

impl<const X: char, R> Ring for HPoly<X, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn inv(&self) -> Option<Self> {
        if self.deg > 0 { 
            return None
        }
        let a = self.coeff.inv()?;
        let inv = Self::new(0, a);
        Some(inv)
    }

    fn is_unit(&self) -> bool {
        self.deg == 0 && self.coeff.is_unit()
    }

    fn normalizing_unit(&self) -> Self {
        let u = self.coeff.normalizing_unit();
        Self::from_const(u)
    }

    fn c_weight(&self) -> f64 {
        self.coeff.c_weight()
    }
}

impl<const X: char, R> HPoly<X, R>
where R: Field, for<'x> &'x R: FieldOps<R> {
    pub fn div_rem(&self, rhs: &Self) -> (Self, Self) { 
        assert!(!rhs.is_zero());

        if self.deg < rhs.deg { 
            return (Self::zero(), self.clone())
        }
        
        let (i, a) = (self.deg, &self.coeff); // ax^i
        let (j, b) = ( rhs.deg,  &rhs.coeff); // bx^j
        
        let k = i - j; // >= 0
        let c = a / b;
        let q = HPoly::new(k, c); // cx^k = (a/b) x^{i-j}.
        let r = Self::zero();
        
        (q, r)
    }
}

#[auto_ops]
impl<const X: char, R> Div<&HPoly<X, R>> for HPoly<X, R>
where R: Field, for<'x> &'x R: FieldOps<R> {
    type Output = Self;

    fn div(self, rhs: &HPoly<X, R>) -> Self {
        self.div_rem(rhs).0
    }
}

#[auto_ops]
impl<const X: char, R> Rem<&HPoly<X, R>> for HPoly<X, R>
where R: Field, for<'x> &'x R: FieldOps<R> {
    type Output = Self;

    fn rem(self, rhs: &HPoly<X, R>) -> Self::Output {
        self.div_rem(rhs).1
    }
}

impl<const X: char, R> EucRingOps<HPoly<X, R>> for HPoly<X, R>
where R: Field, for<'x> &'x R: FieldOps<R> {}

impl<const X: char, R> EucRingOps<HPoly<X, R>> for &HPoly<X, R>
where R: Field, for<'x> &'x R: FieldOps<R> {}

impl<const X: char, R> EucRing for HPoly<X, R>
where R: Field, for<'x> &'x R: FieldOps<R> {}

#[cfg(feature = "tex")]
mod tex {
    use crate::tex::TeX;
    use super::*;

    impl<const X: char, R> TeX for HPoly<X, R>
    where R: Ring + TeX, for<'x> &'x R: RingOps<R> {
        fn tex_math_symbol() -> String { 
            format!("{}[{X}]", R::tex_math_symbol())
        }

        fn tex_string(&self) -> String {
            use crate::util::format::lc;
            let x = fmt_mono(&X.to_string(), &self.deg, false);
            lc([(x, self.coeff.tex_string())].into_iter())
        }
    }
}

#[cfg(test)]
mod tests { 
    use crate::Ratio;

    use super::*;

    #[test]
    fn zero() { 
        type R = i64;
        type P = HPoly<'x', R>;

        let a = P::new(0, R::zero());
        let b = P::new(0, R::one());
        let c = P::new(1, R::zero());

        assert!(a.is_zero());
        assert!(!b.is_zero());
        assert!(c.is_zero());
    }
    
    #[test]
    fn one() { 
        type R = i64;
        type P = HPoly<'x', R>;

        let a = P::new(0, R::zero());
        let b = P::new(0, R::one());
        let c = P::new(1, R::one());

        assert!(!a.is_one());
        assert!(b.is_one());
        assert!(!c.is_one());
    }
    
    #[test]
    fn eq() { 
        type R = i64;
        type P = HPoly<'x', R>;

        let a = P::new(0, R::zero());
        let b = P::new(0, R::one());
        let c = P::new(1, R::zero());
        let d = P::new(1, R::one());

        assert_ne!(a, b);
        assert_ne!(b, c);
        assert_ne!(c, d);
        assert_eq!(a, c);
    }

    #[test]
    fn add() { 
        type R = i64;
        type P = HPoly<'x', R>;

        let o = P::zero();
        let a = P::new(2, 1);
        let b = P::new(2, 3);

        assert_eq!(a + b, P::new(2, 4));
        assert_eq!(a + o, a);
        assert_eq!(o + a, a);
    }

    #[test]
    fn sub() { 
        type R = i64;
        type P = HPoly<'x', R>;

        let o = P::zero();
        let a = P::new(2, 1);
        let b = P::new(2, 3);

        assert_eq!(a - b, P::new(2, -2));
        assert_eq!(a - o, a);
        assert_eq!(o - a, -a);
    }

    #[test]
    fn mul() { 
        type R = i64;
        type P = HPoly<'x', R>;

        let o = P::zero();
        let a = P::new(2, 4);
        let b = P::new(3, 5);

        assert_eq!(a * b, P::new(5, 20));
        assert_eq!(a * o, o);
        assert_eq!(o * a, o);
    }

    #[test]
    fn div() { 
        type R = Ratio<i64>;
        type P = HPoly<'x', R>;

        let o = P::zero();
        let a = P::new(5, R::from(4));
        let b = P::new(3, R::from(5));

        assert_eq!(a / b, P::new(2, R::new(4, 5)));
        assert_eq!(b / a, o);
    }

    #[test]
    fn rem() { 
        type R = Ratio<i64>;
        type P = HPoly<'x', R>;

        let o = P::zero();
        let a = P::new(5, R::from(4));
        let b = P::new(3, R::from(5));

        assert_eq!(a % b, o);
        assert_eq!(b % a, b);
    }

    #[test]
    fn from_str() {
        type R = i64;
        type P = HPoly<'x', R>;

        assert_eq!(P::from_str("0"), Ok(P::zero()));
        assert_eq!(P::from_str("3"), Ok(P::from_const(3)));
        assert_eq!(P::from_str("x"), Ok(P::variable()));
        assert_eq!(P::from_str("x^2"), Ok(P::new(2, 1)));
        // assert_eq!(P::from_str("3x^2"), Ok(P::new(2, 3))); // not supported yet
        assert_eq!(P::from_str("x + 1"), Err(()));
    }

    #[test]
    #[cfg(feature = "tex")]
    fn tex() { 
        use crate::tex::TeX;
        type R = i64;
        type P = HPoly<'x', R>;

        assert_eq!(P::tex_math_symbol(), "\\mathbb{Z}[x]");

        assert_eq!(&P::zero().tex_string(), "0");
        assert_eq!(&P::from_const(3).tex_string(), "3");
        assert_eq!(&P::variable().tex_string(), "x");
        assert_eq!(&P::new(2, 1).tex_string(), "x^2");
        assert_eq!(&P::new(2, 3).tex_string(), "3x^2");
    }
}