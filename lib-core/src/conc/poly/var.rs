//! Univariate monomial `X^d`. With `I = isize` it is a Laurent monomial
//! (allowing `d < 0`).
//!
//! See: <https://en.wikipedia.org/wiki/Monomial>,
//! <https://en.wikipedia.org/wiki/Laurent_polynomial>

use std::fmt::{Display, Debug};
use std::ops::{AddAssign, Mul, MulAssign, DivAssign, SubAssign, Div};
use std::str::FromStr;
use num_traits::{Zero, One, ToPrimitive, Pow, FromPrimitive};
use auto_impl_ops::auto_ops;

use crate::abst::{MathType, IndexType};
use crate::lc::LcKey;
use crate::util::format::superscript;
use crate::util::parse_err::ParseErr;
use super::{Mono, MonoOrd};

/// A univariate monomial `X^d`, where the variable symbol `X` is a const
/// generic and `I` is the exponent type (`usize` or `isize`).
#[derive(Clone, Copy, PartialEq, Eq, Hash, Default, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde_with::DeserializeFromStr))]
pub struct Var<const X: char, I>(
    I
);

impl<const X: char, I> Var<X, I> {
    pub fn var_symbol() -> char { 
        X
    }

    pub fn eval<R>(&self, x: &R) -> R
    where R: Mul<Output = R>, I: Copy, for<'x> &'x R: Pow<I, Output = R> {
        x.pow(self.0)
    }

    fn to_string_u(&self, unicode: bool) -> String
    where I: ToPrimitive { 
        fmt_mono(&X.to_string(), &self.0, unicode)
    }
}

impl<const X: char, I> FromStr for Var<X, I>
where I: Zero + FromStr + FromPrimitive {
    type Err = ParseErr;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == "1" {
            Ok(Self(I::zero()))
        } else if let Some(d) = parse_mono_deg(&X.to_string(), s) {
            Ok(Self(d))
        } else {
            Err(ParseErr::invalid(s, &format!("a monomial in {X}")))
        }
    }
}

impl<const X: char, I> From<I> for Var<X, I> {
    fn from(d: I) -> Self {
        Self(d)
    }
}

impl<const X: char, I> One for Var<X, I>
where I: for<'x >AddAssign<&'x I> + Zero {
    fn one() -> Self {
        Self::from(I::zero()) // x^0 = 1.
    }
}

// note: clippy flags the `+`/`-` below, but exponents add on × and subtract on ÷.
// `#[allow]` doesn't help: `#[auto_ops]` regenerates these impls and drops it.
#[auto_ops]
impl<const X: char, I> MulAssign<&Var<X, I>> for Var<X, I>
where I: for<'x >AddAssign<&'x I> {
    fn mul_assign(&mut self, rhs: &Var<X, I>) {
        self.0 += &rhs.0 // x^i * x^j = x^{i+j}
    }
}

#[auto_ops]
impl<const X: char, I> DivAssign<&Var<X, I>> for Var<X, I>
where I: for<'x >SubAssign<&'x I> {
    fn div_assign(&mut self, rhs: &Var<X, I>) {
        self.0 -= &rhs.0 // x^i / x^j = x^{i-j}
    }
}

impl<const X: char, I> Display for Var<X, I>
where I: ToPrimitive { 
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = self.to_string_u(true);
        f.write_str(&s)
    }
}

impl<const X: char, I> Debug for Var<X, I>
where I: ToPrimitive { 
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Display::fmt(self, f)
    }
}

#[cfg(feature = "serde")]
impl<const X: char, I> serde::Serialize for Var<X, I>
where I: ToPrimitive {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where S: serde::Serializer {
        serializer.serialize_str(&self.to_string_u(false))
    }
}

impl<const X: char, I> MathType for Var<X, I>
where I: IndexType + ToPrimitive { 
    fn math_symbol() -> String {
        format!("{X}")
    }
}

impl<const X: char, I> MonoOrd for Var<X, I> 
where I: IndexType + ToPrimitive {
    fn cmp_lex(&self, other: &Self) -> std::cmp::Ordering {
        I::cmp(&self.0, &other.0)
    }

    fn cmp_grlex(&self, other: &Self) -> std::cmp::Ordering {
        I::cmp(&self.0, &other.0)
    }
}

impl<const X: char, I> LcKey for Var<X, I> 
where I: IndexType + ToPrimitive {}

macro_rules! impl_univar_unsigned {
    ($I:ty) => {
        impl<const X: char> Mono for Var<X, $I> {
            type Deg = $I;

            fn deg(&self) -> Self::Deg {
                self.0
            }

            fn is_unit(&self) -> bool { 
                self.0.is_zero()
            }

            fn inv(&self) -> Option<Self> {
                if self.is_unit() { 
                    Some(Self(0))
                } else { 
                    None
                }
            }

            fn divides(&self, other: &Self) -> bool { 
                self.0 <= other.0
            }
        }
    };
}

macro_rules! impl_univar_signed {
    ($I:ty) => {
        impl<const X: char> Mono for Var<X, $I> {
            type Deg = $I;

            fn deg(&self) -> Self::Deg {
                self.0
            }

            fn is_unit(&self) -> bool { 
                true
            }

            fn inv(&self) -> Option<Self> { // (x^i)^{-1} = x^{-i}
                Some(Self(-self.0))
            }

            fn divides(&self, _other: &Self) -> bool { 
                true
            }
        }
    };
}

impl_univar_unsigned!(usize);
impl_univar_signed!  (isize);

mod tex {
    use crate::util::tex::TeX;
    use super::*;

    impl<const X: char, I> TeX for Var<X, I>
    where I: ToPrimitive {
        fn tex_math_symbol() -> String { 
            String::from(X)
        }
        fn tex_string(&self) -> String {
            self.to_string_u(false)
        }
    }
}

pub(crate) fn fmt_mono<I>(x: &str, d: &I, unicode: bool) -> String
where I: ToPrimitive {
    let d = d.to_isize().unwrap();
    if d.is_zero() { 
        "1".to_string()
    } else if d.is_one() { 
        x.to_string()
    } else if unicode {
        let e = superscript(d); 
        format!("{x}{e}")
    } else { 
        let d = d.to_string();
        if d.len() == 1 { 
            format!("{x}^{d}")
        } else { 
            format!("{x}^{{{d}}}")
        }
    }
}

pub(crate) fn parse_mono_deg<I>(x: &str, s: &str) -> Option<I>
where I: FromStr + FromPrimitive {
    use regex::Regex;

    let p1 = format!(r"^{x}\^(-?[0-9]+)$");
    let p2 = format!(r"^{x}\^\{{(-?[0-9]+)\}}$");

    let r1 = Regex::new(&p1).unwrap();
    let r2 = Regex::new(&p2).unwrap();

    if s == "1" { 
        I::from_i32(0)
    } else if s == x { 
        I::from_i32(1)
    } else if let Some(c) = r1.captures(s) { 
        I::from_str(&c[1]).ok()
    } else if let Some(c) = r2.captures(s) { 
        I::from_str(&c[1]).ok()
    } else { 
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn init() { 
        type M = Var<'X',usize>;
        let x = M::from;

        let d = x(2);

        assert_eq!(d.0, 2);
        assert_eq!(d.deg(), 2);
        assert_eq!(M::var_symbol(), 'X');
    }

    #[test]
    fn from_str() { 
        type M = Var<'X',isize>;
        let x = M::from;
        
        assert_eq!(M::from_str("1"), Ok(M::one()));
        assert_eq!(M::from_str("X"), Ok(x(1)));
        assert_eq!(M::from_str("X^2"), Ok(x(2)));
        assert_eq!(M::from_str("X^{-2}"), Ok(x(-2)));
        assert!(M::from_str("2").is_err());
        assert!(M::from_str("x").is_err());

        // unbraced, beyond one digit and negative
        assert_eq!(M::from_str("X^23"), Ok(x(23)));
        assert_eq!(M::from_str("X^-2"), Ok(x(-2)));
    }

    #[test]
    fn display() { 
        type M = Var<'X', isize>;
        let x = M::from;

        let d = x(0);
        assert_eq!(&d.to_string(), "1");

        let d = x(1);
        assert_eq!(&d.to_string(), "X");

        let d = x(2);
        assert_eq!(&d.to_string(), "X²");

        let d = x(-1);
        assert_eq!(&d.to_string(), "X⁻¹");

        let d = x(-2);
        assert_eq!(&d.to_string(), "X⁻²");
    }

    #[test]
    fn neg_opt_unsigned() { 
        type M = Var<'X',usize>;
        let x = M::from;

        let d = x(0);
        assert_eq!(d.inv(), Some(x(0)));

        let d = x(1);
        assert_eq!(d.inv(), None);
    }

    #[test]
    fn neg_opt_signed() { 
        type M = Var<'X',isize>;
        let x = M::from;

        let d = x(0);
        assert_eq!(d.inv(), Some(x(0)));

        let d = x(1);
        assert_eq!(d.inv(), Some(x(-1)));

        let d = x(-3);
        assert_eq!(d.inv(), Some(x(3)));
    }

    #[test]
    fn eval() { 
        type M = Var<'X', usize>;
        let x = M::from;

        let d = x(0);
        assert_eq!(d.eval(&2), 1);

        let d = x(1);
        assert_eq!(d.eval(&2), 2);

        let d = x(2);
        assert_eq!(d.eval(&2), 4);
    }

    #[test]
    fn ord() { 
        type M = Var<'X', isize>;
        let x = M::from;

        assert!(x(0) < x(1));
        assert!(x(1) < x(2));
        assert!(x(-1) < x(0));
    }

    #[test]
    #[cfg(feature = "serde")]
    fn serialize() { 
        type M = Var<'X', isize>;
        let x = M::from;

        let d = x(0);
        let ser = serde_json::to_string(&d).unwrap();
        let des = serde_json::from_str::<M>(&ser).unwrap();
        
        assert_eq!(&ser, "\"1\"");
        assert_eq!(d, des);

        let d = x(1);
        let ser = serde_json::to_string(&d).unwrap();
        let des = serde_json::from_str::<M>(&ser).unwrap();

        assert_eq!(&ser, "\"X\"");
        assert_eq!(d, des);

        let d = x(2);
        let ser = serde_json::to_string(&d).unwrap();
        let des = serde_json::from_str::<M>(&ser).unwrap();
        
        assert_eq!(&ser, "\"X^2\"");
        assert_eq!(d, des);

        let d = x(21);
        let ser = serde_json::to_string(&d).unwrap();
        let des = serde_json::from_str::<M>(&ser).unwrap();
        
        assert_eq!(&ser, "\"X^{21}\"");
        assert_eq!(d, des);

        let d = x(-2);
        let ser = serde_json::to_string(&d).unwrap();
        let des = serde_json::from_str::<M>(&ser).unwrap();
        
        assert_eq!(&ser, "\"X^{-2}\"");
        assert_eq!(d, des);
    }
}