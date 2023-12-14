use std::fmt::{Display, Debug};
use std::hash::Hash;
use std::ops::{AddAssign, Mul, MulAssign, DivAssign, SubAssign, Div};
use std::str::FromStr;
use num_traits::{Zero, One, ToPrimitive, Pow, FromPrimitive};
use auto_impl_ops::auto_ops;

use crate::Elem;
use crate::lc::Gen;
use crate::util::format::superscript;
use super::Mono;

// `Univar<X, I>` : represents monomials X^d (univar) or ΠX_i^{d_i} (multivar).
// `I` is one of `usize`, `isize`, `MultiDeg<usize>`, `MultiDeg<isize>`.

#[derive(Clone, PartialEq, Eq, Hash, Default, PartialOrd, Ord)]
pub struct Univar<const X: char, I>(
    pub(crate) I
);

impl<const X: char, I> Univar<X, I> {
    pub fn var_symbol() -> char { 
        X
    }
}

impl<const X: char, I> From<I> for Univar<X, I> {
    fn from(d: I) -> Self {
        Self(d)
    }
}

impl<const X: char, I> One for Univar<X, I>
where I: for<'x >AddAssign<&'x I> + Zero {
    fn one() -> Self {
        Self::from(I::zero()) // x^0 = 1.
    }
}

#[auto_ops]
impl<const X: char, I> MulAssign<&Univar<X, I>> for Univar<X, I>
where I: for<'x >AddAssign<&'x I> {
    fn mul_assign(&mut self, rhs: &Univar<X, I>) {
        self.0 += &rhs.0 // x^i * x^j = x^{i+j}
    }
}

#[auto_ops]
impl<const X: char, I> DivAssign<&Univar<X, I>> for Univar<X, I>
where I: for<'x >SubAssign<&'x I> {
    fn div_assign(&mut self, rhs: &Univar<X, I>) {
        self.0 -= &rhs.0 // x^i * x^j = x^{i+j}
    }
}

macro_rules! impl_univar {
    ($I:ty) => {
        impl<const X: char> Univar<X, $I> {
            pub fn eval<R>(&self, x: &R) -> R
            where R: Mul<Output = R>, for<'x, 'y> &'x R: Pow<&'y $I, Output = R> {
                x.pow(&self.0)
            }

            fn fmt_impl(&self, unicode: bool) -> String { 
                fmt_mono(&X.to_string(), self.0, unicode)
            }
        }
        
        impl<const X: char> Display for Univar<X, $I> { 
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let s = self.fmt_impl(true);
                f.write_str(&s)
            }
        }

        impl<const X: char> Debug for Univar<X, $I> { 
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                Display::fmt(self, f)
            }
        }
        
        impl<const X: char> FromStr for Univar<X, $I> {
            type Err = String;
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                parse_mono(&X.to_string(), s).map(Self)
            }
        }

        #[cfg(feature = "serde")]
        impl<const X: char> serde::Serialize for Univar<X, $I> {
            fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
            where S: serde::Serializer {
                serializer.serialize_str(&self.fmt_impl(false))
            }
        }
        
        #[cfg(feature = "serde")]
        impl<'de, const X: char> serde::Deserialize<'de> for Univar<X, $I> {
            fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
            where D: serde::Deserializer<'de> {
                let s = String::deserialize(deserializer)?;
                Self::from_str(&s).map_err(serde::de::Error::custom)
            }
        }        

        impl<const X: char> Elem for Univar<X, $I> { 
            fn math_symbol() -> String {
                format!("{X}")
            }
        }
        
        impl<const X: char> Gen for Univar<X, $I> {}
    }
}

macro_rules! impl_univar_unsigned {
    ($I:ty) => {
        impl_univar!($I);

        impl<const X: char> Mono for Univar<X, $I> {
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
        impl_univar!($I);

        impl<const X: char> Mono for Univar<X, $I> {
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

pub(crate) fn fmt_mono<I>(x: &str, d: I, unicode: bool) -> String
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
        format!("{x}^{d}")
    }
}

pub(crate) fn parse_mono<I, E>(x: &str, s: &str) -> Result<I, String>
where I: FromStr<Err = E> + FromPrimitive, E: ToString {
    if s == "1" { 
        Ok(I::from_i32(0).unwrap())
    } else if s == x { 
        Ok(I::from_i32(1).unwrap())
    } else if let Some(d) = s.strip_prefix(&format!("{x}^")) { 
        match I::from_str(d) { 
            Ok(d) => Ok(d),
            Err(e) => Err(e.to_string())
        }
    } else { 
        Err(format!("Failed to parse: {s}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn init() { 
        type M = Univar<'X',usize>;
        let x = M::from;

        let d = x(2);

        assert_eq!(d.0, 2);
        assert_eq!(d.deg(), 2);
        assert_eq!(M::var_symbol(), 'X');
    }

    #[test]
    fn from_str() { 
        type M = Univar<'X',isize>;
        let x = M::from;
        
        assert_eq!(M::from_str("1"), Ok(M::one()));
        assert_eq!(M::from_str("X"), Ok(x(1)));
        assert_eq!(M::from_str("X^2"), Ok(x(2)));
        assert_eq!(M::from_str("X^-2"), Ok(x(-2)));
        assert!(M::from_str("2").is_err());
        assert!(M::from_str("x").is_err());
    }

    #[test]
    fn display() { 
        type M = Univar<'X', isize>;
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
        type M = Univar<'X',usize>;
        let x = M::from;

        let d = x(0);
        assert_eq!(d.inv(), Some(x(0)));

        let d = x(1);
        assert_eq!(d.inv(), None);
    }

    #[test]
    fn neg_opt_signed() { 
        type M = Univar<'X',isize>;
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
        type M = Univar<'X', usize>;
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
        type M = Univar<'X', isize>;
        let x = M::from;

        assert!(x(0) < x(1));
        assert!(x(1) < x(2));
        assert!(x(-1) < x(0));
    }

    #[test]
    #[cfg(feature = "serde")]
    fn serialize() { 
        type M = Univar<'X', isize>;
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

        let d = x(-2);
        let ser = serde_json::to_string(&d).unwrap();
        let des = serde_json::from_str::<M>(&ser).unwrap();
        
        assert_eq!(&ser, "\"X^-2\"");
        assert_eq!(d, des);
    }
}