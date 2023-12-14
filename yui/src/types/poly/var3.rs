use core::panic;
use std::fmt::{Display, Debug};
use std::hash::Hash;
use std::ops::{AddAssign, Mul, MulAssign, DivAssign, SubAssign, Div, Add};
use std::str::FromStr;
use itertools::Itertools;
use num_traits::{Zero, One, Pow, FromPrimitive, ToPrimitive};
use auto_impl_ops::auto_ops;

use crate::{Elem, ElemBase};
use crate::lc::Gen;

use super::Mono;
use super::var::{fmt_mono, parse_mono};

// `Var3<X, Y, Z, I>` : represents trivariant monomials X^i Y^j Z^k.
// `I` is either `usize` or `isize`.

#[derive(Clone, PartialEq, Eq, Hash, Default)]
pub struct Var3<const X: char, const Y: char, const Z: char, I>(
    I, I, I
);

impl<const X: char, const Y: char, const Z: char, I> Var3<X, Y, Z, I> {
    pub fn var_symbol(i: usize) -> char { 
        assert!(i < 3);
        match i { 
            0 => X,
            1 => Y,
            2 => Z,
            _ => panic!()
        }
    }

    pub fn deg_for(&self, i: usize) -> I
    where I: Copy { 
        assert!(i < 3);
        match i { 
            0 => self.0,
            1 => self.1,
            2 => self.2,
            _ => panic!()
        }
    }

    pub fn total_deg(&self) -> I
    where I: Copy + for<'x> Add<&'x I, Output = I> { 
        self.0 + &self.1 + &self.2
    }

    pub fn eval<R>(&self, x: &R, y: &R, z: &R) -> R
    where R: Mul<Output = R>, for<'x, 'y> &'x R: Pow<&'y I, Output = R> {
        x.pow(&self.0) * y.pow(&self.1) * z.pow(&self.2)
    }

    fn fmt_impl(&self, unicode: bool) -> String
    where I: ToPrimitive { 
        let Var3(d0, d1, d2) = self;
        let s = [(X, d0), (Y, d1), (Z, d2)].into_iter().map(|(x, d)|
            fmt_mono(&x.to_string(), d, unicode)
        ).filter(|s| s != "1").join("");
        if s == "" { 
            "1".to_string()
        } else { 
            s
        }
    }
}

impl<const X: char, const Y: char, const Z: char, I> From<(I, I, I)> for Var3<X, Y, Z, I> {
    fn from(d: (I, I, I)) -> Self {
        Self(d.0, d.1, d.2)
    }
}

impl<const X: char, const Y: char, const Z: char, I> FromStr for Var3<X, Y, Z, I>
where I: Zero + FromStr + FromPrimitive + Debug, <I as FromStr>::Err: ToString {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use regex::Regex;
        if s == "1" { 
            return Ok(Self(I::zero(), I::zero(), I::zero()))
        }

        let p1 = format!(r"{X}(\^-?[0-9]+)?");
        let p2 = format!(r"{Y}(\^-?[0-9]+)?");
        let p3 = format!(r"{Z}(\^-?[0-9]+)?");
        let p_all = format!(r"^({p1})?\s?({p2})?\s?({p3})?$");

        if !Regex::new(&p_all).unwrap().is_match(&s) { 
            return Err(format!("Failed to parse: {s}"))
        }

        let r1 = Regex::new(&p1).unwrap();
        let r2 = Regex::new(&p2).unwrap();
        let r3 = Regex::new(&p3).unwrap();

        let d = [(X, r1), (Y, r2), (Z, r3)].into_iter().map(|(x, r)| { 
            if let Some(c) = r.captures(s) { 
                parse_mono(&x.to_string(), &c[0])
            } else { 
                Ok(I::zero())
            }
        }).collect::<Result<Vec<_>, _>>()?;
        let [d1, d2, d3] = d.try_into().unwrap();

        Ok(Self(d1, d2, d3))
    }
}

impl<const X: char, const Y: char, const Z: char, I> One for Var3<X, Y, Z, I>
where I: for<'x >AddAssign<&'x I> + Zero {
    fn one() -> Self {
        Self::from((I::zero(), I::zero(), I::zero())) // x^0 = 1.
    }
}

#[auto_ops]
impl<const X: char, const Y: char, const Z: char, I> MulAssign<&Var3<X, Y, Z, I>> for Var3<X, Y, Z, I>
where I: for<'x >AddAssign<&'x I> {
    fn mul_assign(&mut self, rhs: &Var3<X, Y, Z, I>) {
        self.0 += &rhs.0; // x^i * x^j = x^{i+j}
        self.1 += &rhs.1;
        self.2 += &rhs.2;
    }
}

#[auto_ops]
impl<const X: char, const Y: char, const Z: char, I> DivAssign<&Var3<X, Y, Z, I>> for Var3<X, Y, Z, I>
where I: for<'x >SubAssign<&'x I> {
    fn div_assign(&mut self, rhs: &Var3<X, Y, Z, I>) {
        self.0 -= &rhs.0; // x^i / x^j = x^{i-j}
        self.1 -= &rhs.1;
        self.2 -= &rhs.2;
    }
}

impl<const X: char, const Y: char, const Z: char, I> PartialOrd for Var3<X, Y, Z, I>
where I: Copy + Eq + Ord + for<'x> Add<&'x I, Output = I> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(Ord::cmp(self, other))
    }
}

impl<const X: char, const Y: char, const Z: char, I> Ord for Var3<X, Y, Z, I>
where I: Copy + Eq + Ord + for<'x> Add<&'x I, Output = I> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::*;

        Ord::cmp(&self.total_deg(), &other.total_deg())
        .then_with(|| 
            Ord::cmp(&self.0, &other.0)
        ).then_with(|| 
            Ord::cmp(&self.1, &other.1)
        ).then_with(|| 
            Ord::cmp(&self.2, &other.2)
        )
    }
}

impl<const X: char, const Y: char, const Z: char, I> Display for Var3<X, Y, Z, I>
where I: ToPrimitive { 
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = self.fmt_impl(true);
        f.write_str(&s)
    }
}

impl<const X: char, const Y: char, const Z: char, I> Debug for Var3<X, Y, Z, I>
where I: ToPrimitive { 
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Display::fmt(self, f)
    }
}

#[cfg(feature = "serde")]
impl<const X: char, const Y: char, const Z: char, I> serde::Serialize for Var3<X, Y, Z, I>
where I: ToPrimitive {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where S: serde::Serializer {
        serializer.serialize_str(&self.fmt_impl(false))
    }
}

#[cfg(feature = "serde")]
impl<'de, const X: char, const Y: char, const Z: char, I> serde::Deserialize<'de> for Var3<X, Y, Z, I>
where Self: FromStr, <Self as FromStr>::Err: Display {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where D: serde::Deserializer<'de> {
        let s = String::deserialize(deserializer)?;
        Self::from_str(&s).map_err(serde::de::Error::custom)
    }
}        

impl<const X: char, const Y: char, const Z: char, I> Elem for Var3<X, Y, Z, I>
where I: ElemBase + ToPrimitive { 
    fn math_symbol() -> String {
        format!("{X}, {Y}, {Z}")
    }
}
        
impl<const X: char, const Y: char, const Z: char, I> Gen for Var3<X, Y, Z, I>
where I: ElemBase + Copy + Hash + Ord + for<'x> Add<&'x I, Output = I> + ToPrimitive {}

macro_rules! impl_trivar_unsigned {
    ($I:ty) => {
        impl<const X: char, const Y: char, const Z: char> Mono for Var3<X, Y, Z, $I> {
            type Deg = ($I, $I, $I);

            fn deg(&self) -> Self::Deg {
                (self.0, self.1, self.2)
            }

            fn is_unit(&self) -> bool { 
                self.0.is_zero() && self.1.is_zero() && self.2.is_zero()
            }

            fn inv(&self) -> Option<Self> { // (x^i)^{-1} = x^{-i}
                if self.is_unit() {
                    Some(Self(0, 0, 0))
                } else { 
                    None
                }
            }

            fn divides(&self, other: &Self) -> bool { 
                self.0 <= other.0 && self.1 <= other.1 && self.2 <= other.2
            }
        }
    };
}

macro_rules! impl_trivar_signed {
    ($I:ty) => {
        impl<const X: char, const Y: char, const Z: char> Mono for Var3<X, Y, Z, $I> {
            type Deg = ($I, $I, $I);

            fn deg(&self) -> Self::Deg {
                (self.0, self.1, self.2)
            }

            fn is_unit(&self) -> bool { 
                true
            }

            fn inv(&self) -> Option<Self> { // (x^i)^{-1} = x^{-i}
                Some(Self(-self.0, -self.1, -self.2))
            }

            fn divides(&self, _other: &Self) -> bool { 
                true
            }
        }
    };
}

impl_trivar_unsigned!(usize);
impl_trivar_signed!  (isize);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn var_symbol() { 
        type M = Var3<'X','Y','Z',usize>;
        
        assert_eq!(M::var_symbol(0), 'X');
        assert_eq!(M::var_symbol(1), 'Y');
        assert_eq!(M::var_symbol(2), 'Z');
    }

    #[test]
    fn init() { 
        type M = Var3<'X','Y','Z',usize>;
        let xyz = |i, j, k| M::from((i, j, k));

        let d = xyz(2, 3, 1);

        assert_eq!(d.0, 2);
        assert_eq!(d.1, 3);
        assert_eq!(d.2, 1);
    }

    #[test]
    fn from_str() { 
        type M = Var3<'X','Y','Z',usize>;
        let xyz = |i, j, k| M::from((i, j, k));
        
        assert_eq!(M::from_str("1"), Ok(M::one()));
        assert_eq!(M::from_str("X"), Ok(xyz(1, 0, 0)));
        assert_eq!(M::from_str("Y"), Ok(xyz(0, 1, 0)));
        assert_eq!(M::from_str("Z"), Ok(xyz(0, 0, 1)));
        assert_eq!(M::from_str("X^2"), Ok(xyz(2, 0, 0)));
        assert_eq!(M::from_str("Y^2"), Ok(xyz(0, 2, 0)));
        assert_eq!(M::from_str("Z^2"), Ok(xyz(0, 0, 2)));
        assert_eq!(M::from_str("XYZ"), Ok(xyz(1, 1, 1)));
        assert_eq!(M::from_str("X^2Y^3Z"), Ok(xyz(2, 3, 1)));
        assert!(M::from_str("2").is_err());
    }

    #[test]
    fn display() { 
        type M = Var3<'X','Y','Z',usize>;
        let xyz = |i, j, k| M::from((i, j, k));

        let d = xyz(0, 0, 0);
        assert_eq!(&d.to_string(), "1");

        let d = xyz(1, 0, 0);
        assert_eq!(&d.to_string(), "X");

        let d = xyz(2, 0, 0);
        assert_eq!(&d.to_string(), "X²");

        let d = xyz(0, 1, 0);
        assert_eq!(&d.to_string(), "Y");

        let d = xyz(0, 2, 0);
        assert_eq!(&d.to_string(), "Y²");

        let d = xyz(1, 1, 0);
        assert_eq!(&d.to_string(), "XY");

        let d = xyz(2, 3, 1);
        assert_eq!(&d.to_string(), "X²Y³Z");
    }

    #[test]
    fn neg_opt_unsigned() { 
        type M = Var3<'X','Y','Z',usize>;
        let xyz = |i, j, k| M::from((i, j, k));

        let d = xyz(0, 0, 0);
        assert_eq!(d.inv(), Some(xyz(0, 0, 0)));

        let d = xyz(1, 0, 0);
        assert_eq!(d.inv(), None);

        let d = xyz(0, 1, 0);
        assert_eq!(d.inv(), None);

        let d = xyz(0, 0, 1);
        assert_eq!(d.inv(), None);
    }

    #[test]
    fn neg_opt_signed() { 
        type M = Var3<'X','Y','Z',isize>;
        let xyz = |i, j,k| M::from((i, j, k));

        let d = xyz(0, 0, 0);
        assert_eq!(d.inv(), Some(xyz(0, 0, 0)));

        let d = xyz(1, 0, 0);
        assert_eq!(d.inv(), Some(xyz(-1, 0, 0)));

        let d = xyz(0, 1, 0);
        assert_eq!(d.inv(), Some(xyz(0, -1, 0)));

        let d = xyz(2, 3, 1);
        assert_eq!(d.inv(), Some(xyz(-2, -3, -1)));
    }

    #[test]
    fn eval() { 
        type M = Var3<'X','Y','Z',usize>;
        let xyz = |i, j,k| M::from((i, j, k));

        let d = xyz(0, 0, 0);
        assert_eq!(d.eval::<i32>(&2, &3, &4), 1);

        let d = xyz(1, 0, 0);
        assert_eq!(d.eval::<i32>(&2, &3, &4), 2);

        let d = xyz(0, 1, 0);
        assert_eq!(d.eval::<i32>(&2, &3, &4), 3);

        let d = xyz(1, 1, 1);
        assert_eq!(d.eval::<i32>(&2, &3, &4), 24);

        let d = xyz(2, 3, 1);
        assert_eq!(d.eval::<i32>(&2, &3, &4), 432);
    }

    #[test]
    fn ord() { 
        type M = Var3<'X','Y','Z',usize>;
        let xyz = |i, j,k| M::from((i, j, k));

        assert!(xyz(2, 1, 1) < xyz(1, 3, 1));
        assert!(xyz(1, 2, 1) < xyz(2, 1, 1));
        assert!(xyz(1, 2, 1) < xyz(1, 3, 1));
        assert!(xyz(1, 2, 2) < xyz(2, 1, 2));
        assert!(xyz(1, 2, 1) < xyz(1, 2, 2));
    }
}