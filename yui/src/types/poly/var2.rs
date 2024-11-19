use core::panic;
use std::fmt::{Display, Debug};
use std::hash::Hash;
use std::ops::{AddAssign, Mul, MulAssign, DivAssign, SubAssign, Div, Add};
use std::str::FromStr;
use num_traits::{Zero, One, Pow, FromPrimitive, ToPrimitive};
use auto_impl_ops::auto_ops;

use crate::{Elem, ElemBase};
use crate::lc::{Gen, OrdForDisplay};

use super::{Mono, MonoOrd};
use super::var::parse_mono_deg;
use super::mvar::fmt_mono_n;

// `Var2<X, Y, I>` : represents bivariant monomials X^i Y^j.
// `I` is either `usize` or `isize`.

#[derive(Clone, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde_with::DeserializeFromStr))]
pub struct Var2<const X: char, const Y: char, I> (
    I, I
);

impl<const X: char, const Y: char, I> Var2<X, Y, I> {
    pub fn var_symbol(i: usize) -> char { 
        assert!(i < 2);
        match i { 
            0 => X,
            1 => Y,
            _ => panic!()
        }
    }

    pub fn deg_for(&self, i: usize) -> I
    where I: Copy { 
        assert!(i < 2);
        match i { 
            0 => self.0,
            1 => self.1,
            _ => panic!()
        }
    }

    pub fn total_deg(&self) -> I
    where I: Copy + for<'x> Add<&'x I, Output = I> { 
        self.0 + &self.1
    }

    pub fn eval<R>(&self, x: &R, y: &R) -> R
    where R: Mul<Output = R>, for<'x, 'y> &'x R: Pow<&'y I, Output = R> {
        x.pow(&self.0) * y.pow(&self.1)
    }

    fn to_string_u(&self, unicode: bool) -> String
    where I: ToPrimitive { 
        let Var2(d0, d1) = self;
        let seq = [(X, d0), (Y, d1)];
        fmt_mono_n(seq, unicode)
    }
}

impl<const X: char, const Y: char, I> From<(I, I)> for Var2<X, Y, I> {
    fn from(d: (I, I)) -> Self {
        Self(d.0, d.1)
    }
}

impl<const X: char, const Y: char, I> FromStr for Var2<X, Y, I>
where I: Zero + AddAssign + FromStr + FromPrimitive {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use regex::Regex;

        if s == "1" { 
            return Ok(Self(I::zero(), I::zero()))
        }

        let p = format!(r"({X}|{Y})(\^\{{?-?[0-9]+\}}?)?");
        let p_all = format!(r"^({p}\s?)+$");

        let r = Regex::new(&p).unwrap();
        let r_all = Regex::new(&p_all).unwrap();

        if !r_all.is_match(&s) { 
            return Err(format!("Failed to parse: {s}"))
        }

        let mut deg = (I::zero(), I::zero());
        
        for c in r.captures_iter(&s) {
            let x = &c[1];
            let i = parse_mono_deg(x, &c[0]).unwrap();
            if x.starts_with(X) { 
                deg.0 += i;
            } else { 
                deg.1 += i;
            }
        };

        Ok(Self::from(deg))
    }
}

impl<const X: char, const Y: char, I> One for Var2<X, Y, I>
where I: for<'x >AddAssign<&'x I> + Zero {
    fn one() -> Self {
        Self::from((I::zero(), I::zero())) // x^0 = 1.
    }
}

#[auto_ops]
impl<const X: char, const Y: char, I> MulAssign<&Var2<X, Y, I>> for Var2<X, Y, I>
where I: for<'x >AddAssign<&'x I> {
    fn mul_assign(&mut self, rhs: &Var2<X, Y, I>) {
        self.0 += &rhs.0; // x^i * x^j = x^{i+j}
        self.1 += &rhs.1; // x^i * x^j = x^{i+j}
    }
}

#[auto_ops]
impl<const X: char, const Y: char, I> DivAssign<&Var2<X, Y, I>> for Var2<X, Y, I>
where I: for<'x >SubAssign<&'x I> {
    fn div_assign(&mut self, rhs: &Var2<X, Y, I>) {
        self.0 -= &rhs.0; // x^i * x^j = x^{i+j}
        self.1 -= &rhs.1;
    }
}

impl<const X: char, const Y: char, I> MonoOrd for Var2<X, Y, I>
where I: Copy + Eq + Ord + for<'x> Add<&'x I, Output = I> {
    fn cmp_lex(&self, other: &Self) -> std::cmp::Ordering {
        // must have x_0 > x_1
        I::cmp(&self.0, &other.0).then_with(|| 
            I::cmp(&self.1, &other.1)
        )
    }

    fn cmp_grlex(&self, other: &Self) -> std::cmp::Ordering {
        I::cmp(&self.total_deg(), &other.total_deg()).then_with(|| 
            Self::cmp_lex(self, other)
        )
    }
}

impl<const X: char, const Y: char, I> OrdForDisplay for Var2<X, Y, I>
where I: Copy + Eq + Ord + for<'x> Add<&'x I, Output = I> {
    fn cmp_for_display(&self, other: &Self) -> std::cmp::Ordering {
        Self::cmp_lex(self, other)
    }
}

impl<const X: char, const Y: char, I> Display for Var2<X, Y, I>
where I: ToPrimitive { 
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = self.to_string_u(true);
        f.write_str(&s)
    }
}

impl<const X: char, const Y: char, I> Debug for Var2<X, Y, I>
where I: ToPrimitive { 
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Display::fmt(self, f)
    }
}

#[cfg(feature = "serde")]
impl<const X: char, const Y: char, I> serde::Serialize for Var2<X, Y, I>
where I: ToPrimitive {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where S: serde::Serializer {
        serializer.serialize_str(&self.to_string_u(false))
    }
}

impl<const X: char, const Y: char, I> Elem for Var2<X, Y, I>
where I: ElemBase + ToPrimitive { 
    fn math_symbol() -> String {
        format!("{X}, {Y}")
    }
}
        
impl<const X: char, const Y: char, I> Gen for Var2<X, Y, I>
where I: ElemBase + Copy + Hash + Ord + for<'x> Add<&'x I, Output = I> + ToPrimitive {}

macro_rules! impl_bivar_unsigned {
    ($I:ty) => {
        impl<const X: char, const Y: char> Mono for Var2<X, Y, $I> {
            type Deg = ($I, $I);

            fn deg(&self) -> Self::Deg {
                (self.0, self.1)
            }

            fn is_unit(&self) -> bool { 
                self.0.is_zero() && self.1.is_zero()
            }

            fn inv(&self) -> Option<Self> { // (x^i)^{-1} = x^{-i}
                if self.is_unit() {
                    Some(Self(0, 0))
                } else { 
                    None
                }
            }

            fn divides(&self, other: &Self) -> bool { 
                self.0 <= other.0 && self.1 <= other.1
            }
        }
    };
}

macro_rules! impl_bivar_signed {
    ($I:ty) => {
        impl<const X: char, const Y: char> Mono for Var2<X, Y, $I> {
            type Deg = ($I, $I);

            fn deg(&self) -> Self::Deg {
                (self.0, self.1)
            }

            fn is_unit(&self) -> bool { 
                true
            }

            fn inv(&self) -> Option<Self> { // (x^i)^{-1} = x^{-i}
                Some(Self(-self.0, -self.1))
            }

            fn divides(&self, _other: &Self) -> bool { 
                true
            }
        }
    };
}

impl_bivar_unsigned!(usize);
impl_bivar_signed!  (isize);

#[cfg(feature = "tex")] 
mod tex {
    use crate::TeX;
    use super::*;

    impl<const X: char, const Y: char, I> TeX for Var2<X, Y, I>
    where I: ToPrimitive {
        fn tex_math_symbol() -> String { 
            format!("{},{}", X, Y)
        }
        fn tex_string(&self) -> String {
            self.to_string_u(false)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn var_symbol() { 
        type M = Var2<'X','Y',usize>;
        
        assert_eq!(M::var_symbol(0), 'X');
        assert_eq!(M::var_symbol(1), 'Y');
    }

    #[test]
    fn init() { 
        type M = Var2<'X','Y',usize>;
        let xy = |i, j| M::from((i, j));

        let d = xy(2, 3);

        assert_eq!(d.0, 2);
        assert_eq!(d.1, 3);
    }

    #[test]
    fn from_str() { 
        type M = Var2<'X','Y', isize>;
        let xy = |i, j| M::from((i, j));
        
        assert_eq!(M::from_str("1"), Ok(M::one()));
        assert_eq!(M::from_str("X"), Ok(xy(1, 0)));
        assert_eq!(M::from_str("Y"), Ok(xy(0, 1)));
        assert_eq!(M::from_str("X^2"), Ok(xy(2, 0)));
        assert_eq!(M::from_str("Y^2"), Ok(xy(0, 2)));
        assert_eq!(M::from_str("XY"), Ok(xy(1, 1)));
        assert_eq!(M::from_str("X^2Y^3"), Ok(xy(2, 3)));
        assert_eq!(M::from_str(r"X^{21}Y^{32}"), Ok(xy(21, 32)));
        assert_eq!(M::from_str(r"X^{-2}Y^{-3}"), Ok(xy(-2, -3)));
        assert!(M::from_str("2").is_err());
    }

    #[test]
    fn display() { 
        type M = Var2<'X','Y',usize>;
        let xy = |i, j| M::from((i, j));

        let d = xy(0, 0);
        assert_eq!(&d.to_string(), "1");

        let d = xy(1, 0);
        assert_eq!(&d.to_string(), "X");

        let d = xy(2, 0);
        assert_eq!(&d.to_string(), "X²");

        let d = xy(0, 1);
        assert_eq!(&d.to_string(), "Y");

        let d = xy(0, 2);
        assert_eq!(&d.to_string(), "Y²");

        let d = xy(1, 1);
        assert_eq!(&d.to_string(), "XY");

        let d = xy(2, 3);
        assert_eq!(&d.to_string(), "X²Y³");
    }

    #[test]
    fn neg_opt_unsigned() { 
        type M = Var2<'X','Y',usize>;
        let xy = |i, j| M::from((i, j));

        let d = xy(0, 0);
        assert_eq!(d.inv(), Some(xy(0, 0)));

        let d = xy(1, 0);
        assert_eq!(d.inv(), None);

        let d = xy(0, 1);
        assert_eq!(d.inv(), None);

        let d = xy(1, 1);
        assert_eq!(d.inv(), None);
    }

    #[test]
    fn neg_opt_signed() { 
        type M = Var2<'X','Y',isize>;
        let xy = |i, j| M::from((i, j));

        let d = xy(0, 0);
        assert_eq!(d.inv(), Some(xy(0, 0)));

        let d = xy(1, 0);
        assert_eq!(d.inv(), Some(xy(-1, 0)));

        let d = xy(0, 1);
        assert_eq!(d.inv(), Some(xy(0, -1)));

        let d = xy(2, 3);
        assert_eq!(d.inv(), Some(xy(-2, -3)));
    }

    #[test]
    fn eval() { 
        type M = Var2<'X','Y',usize>;
        let xy = |i, j| M::from((i, j));

        let d = xy(0, 0);
        assert_eq!(d.eval::<i32>(&2, &3), 1);

        let d = xy(1, 0);
        assert_eq!(d.eval::<i32>(&2, &3), 2);

        let d = xy(0, 1);
        assert_eq!(d.eval::<i32>(&2, &3), 3);

        let d = xy(1, 1);
        assert_eq!(d.eval::<i32>(&2, &3), 6);

        let d = xy(2, 3);
        assert_eq!(d.eval::<i32>(&2, &3), 108);
    }

    #[test]
    fn cmp_lex() { 
        type M = Var2<'X','Y',usize>;
        let xy = |i, j| M::from((i, j));

        // x^2 y > x y^2 > x > y^2 > y > 1
        assert!(Var2::cmp_lex(&xy(2, 1), &xy(1, 2)).is_gt());
        assert!(Var2::cmp_lex(&xy(1, 2), &xy(1, 0)).is_gt());
        assert!(Var2::cmp_lex(&xy(1, 0), &xy(0, 2)).is_gt());
        assert!(Var2::cmp_lex(&xy(0, 2), &xy(0, 1)).is_gt());
        assert!(Var2::cmp_lex(&xy(0, 1), &xy(0, 0)).is_gt());
    }

    #[test]
    fn cmp_grlex() { 
        type M = Var2<'X','Y',usize>;
        let xy = |i, j| M::from((i, j));

        // x^2 y > x y^2 > y^2 > x > y > 1
        assert!(Var2::cmp_grlex(&xy(2, 1), &xy(1, 2)).is_gt());
        assert!(Var2::cmp_grlex(&xy(1, 2), &xy(0, 2)).is_gt());
        assert!(Var2::cmp_grlex(&xy(0, 2), &xy(1, 0)).is_gt());
        assert!(Var2::cmp_grlex(&xy(1, 0), &xy(0, 1)).is_gt());
        assert!(Var2::cmp_grlex(&xy(0, 1), &xy(0, 0)).is_gt());
    }
}