use core::panic;
use std::fmt::{Display, Debug};
use std::ops::{AddAssign, Mul, MulAssign, DivAssign, SubAssign, Div, Add};
use std::str::FromStr;
use num_traits::{Zero, One, Pow};
use auto_impl_ops::auto_ops;

use crate::Elem;
use crate::lc::Gen;

use super::Mono;
use super::univar::fmt_mono;

// `BiVar<X, Y, I>` : represents bivariant monomials X^i Y^j.
// `I` is either `usize` or `isize`.

#[derive(Clone, PartialEq, Eq, Hash, Default)]
pub struct BiVar<const X: char, const Y: char, I>(
    I, 
    I
);

impl<const X: char, const Y: char, I> BiVar<X, Y, I> {
    pub fn var_symbol(i: usize) -> char { 
        assert!(i < 2);
        match i { 
            0 => X,
            1 => Y,
            _ => panic!()
        }
    }

    pub fn eval<R>(&self, x: &R, y: &R) -> R
    where R: Mul<Output = R>, for<'x, 'y> &'x R: Pow<&'y I, Output = R> {
        x.pow(&self.0) * y.pow(&self.1)
    }
}

impl<const X: char, const Y: char, I> BiVar<X, Y, I>
where I: Clone {
    pub fn deg_for(&self, i: usize) -> I { 
        assert!(i < 2);
        match i { 
            0 => self.0.clone(),
            1 => self.1.clone(),
            _ => panic!()
        }
    }
}

impl<const X: char, const Y: char, I> BiVar<X, Y, I>
where for<'x> &'x I: Add<&'x I, Output = I> {
    pub fn total_deg(&self) -> I { 
        &self.0 + &self.1
    }
}

impl<const X: char, const Y: char, I> From<(I, I)> for BiVar<X, Y, I> {
    fn from(d: (I, I)) -> Self {
        Self(d.0, d.1)
    }
}

impl<const X: char, const Y: char, I> FromStr for BiVar<X, Y, I>
where I: Zero + One {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() == 1 { 
            let x = s.chars().next().unwrap();
            if x == '1' { 
                return Ok(Self(I::zero(), I::zero()));
            } else if x == X { 
                return Ok(Self(I::one(),  I::zero()));
            } else if x == Y { 
                return Ok(Self(I::zero(),  I::one()));
            }
        }
        
        // TODO support more complex format. 
        Err(())
    }
}

impl<const X: char, const Y: char, I> One for BiVar<X, Y, I>
where I: for<'x >AddAssign<&'x I> + Zero {
    fn one() -> Self {
        Self::from((I::zero(), I::zero())) // x^0 = 1.
    }
}

#[auto_ops]
impl<const X: char, const Y: char, I> MulAssign<&BiVar<X, Y, I>> for BiVar<X, Y, I>
where I: for<'x >AddAssign<&'x I> {
    fn mul_assign(&mut self, rhs: &BiVar<X, Y, I>) {
        self.0 += &rhs.0; // x^i * x^j = x^{i+j}
        self.1 += &rhs.1; // x^i * x^j = x^{i+j}
    }
}

#[auto_ops]
impl<const X: char, const Y: char, I> DivAssign<&BiVar<X, Y, I>> for BiVar<X, Y, I>
where I: for<'x >SubAssign<&'x I> {
    fn div_assign(&mut self, rhs: &BiVar<X, Y, I>) {
        self.0 -= &rhs.0; // x^i * x^j = x^{i+j}
        self.1 -= &rhs.1;
    }
}

impl<const X: char, const Y: char, I> PartialOrd for BiVar<X, Y, I>
where I: Eq + Ord, for<'x> &'x I: Add<&'x I, Output = I> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(Ord::cmp(self, other))
    }
}

impl<const X: char, const Y: char, I> Ord for BiVar<X, Y, I>
where I: Eq + Ord, for<'x> &'x I: Add<&'x I, Output = I> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::*;

        Ord::cmp(&self.total_deg(), &other.total_deg())
        .then_with(|| 
            Ord::cmp(&self.0, &other.0)
        ).then_with(|| 
            Ord::cmp(&self.1, &other.1)
        )
    }
}

macro_rules! impl_bivar {
    ($I:ty) => {
        impl<const X: char, const Y: char> Display for BiVar<X, Y, $I> { 
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let BiVar(d0, d1) = self;
                let x = fmt_mono(X.to_string(), *d0);
                let y = fmt_mono(Y.to_string(), *d1);
                
                match (x.as_str(), y.as_str()) {
                    ("1", "1") => write!(f, "1"),
                    ( _ , "1") => write!(f, "{x}"),
                    ("1",  _ ) => write!(f, "{y}"),
                    _          => write!(f, "{x}{y}")
                }
            }
        }
        
        impl<const X: char, const Y: char> Debug for BiVar<X, Y, $I> { 
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                Display::fmt(self, f)
            }
        }

        impl<const X: char, const Y: char> Elem for BiVar<X, Y, $I> { 
            fn math_symbol() -> String {
                format!("{X}, {Y}")
            }
        }
                
        impl<const X: char, const Y: char> Gen for BiVar<X, Y, $I> {}
    }
}

macro_rules! impl_bivar_unsigned {
    ($I:ty) => {
        impl_bivar!($I);

        impl<const X: char, const Y: char> Mono for BiVar<X, Y, $I> {
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
        impl_bivar!($I);

        impl<const X: char, const Y: char> Mono for BiVar<X, Y, $I> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn var_symbol() { 
        type M = BiVar<'X','Y',usize>;
        
        assert_eq!(M::var_symbol(0), 'X');
        assert_eq!(M::var_symbol(1), 'Y');
    }

    #[test]
    fn init() { 
        type M = BiVar<'X','Y',usize>;
        let xy = |i, j| M::from((i, j));

        let d = xy(2, 3);

        assert_eq!(d.0, 2);
        assert_eq!(d.1, 3);
    }

    #[test]
    fn from_str() { 
        type M = BiVar<'X','Y',usize>;
        let xy = |i, j| BiVar::<'X','Y',usize>(i, j);
        
        assert_eq!(M::from_str("1"), Ok(M::one()));
        assert_eq!(M::from_str("X"), Ok(xy(1, 0)));
        assert_eq!(M::from_str("Y"), Ok(xy(0, 1)));
        assert_eq!(M::from_str("2"), Err(()));
    }

    #[test]
    fn display() { 
        type M = BiVar<'X','Y',usize>;
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
        type M = BiVar<'X','Y',usize>;
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
        type M = BiVar<'X','Y',isize>;
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
        type M = BiVar<'X','Y',usize>;
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
    fn ord() { 
        type M = BiVar<'X','Y',usize>;
        let xy = |i, j| M::from((i, j));

        assert!(xy(2, 1) < xy(1, 3));
        assert!(xy(1, 2) < xy(2, 1));
        assert!(xy(1, 2) < xy(1, 3));
    }
}