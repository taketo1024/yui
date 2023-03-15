use std::fmt::Display;
use std::ops::{AddAssign, Mul, MulAssign};
use std::str::FromStr;
use num_traits::{Zero, One, Pow};
use auto_impl_ops::auto_ops;

use yui_core::Elem;
use yui_lin_comb::FreeGen;
use yui_utils::map;

use crate::{PolyDeg, PolyGen, MDegree};
use super::mono::fmt_mono;

// `Mono2<X, Y, I>` : a struct representing X^i Y^j.
// `I` is one of `usize`, `isize`.

#[derive(Clone, PartialEq, Eq, Hash, Default, Debug, PartialOrd, Ord)]
pub struct Mono2<const X: char, const Y: char, I>(I, I);

impl<const X: char, const Y: char, I> Mono2<X, Y, I> {
    pub fn new(d0: I, d1: I) -> Self { 
        Self(d0, d1)
    }
}

impl<const X: char, const Y: char, I> From<(I, I)> for Mono2<X, Y, I> {
    fn from(d: (I, I)) -> Self {
        Self(d.0, d.1)
    }
}

impl<const X: char, const Y: char, I> From<MDegree<I>> for Mono2<X, Y, I>
where I: Zero + Clone {
    fn from(mdeg: MDegree<I>) -> Self {
        assert!(mdeg.len() <= 2);
        let d0 = mdeg.deg(0);
        let d1 = mdeg.deg(1);
        Self(d0, d1)
    }
}

impl<const X: char, const Y: char, I> FromStr for Mono2<X, Y, I>
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

impl<const X: char, const Y: char, I> Mono2<X, Y, I>
where I: Zero + Clone {
    pub fn deg(&self, i: usize) -> I { 
        match i { 
            0 => self.0.clone(),
            1 => self.1.clone(),
            _ => panic!()
        }
    }
}

impl<const X: char, const Y: char, I> One for Mono2<X, Y, I>
where I: for<'x >AddAssign<&'x I> + Zero {
    fn one() -> Self {
        Self::from((I::zero(), I::zero())) // x^0 = 1.
    }
}

#[auto_ops]
impl<const X: char, const Y: char, I> MulAssign<&Mono2<X, Y, I>> for Mono2<X, Y, I>
where I: for<'x >AddAssign<&'x I> {
    fn mul_assign(&mut self, rhs: &Mono2<X, Y, I>) {
        self.0 += &rhs.0; // x^i * x^j = x^{i+j}
        self.1 += &rhs.1; // x^i * x^j = x^{i+j}
    }
}

macro_rules! impl_poly_gen {
    ($I:ty) => {
        impl<const X: char, const Y: char> Display for Mono2<X, Y, $I> { 
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let Mono2(d0, d1) = self;
                let x = fmt_mono(X.to_string(), *d0 as isize);
                let y = fmt_mono(Y.to_string(), *d1 as isize);
                
                match (x.as_str(), y.as_str()) {
                    ("1", "1") => write!(f, "1"),
                    ( _ , "1") => write!(f, "{x}"),
                    ("1",  _ ) => write!(f, "{y}"),
                    _          => write!(f, "{x}{y}")
                }
            }
        }
        
        impl<const X: char, const Y: char> Elem for Mono2<X, Y, $I> { 
            fn set_symbol() -> String {
                format!("{X}, {Y}")
            }
        }
                
        impl<const X: char, const Y: char> FreeGen for Mono2<X, Y, $I> {}

        impl<const X: char, const Y: char> PolyGen for Mono2<X, Y, $I> {
            type Degree = MDegree<$I>;

            fn degree(&self) -> Self::Degree {
                MDegree::new(map!{ 0 => self.0, 1 => self.1 })
            }

            fn is_unit(&self) -> bool { 
                self.0.is_add_unit() && self.1.is_add_unit()
            }

            fn inv(&self) -> Option<Self> { // (x^i)^{-1} = x^{-i}
                if let (Some(d0), Some(d1)) = (self.0.add_inv(), self.1.add_inv()) {
                    Some(Self(d0, d1))
                } else { 
                    None
                }
            }
        }
    };
}

impl_poly_gen!(usize);
impl_poly_gen!(isize);

impl<const X: char, const Y: char, I> Mono2<X, Y, I> {
    pub fn eval<R>(&self, x: &R, y: &R) -> R
    where R: Mul<Output = R>, for<'x, 'y> &'x R: Pow<&'y I, Output = R> {
        x.pow(&self.0) * y.pow(&self.1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
 
    #[test]
    fn init() { 
        type M = Mono2<'X', 'Y', usize>;
        let d = M::new(2, 3);

        assert_eq!(d.0, 2);
        assert_eq!(d.1, 3);
    }

    #[test]
    fn from_str() { 
        type M = Mono2<'X', 'Y', usize>;
        
        assert_eq!(M::from_str("1"), Ok(M::one()));
        assert_eq!(M::from_str("X"), Ok(M::new(1, 0)));
        assert_eq!(M::from_str("Y"), Ok(M::new(0, 1)));
        assert_eq!(M::from_str("2"), Err(()));
    }

    #[test]
    fn display() { 
        type M = Mono2<'X', 'Y', usize>;
        let d = M::new(0, 0);
        assert_eq!(&d.to_string(), "1");

        let d = M::new(1, 0);
        assert_eq!(&d.to_string(), "X");

        let d = M::new(2, 0);
        assert_eq!(&d.to_string(), "X²");

        let d = M::new(0, 1);
        assert_eq!(&d.to_string(), "Y");

        let d = M::new(0, 2);
        assert_eq!(&d.to_string(), "Y²");

        let d = M::new(1, 1);
        assert_eq!(&d.to_string(), "XY");

        let d = M::new(2, 3);
        assert_eq!(&d.to_string(), "X²Y³");
    }

    #[test]
    fn add_inv_unsigned() { 
        type M = Mono2<'X', 'Y', usize>;

        let d = M::new(0, 0);
        assert_eq!(d.inv(), Some(M::new(0, 0)));

        let d = M::new(1, 0);
        assert_eq!(d.inv(), None);

        let d = M::new(0, 1);
        assert_eq!(d.inv(), None);

        let d = M::new(1, 1);
        assert_eq!(d.inv(), None);
    }

    #[test]
    fn add_inv_signed() { 
        type M = Mono2<'X', 'Y', isize>;

        let d = M::new(0, 0);
        assert_eq!(d.inv(), Some(M::new(0, 0)));

        let d = M::new(1, 0);
        assert_eq!(d.inv(), Some(M::new(-1, 0)));

        let d = M::new(0, 1);
        assert_eq!(d.inv(), Some(M::new(0, -1)));

        let d = M::new(2, 3);
        assert_eq!(d.inv(), Some(M::new(-2, -3)));
    }

    #[test]
    fn eval() { 
        type M = Mono2<'X', 'Y', usize>;

        let d = M::new(0, 0);
        assert_eq!(d.eval::<i32>(&2, &3), 1);

        let d = M::new(1, 0);
        assert_eq!(d.eval::<i32>(&2, &3), 2);

        let d = M::new(0, 1);
        assert_eq!(d.eval::<i32>(&2, &3), 3);

        let d = M::new(1, 1);
        assert_eq!(d.eval::<i32>(&2, &3), 6);

        let d = M::new(2, 3);
        assert_eq!(d.eval::<i32>(&2, &3), 108);
    }
}