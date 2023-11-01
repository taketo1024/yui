use std::fmt::Display;
use std::ops::{AddAssign, Mul, MulAssign, DivAssign, SubAssign, Div};
use std::str::FromStr;
use num_traits::{Zero, One};
use auto_impl_ops::auto_ops;

use yui_core::Elem;
use yui_lin_comb::Gen;
use yui_utils::superscript;

use crate::{VarDeg, PolyGen};

// `Univar<X, I>` : a struct representing X^d (univar) or ΠX_i^{d_i} (multivar).
// `I` is one of `usize`, `isize`, `MultiDeg<usize>`, `MultiDeg<isize>`.

#[derive(Clone, PartialEq, Eq, Hash, Default, Debug, PartialOrd, Ord)]
pub struct Univar<const X: char, I>(pub(crate) I);

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

pub(crate) fn fmt_mono(x: String, d: isize) -> String { 
    if d.is_zero() { 
        "1".to_string()
    } else if d.is_one() { 
        x
    } else {
        let e = superscript(d); 
        format!("{x}{e}")
    }
}

macro_rules! impl_univar {
    ($I:ty) => {
        impl<const X: char> Elem for Univar<X, $I> { 
            fn math_symbol() -> String {
                format!("{X}")
            }
        }

        impl<const X: char> Display for Univar<X, $I> { 
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let x = X.to_string();
                let d = self.deg() as isize;
                write!(f, "{}", fmt_mono(x, d))
            }
        }

        impl<const X: char> FromStr for Univar<X, $I> {
            type Err = ();
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                if s.len() == 1 { 
                    let x = s.chars().next().unwrap();
                    if x == '1' { 
                        return Ok(Self(0))
                    } else if x == X {
                        return Ok(Self(1))
                    }
                }

                // TODO support more complex format. 
                Err(())
            }
        }
    };
}

impl_univar!(usize);
impl_univar!(isize);

macro_rules! impl_poly_gen {
    ($I:ty) => {
        impl<const X: char> Gen for Univar<X, $I> {}

        impl<const X: char> PolyGen for Univar<X, $I> {
            type Deg = $I;

            fn deg(&self) -> Self::Deg {
                self.0.clone()
            }

            fn is_unit(&self) -> bool { 
                self.deg().is_negatable()
            }

            fn inv(&self) -> Option<Self> { // (x^i)^{-1} = x^{-i}
                if let Some(i) = self.deg().neg_opt() { 
                    Some(Self(i))
                } else { 
                    None
                }
            }

            fn divides(&self, other: &Self) -> bool { 
                self.0.strict_leq(&other.0)
            }
        }
    };
}

impl_poly_gen!(usize);
impl_poly_gen!(isize);

pub(crate) use impl_poly_gen;