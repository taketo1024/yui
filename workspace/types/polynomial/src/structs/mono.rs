use std::fmt::Display;
use std::ops::{AddAssign, Mul, MulAssign};
use std::str::FromStr;
use num_traits::{Zero, One};
use auto_impl_ops::auto_ops;

use yui_core::Elem;
use yui_lin_comb::FreeGen;
use yui_utils::{subscript, superscript};

use crate::{PolyDeg, PolyGen, MDegree};

// `Mono<X, I>` : a struct representing X^d (univar) or ΠX_i^{d_i} (multivar).
// `I` is one of `usize`, `isize`, `MDegree<usize>`, `MDegree<isize>`.

#[derive(Clone, PartialEq, Eq, Hash, Default, Debug, PartialOrd, Ord)]
pub struct Mono<const X: char, I>(I);

impl<const X: char, I> From<I> for Mono<X, I> {
    fn from(d: I) -> Self {
        Self(d)
    }
}

impl<const X: char, I> One for Mono<X, I>
where I: for<'x >AddAssign<&'x I> + Zero {
    fn one() -> Self {
        Self::from(I::zero()) // x^0 = 1.
    }
}

#[auto_ops]
impl<const X: char, I> MulAssign<&Mono<X, I>> for Mono<X, I>
where I: for<'x >AddAssign<&'x I> {
    fn mul_assign(&mut self, rhs: &Mono<X, I>) {
        self.0 += &rhs.0 // x^i * x^j = x^{i+j}
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

// impls for univar-type.
macro_rules! impl_mono_univar {
    ($I:ty) => {
        impl<const X: char> Elem for Mono<X, $I> { 
            fn set_symbol() -> String {
                format!("{X}")
            }
        }

        impl<const X: char> Display for Mono<X, $I> { 
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let x = X.to_string();
                let d = self.degree() as isize;
                write!(f, "{}", fmt_mono(x, d))
            }
        }

        impl<const X: char> FromStr for Mono<X, $I> {
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

impl_mono_univar!(usize);
impl_mono_univar!(isize);

// impls for multivar-type.
macro_rules! impl_mono_multivar {
    ($I:ty) => {
        impl<const X: char> Elem for Mono<X, MDegree<$I>> { 
            fn set_symbol() -> String {
                format!("{X}")
            }
        }

        impl<const X: char> Display for Mono<X, MDegree<$I>> { 
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let d = self.degree();
                let s = d.iter().map(|(&i, &d_i)| {
                    let x = format!("{X}{}", subscript(i as isize));
                    let d = d_i as isize;
                    fmt_mono(x, d)
                }).filter(|s| s != "1").collect::<Vec<_>>().join("");
                if s.is_empty() { 
                    write!(f, "1")
                } else { 
                    write!(f, "{s}")
                }
            }
        }

        impl<const X: char> FromStr for Mono<X, MDegree<$I>> {
            type Err = ();
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                if s == "1" { 
                    return Ok(Mono::one())
                }

                let r = regex::Regex::new(&format!("^{X}([0-9]+)$")).unwrap();
                if let Some(m) = r.captures(&s) { 
                    let i = usize::from_str(&m[1]).unwrap();
                    let mdeg = MDegree::from((i, 1));
                    return Ok(Mono(mdeg))
                }

                // TODO support more complex format. 
                Err(())
            }
        }
    };
}

impl_mono_multivar!(usize);
impl_mono_multivar!(isize);

macro_rules! impl_poly_gen {
    ($I:ty) => {
        impl<const X: char> FreeGen for Mono<X, $I> {}

        impl<const X: char> PolyGen for Mono<X, $I> {
            type Degree = $I;

            fn degree(&self) -> Self::Degree {
                self.0.clone()
            }

            fn is_unit(&self) -> bool { 
                self.degree().is_add_unit()
            }

            fn inv(&self) -> Option<Self> { // (x^i)^{-1} = x^{-i}
                if let Some(i) = self.degree().add_inv() { 
                    Some(Self(i))
                } else { 
                    None
                }
            }
        }
    };
}

impl_poly_gen!(usize);
impl_poly_gen!(isize);
impl_poly_gen!(MDegree<usize>);
impl_poly_gen!(MDegree<isize>);