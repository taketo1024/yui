use std::fmt::Display;
use std::str::FromStr;
use num_traits::One;

use yui_core::Elem;
use yui_lin_comb::Gen;
use yui_utils::subscript;

use crate::{PolyGen, VarDeg, MultiDeg, Univar, fmt_mono};

use super::univar::impl_poly_gen;

pub type MultiVar<const X: char, I> = Univar<X, MultiDeg<I>>;

// impls for multivar-type.
macro_rules! impl_multivar {
    ($I:ty) => {
        impl<const X: char> Elem for MultiVar<X, $I> { 
            fn math_symbol() -> String {
                format!("{X}")
            }
        }

        impl<const X: char> Display for MultiVar<X, $I> { 
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let d = self.deg();
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

        impl<const X: char> FromStr for MultiVar<X, $I> {
            type Err = ();
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                if s == "1" { 
                    return Ok(Univar::one())
                }

                let r = regex::Regex::new(&format!("^{X}([0-9]+)$")).unwrap();
                if let Some(m) = r.captures(&s) { 
                    let i = usize::from_str(&m[1]).unwrap();
                    let mdeg = MultiDeg::from((i, 1));
                    return Ok(Univar(mdeg))
                }

                // TODO support more complex format. 
                Err(())
            }
        }
    };
}

impl_multivar!(usize);
impl_multivar!(isize);

impl_poly_gen!(MultiDeg<usize>);
impl_poly_gen!(MultiDeg<isize>);