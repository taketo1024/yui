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

#[cfg(test)]
mod tests { 
    use num_traits::Zero;

    use super::*;

    fn mono<I, Itr>(itr: Itr) -> MultiVar<'X', I> 
    where I: Zero, Itr: IntoIterator<Item = (usize, I)> {
        let mdeg = MultiDeg::from_iter(itr);
        MultiVar::from(mdeg)
    }

    #[test]
    fn is_divisible() { 
        let one = mono([]);
        let d1 = mono::<usize, _>([(0, 1), (1, 2), (2, 3)]);
        let d2 = mono::<usize, _>([        (1, 2), (2, 1)]);
        let d3 = mono::<usize, _>([(0, 1), (1, 3)        ]);

        assert!(one.divides(&d1));
        assert!(d1.divides(&d1));
        assert!(d2.divides(&d1));
        assert!(!d1.divides(&d2));
        assert!(!d3.divides(&d1));
        assert!(!d1.divides(&d3));
    }

    #[test]
    fn div() { 
        let one = mono([]);
        let d1 = mono::<usize, _>([(0, 1), (1, 2), (2, 3)]);
        let d2 = mono::<usize, _>([        (1, 2), (2, 1)]);

        assert_eq!(&d1 / &one, d1);
        assert_eq!(&d1 / &d1, mono([]));
        assert_eq!(&d1 / &d2, mono([(0, 1), (2, 2)]));
    }

    #[test]
    fn is_divisible_isize() { 
        let one = mono([]);
        let d1 = mono::<isize, _>([(0, 1), (1, 2), (2, 3)]);
        let d2 = mono::<isize, _>([        (1, 2), (2, 1)]);
        let d3 = mono::<isize, _>([(0, 1), (1, 3)        ]);

        assert!(one.divides(&d1));
        assert!(d1.divides(&d1));
        assert!(d2.divides(&d1));
        assert!(d1.divides(&d2));
        assert!(d3.divides(&d1));
        assert!(d1.divides(&d3));
    }

    #[test]
    fn div_isize() { 
        let one = mono::<isize, _>([]);
        let d1 = mono::<isize, _>([(0, 1), (1, 2), (2, 3)]);
        let d2 = mono::<isize, _>([        (1, 2), (2, 1)]);

        assert_eq!(&d1 / &one, d1);
        assert_eq!(&one / &d1, d1.inv().unwrap());
        assert_eq!(&d1 / &d2, mono([(0, 1), (2, 2)]));
        assert_eq!(&d2 / &d1, mono([(0, -1), (2, -2)]));
    }
}