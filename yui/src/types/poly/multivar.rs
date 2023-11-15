use std::fmt::{Display, Debug};
use std::str::FromStr;
use num_traits::{Zero, One};

use crate::Elem;
use crate::lc::Gen;
use crate::format::subscript;
use super::{Mono, MultiDeg, Univar, fmt_mono};

pub type MultiVar<const X: char, I> = Univar<X, MultiDeg<I>>;

// impls for multivar-type.
macro_rules! impl_multivar {
    ($I:ty) => {
        impl<const X: char> MultiVar<X, $I> {
            pub fn deg_for(&self, i: usize) -> $I {
                self.0[i]
            }

            pub fn total_deg(&self) -> $I {
                self.0.total()
            }
        }

        impl<const X: char> From<(usize, $I)> for MultiVar<X, $I> {
            fn from(value: (usize, $I)) -> Self {
                Self::from_iter([value])
            }
        }
        
        impl<const X: char, const N: usize> From<[$I; N]> for MultiVar<X, $I> {
            fn from(degs: [$I; N]) -> Self {
                Self::from(MultiDeg::from(degs))
            }
        }

        impl<const X: char> FromIterator<(usize, $I)> for MultiVar<X, $I> {
            fn from_iter<T: IntoIterator<Item = (usize, $I)>>(iter: T) -> Self {
                Self::from(MultiDeg::from_iter(iter))
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

        impl<const X: char> Debug for MultiVar<X, $I> { 
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                Display::fmt(self, f)
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
                    let m = MultiVar::from(mdeg);
                    return Ok(m)
                }

                // TODO support more complex format. 
                Err(())
            }
        }

        impl<const X: char> Elem for MultiVar<X, $I> { 
            fn math_symbol() -> String {
                format!("{X}")
            }
        }

        impl<const X: char> Gen for MultiVar<X, $I> {}
    };
}

macro_rules! impl_multivar_unsigned {
    ($I:ty) => {
        impl_multivar!($I);

        impl<const X: char> Mono for MultiVar<X, $I> {
            type Deg = MultiDeg<$I>;

            fn deg(&self) -> Self::Deg {
                self.0.clone()
            }

            fn is_unit(&self) -> bool { 
                self.0.is_zero()
            }

            fn inv(&self) -> Option<Self> {
                if self.is_unit() { 
                    Some(Self(MultiDeg::zero()))
                } else { 
                    None
                }
            }

            fn divides(&self, other: &Self) -> bool { 
                self.0.all_leq(&other.0)
            }
        }

        impl<const X: char> MultiVar<X, $I> {
            pub fn generate(n: usize, tot_deg: usize) -> Vec<Self> { 
                type Degs = Vec<(usize, usize)>;
                fn generate(n: usize, tot_deg: usize, i: usize, res: &mut Vec<Degs>, prev: Degs) {
                    if i < n - 1 { 
                        for d_i in (0..=tot_deg).rev() { 
                            let mut curr = prev.clone();
                            curr.push((i, d_i));
                            
                            let rem = tot_deg - d_i;
                            generate(n, rem, i + 1, res, curr)
                        }
                    } else { 
                        let mut curr = prev;
                        curr.push((i, tot_deg));
                        res.push(curr);
                    }
                }
            
                let mut res = vec![];
                generate(n, tot_deg, 0, &mut res, vec![]);
            
                res.into_iter().map(|d| {
                    Self::from_iter(d)
                }).collect()
            }
        }
    };
}

macro_rules! impl_multivar_signed {
    ($I:ty) => {
        impl_multivar!($I);

        impl<const X: char> Mono for MultiVar<X, $I> {
            type Deg = MultiDeg<$I>;

            fn deg(&self) -> Self::Deg {
                self.0.clone()
            }

            fn is_unit(&self) -> bool { 
                true
            }

            fn inv(&self) -> Option<Self> { // (x^i)^{-1} = x^{-i}
                Some(Self(-&self.0))
            }

            fn divides(&self, _other: &Self) -> bool { 
                true
            }
        }
    };
}

impl_multivar_unsigned!(usize);
impl_multivar_signed!  (isize);

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn display() { 
        type M = MultiVar<'X', usize>;

        assert_eq!(format!("{}", M::from([])), "1");
        assert_eq!(format!("{}", M::from([0])), "1");
        assert_eq!(format!("{}", M::from([1])), "X₀");
        assert_eq!(format!("{}", M::from([3])), "X₀³");
        assert_eq!(format!("{}", M::from([1,0,3])), "X₀X₂³");
    }

    #[test]
    fn from_pair() { 
        type M = MultiVar<'X', usize>;

        let d = M::from((2, 3)); // X₂³
        assert_eq!(d.0, MultiDeg::from_iter([(2, 3)]));
    }

    #[test]
    fn from_arr() { 
        type M = MultiVar<'X', usize>;

        let d = M::from([1,0,3]); // X₀X₂³
        assert_eq!(d.0, MultiDeg::from_iter([(0, 1), (2, 3)]));
    }

    #[test]
    fn from_iter() { 
        type M = MultiVar<'X', usize>;

        let d = M::from_iter([(0, 1), (2, 3)]); // X₀X₂³
        assert_eq!(d.0, MultiDeg::from_iter([(0, 1), (2, 3)]));
    }

    #[test]
    fn deg_for() { 
        type M = MultiVar<'X', usize>;

        let d = M::from([1,0,3]); // X₀X₂³
        assert_eq!(d.deg_for(0), 1);
        assert_eq!(d.deg_for(1), 0);
        assert_eq!(d.deg_for(2), 3);
        assert_eq!(d.deg_for(3), 0);
    }

    #[test]
    fn total_deg() { 
        type M = MultiVar<'X', usize>;

        let d = M::from([1,0,3]); // X₀X₂³
        assert_eq!(d.total_deg(), 4);
    }

    #[test]
    fn is_divisible() { 
        type M = MultiVar<'X', usize>;

        let one = M::from([]);
        let d1 = M::from([1,2,3]);
        let d2 = M::from([0,2,1]);
        let d3 = M::from([1,3,0]);

        assert!(one.divides(&d1));
        assert!(d1.divides(&d1));
        assert!(d2.divides(&d1));
        assert!(!d1.divides(&d2));
        assert!(!d3.divides(&d1));
        assert!(!d1.divides(&d3));
    }

    #[test]
    fn div() { 
        type M = MultiVar<'X', usize>;

        let one = M::from([]);
        let d1 = M::from([1,2,3]);
        let d2 = M::from([0,2,1]);

        assert_eq!(&d1 / &one, d1);
        assert_eq!(&d1 / &d1, M::from([]));
        assert_eq!(&d1 / &d2, M::from([1,0,2]));
    }

    #[test]
    fn is_divisible_isize() { 
        type M = MultiVar<'X', isize>;

        let one = M::from([]);
        let d1 = M::from([1,2,3]);
        let d2 = M::from([0,2,1]);
        let d3 = M::from([1,3]);

        assert!(one.divides(&d1));
        assert!(d1.divides(&d1));
        assert!(d2.divides(&d1));
        assert!(d1.divides(&d2));
        assert!(d3.divides(&d1));
        assert!(d1.divides(&d3));
    }

    #[test]
    fn div_isize() { 
        type M = MultiVar<'X', isize>;

        let one = M::from([]);
        let d1 = M::from([1,2,3]);
        let d2 = M::from([0,2,1]);

        assert_eq!(&d1 / &one, d1);
        assert_eq!(&one / &d1, d1.inv().unwrap());
        assert_eq!(&d1 / &d2, M::from([1,0,2]));
        assert_eq!(&d2 / &d1, M::from([-1,0,-2]));
    }

    #[test]
    fn gen_mons() { 
        type M = MultiVar<'X', usize>;

        let n = 3;
        let tot = 5;
        let mons = M::generate(n, tot);
        assert_eq!(mons.len(), 21);
    }
}