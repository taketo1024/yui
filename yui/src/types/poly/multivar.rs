use std::fmt::{Display, Debug};
use std::str::FromStr;
use num_traits::{Zero, One};
use itertools::Itertools;

use crate::Elem;
use crate::lc::Gen;
use crate::util::format::subscript;
use super::{Mono, MultiDeg, Univar};
use super::univar::{fmt_mono, parse_mono};

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

            fn fmt_impl(&self, unicode: bool) -> String { 
                let deg = self.deg();
                let s = deg.iter().map(|(&i, &d)| {
                    let x = if unicode { 
                        format!("{X}{}", subscript(i))
                    } else { 
                        format!("{X}_{}", i)
                    };
                    fmt_mono(&x, d, unicode)
                }).join("");
        
                if s.is_empty() { 
                    "1".to_string()
                } else { 
                    s
                }
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
                let s = self.fmt_impl(true);
                f.write_str(&s)
            }
        }

        impl<const X: char> Debug for MultiVar<X, $I> { 
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                Display::fmt(self, f)
            }
        }

        impl<const X: char> FromStr for MultiVar<X, $I> {
            type Err = String;
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                use regex::Regex;

                if s == "1" { 
                    return Ok(Univar::one())
                }

                // TODO validate s
        
                let p = format!(r"({X}_([0-9]+))(\^-?[0-9]+)?");
                let p_all = format!(r"^({p}\s?)+$");

                if !Regex::new(&p_all).unwrap().is_match(&s) { 
                    return Err(format!("Failed to parse: {s}"))
                }

                let r = Regex::new(&p).unwrap();
                let mut degs = vec![];
                
                for c in r.captures_iter(&s) {
                    let x = &c[1];
                    let i = usize::from_str(&c[2]).map_err(|e| e.to_string())?;
                    let d = parse_mono(x, &c[0])?;
                    degs.push((i, d));
                };
        
                let mvar = MultiVar::from_iter(degs);
                Ok(mvar)
            }
        }

        impl<const X: char> serde::Serialize for MultiVar<X, $I> {
            fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
            where S: serde::Serializer {
                serializer.serialize_str(&self.fmt_impl(false))
            }
        }
        
        impl<'de, const X: char> serde::Deserialize<'de> for MultiVar<X, $I> {
            fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
            where D: serde::Deserializer<'de> {
                let s = String::deserialize(deserializer)?;
                Self::from_str(&s).map_err(serde::de::Error::custom)
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
            pub fn generate(n: usize, tot_deg: usize) -> impl Iterator<Item = Self> { 
                use dinglebit_combinatorics::Combination as C;
                use crate::algo::rep_comb;
        
                // MEMO: 
                // Mathematically this restriction is unnecessary, 
                // but (n, tot_deg) = (0, 0) will fail. 
                assert!(n > 0);
        
                let c = C::new(n + tot_deg - 1, n - 1);
                c.into_iter().map(move |mut list| {
                    list.push(n + tot_deg - 1); // the right-end wall
                    Self::from_iter( rep_comb(&list) )
                })
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
    use itertools::Itertools;
    use num_integer::binomial;

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
        let mons = M::generate(n, tot).collect_vec();

        assert_eq!(mons.len(), binomial(n + tot - 1, n - 1));
        assert!(mons.iter().all(|x| x.total_deg() == tot));
        assert!(mons.iter().all_unique());

        assert_eq!(M::generate(n, 0).collect_vec(), vec![M::one()]); // deg(1) = 0, using n-vars.
    }

    #[test]
    fn from_str() { 
        type M = MultiVar<'X', isize>;

        let s = "1";
        assert_eq!(M::from_str(s), Ok(M::one()));

        let s = "X_0";
        assert_eq!(M::from_str(s), Ok(M::from((0, 1))));

        let s = "X_1";
        assert_eq!(M::from_str(s), Ok(M::from((1, 1))));

        let s = "X_1^2";
        assert_eq!(M::from_str(s), Ok(M::from((1, 2))));

        let s = "X_0X_2^3";
        assert_eq!(M::from_str(s), Ok(M::from_iter([(0, 1), (2, 3)])));

        let s = "X_0^-1";
        assert_eq!(M::from_str(s), Ok(M::from((0, -1))));

        let s = "X_0^-1X_2^4";
        assert_eq!(M::from_str(s), Ok(M::from_iter([(0, -1), (2, 4)])));

        let s = "2";
        assert!(M::from_str(s).is_err());

        let s = "Y_0";
        assert!(M::from_str(s).is_err());
    }

    #[test]
    fn serialize() { 
        type M = MultiVar<'X', isize>;

        let d = M::from([]);
        let ser = serde_json::to_string(&d).unwrap();
        let des = serde_json::from_str::<M>(&ser).unwrap();
        
        assert_eq!(&ser, "\"1\"");
        assert_eq!(d, des);

        let d = M::from([1]);
        let ser = serde_json::to_string(&d).unwrap();
        let des = serde_json::from_str::<M>(&ser).unwrap();

        assert_eq!(&ser, "\"X_0\"");
        assert_eq!(d, des);

        let d = M::from([2]);
        let ser = serde_json::to_string(&d).unwrap();
        let des = serde_json::from_str::<M>(&ser).unwrap();
        
        assert_eq!(&ser, "\"X_0^2\"");
        assert_eq!(d, des);

        let d = M::from([-1, 0, 3]);
        let ser = serde_json::to_string(&d).unwrap();
        let des = serde_json::from_str::<M>(&ser).unwrap();
        
        assert_eq!(&ser, "\"X_0^-1X_2^3\"");
        assert_eq!(d, des);
    }
}