use std::fmt::{Display, Debug};
use std::ops::{AddAssign, MulAssign, Mul, Div, DivAssign, SubAssign, Add};
use std::str::FromStr;
use std::hash::Hash;
use num_traits::{Zero, One, ToPrimitive, FromPrimitive};
use itertools::Itertools;
use auto_impl_ops::auto_ops;

use crate::{Elem, ElemBase};
use crate::lc::{Gen, OrdForDisplay};
use crate::util::format::subscript;
use super::{Mono, MultiDeg, MonoOrd};
use super::var::{fmt_mono, parse_mono_deg};

#[derive(Clone, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde_with::DeserializeFromStr))]
pub struct MultiVar<const X: char, I> (
    MultiDeg<I>
);

impl<const X: char, I> MultiVar<X, I> {
    pub fn var_symbol() -> char { 
        X
    }

    pub fn deg_for(&self, i: usize) -> I
    where I: Copy {
        self.0[i]
    }

    pub fn total_deg(&self) -> I
    where I: Zero + for<'x> Add<&'x I, Output = I> {
        self.0.total()
    }

    fn to_string_u(&self, unicode: bool) -> String
    where I: ToPrimitive { 
        let seq = self.0.iter().map(|(&i, d)| {
            let x = if unicode { 
                format!("{X}{}", subscript(i))
            } else { 
                format!("{X}_{}", i)
            };
            (x, d)
        });
        fmt_mono_n(seq, unicode)
    }
}

impl<const X: char, I> From<MultiDeg<I>> for MultiVar<X, I> {
    fn from(d: MultiDeg<I>) -> Self {
        Self(d)
    }
}

impl<const X: char, I> From<(usize, I)> for MultiVar<X, I>
where I: Zero {
    fn from(value: (usize, I)) -> Self {
        Self::from_iter([value])
    }
}

impl<const X: char, const N: usize, I> From<[I; N]> for MultiVar<X, I>
where I: Zero {
    fn from(degs: [I; N]) -> Self {
        Self::from(MultiDeg::from(degs))
    }
}

impl<const X: char, I> FromIterator<(usize, I)> for MultiVar<X, I>
where I: Zero {
    fn from_iter<T: IntoIterator<Item = (usize, I)>>(iter: T) -> Self {
        Self::from(MultiDeg::from_iter(iter))
    }
}

impl<const X: char, I> FromStr for MultiVar<X, I>
where I: Zero + FromStr + FromPrimitive {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use regex::Regex;

        if s == "1" { 
            return Ok(MultiVar::from(MultiDeg::empty()))
        }

        // TODO must support braced indices.
        
        let p = format!(r"({X}_([0-9]+))(\^\{{?-?[0-9]+\}}?)?");
        let p_all = format!(r"^({p}\s?)+$");

        if !Regex::new(&p_all).unwrap().is_match(&s) { 
            return Err(format!("Failed to parse: {s}"))
        }

        let r = Regex::new(&p).unwrap();
        let mut degs = vec![];
        
        for c in r.captures_iter(&s) {
            let x = &c[1];
            let i = usize::from_str(&c[2]).map_err(|e| e.to_string())?;
            if let Some(d) = parse_mono_deg(x, &c[0]) { 
                degs.push((i, d));
            }
        };

        let mvar = MultiVar::from_iter(degs);
        Ok(mvar)
    }
}

#[auto_ops]
impl<const X: char, I> MulAssign<&MultiVar<X, I>> for MultiVar<X, I>
where I: Zero + for<'x> AddAssign<&'x I> {
    fn mul_assign(&mut self, rhs: &MultiVar<X, I>) {
        self.0 += &rhs.0 // x^i * x^j = x^{i+j}
    }
}

#[auto_ops]
impl<const X: char, I> DivAssign<&MultiVar<X, I>> for MultiVar<X, I>
where I: Zero + for<'x> SubAssign<&'x I> {
    fn div_assign(&mut self, rhs: &MultiVar<X, I>) {
        self.0 -= &rhs.0 // x^i * x^j = x^{i+j}
    }
}

impl<const X: char, I> One for MultiVar<X, I>
where I: Zero + for<'x> AddAssign<&'x I> {
    fn one() -> Self {
        Self::from(MultiDeg::zero()) // x^0 = 1.
    }
}

impl<const X: char, I> Display for MultiVar<X, I>
where I: ToPrimitive { 
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = self.to_string_u(true);
        f.write_str(&s)
    }
}

impl<const X: char, I> Debug for MultiVar<X, I>
where I: ToPrimitive { 
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Display::fmt(self, f)
    }
}

impl<const X: char, I> MonoOrd for MultiVar<X, I>
where I: Zero + Ord + for<'x> Add<&'x I, Output = I> {
    fn cmp_lex(&self, other: &Self) -> std::cmp::Ordering { 
        MultiDeg::cmp_lex(&self.0, &other.0)
    }

    fn cmp_grlex(&self, other: &Self) -> std::cmp::Ordering { 
        MultiDeg::cmp_grlex(&self.0, &other.0)
    }
}

impl<const X: char, I> OrdForDisplay for MultiVar<X, I>
where I: Zero + Ord + for<'x> Add<&'x I, Output = I> {
    fn cmp_for_display(&self, other: &Self) -> std::cmp::Ordering {
        Self::cmp_lex(self, other)
    }
}


#[cfg(feature = "serde")]
impl<const X: char, I> serde::Serialize for MultiVar<X, I>
where I: ToPrimitive {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where S: serde::Serializer {
        serializer.serialize_str(&self.to_string_u(false))
    }
}

impl<const X: char, I> Elem for MultiVar<X, I>
where I: ElemBase + ToPrimitive { 
    fn math_symbol() -> String {
        format!("{X}")
    }
}

impl<const X: char, I> Gen for MultiVar<X, I>
where I: ElemBase + Zero + Ord + Hash + ToPrimitive + for<'x> Add<&'x I, Output = I> {}

macro_rules! impl_multivar_unsigned {
    ($I:ty) => {
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

pub(crate) fn fmt_mono_n<'a, X, I, S>(seq: S, unicode: bool) -> String
where X: ToString, I: 'a + ToPrimitive, S: IntoIterator<Item = (X, &'a I)> { 
    let s = seq.into_iter().map(|(x, d)| {
        let m = fmt_mono(&x.to_string(), d, unicode);
        if m == "1" { 
            "".to_string()
        } else {
            m
        }
    }).join("");

    if s.is_empty() { 
        "1".to_string()
    } else { 
        s
    }
}

#[cfg(feature = "tex")]
mod tex {
    use crate::tex::TeX;
    use super::*;

    impl<const X: char, I> TeX for MultiVar<X, I>
    where I: ToPrimitive {
        fn tex_math_symbol() -> String { 
            format!("{X}_1,\\ldots")
        }
        fn tex_string(&self) -> String {
            self.to_string_u(false)
        }
    }
}

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

        let s = "X_0^{-1}";
        assert_eq!(M::from_str(s), Ok(M::from((0, -1))));

        let s = "X_0^{-1}X_2^{12}";
        assert_eq!(M::from_str(s), Ok(M::from_iter([(0, -1), (2, 12)])));

        let s = "2";
        assert!(M::from_str(s).is_err());

        let s = "Y_0";
        assert!(M::from_str(s).is_err());
    }

    #[test]
    #[cfg(feature = "serde")]
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
        
        assert_eq!(&ser, "\"X_0^{-1}X_2^3\"");
        assert_eq!(d, des);
    }
}