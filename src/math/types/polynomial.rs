
use std::collections::BTreeMap;
use std::fmt::Display;
use std::ops::{Mul, Add, AddAssign};
use itertools::Itertools;
use num_traits::{Zero, One};
use auto_impl_ops::auto_ops;

use crate::math::traits::{Ring, RingOps, AlgBase};
use super::lin_comb::{LinComb, FreeGenerator};
use crate::utils::format::{subscript, superscript};

pub type Poly  <const X: char, R> = PolyBase<UPolyGen <X>, R>; // univar
pub type LPoly <const X: char, R> = PolyBase<ULPolyGen<X>, R>; // univar, Laurent
pub type MPoly <const X: char, R> = PolyBase<MPolyGen <X>, R>; // multivar
pub type MLPoly<const X: char, R> = PolyBase<MLPolyGen<X>, R>; // multivar, Laurent

pub trait PolyGen: FreeGenerator + Mul<Output = Self> + One {
    type Degree;
    fn degree(&self) -> Self::Degree;
}

#[derive(Clone, Default, PartialEq, Eq, Debug, Hash, PartialOrd, Ord)]
pub struct MDegree<Deg>(BTreeMap<usize, Deg>); // x_0^{d_0} ... x_n^{d_n} <--> [ 0 => d_0, ..., n => d_n ]

impl<Deg> MDegree<Deg>
where Deg: Zero {
    pub fn from_vec(degree: Vec<Deg>) -> Self { 
        Self::from_iter(
            degree.into_iter().enumerate().filter(|(_, d)| 
                !d.is_zero()
            )
        )
    }
}

impl<Deg> FromIterator<(usize, Deg)> for MDegree<Deg> {
    fn from_iter<T: IntoIterator<Item = (usize, Deg)>>(iter: T) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<Deg> Zero for MDegree<Deg>
where for<'x> Deg: Clone + AddAssign<&'x Deg> + Zero {
    fn zero() -> Self {
        Self(BTreeMap::new())
    }

    fn is_zero(&self) -> bool {
        self.0.is_empty()
    }
}

#[auto_ops]
impl<Deg> AddAssign<&MDegree<Deg>> for MDegree<Deg>
where for<'x> Deg: Clone + AddAssign<&'x Deg> + Zero {
    fn add_assign(&mut self, rhs: &MDegree<Deg>) {
        let data = &mut self.0;
        for (i, d) in rhs.0.iter() { 
            if let Some(d_i) = data.get_mut(i) { 
                d_i.add_assign(d);
            } else { 
                data.insert(i.clone(), d.clone());
            }
        }
        data.retain(|_, v| !v.is_zero())
    }
}

macro_rules! impl_poly_gen {
    ($struct:ident, $deg_type:ty) => {
        #[derive(Clone, PartialEq, Eq, Hash, Default, Debug, PartialOrd, Ord)]
        pub struct $struct<const X: char>($deg_type);

        impl<const X: char> From<$deg_type> for $struct<X> {
            fn from(d: $deg_type) -> Self {
                Self(d)
            }
        }

        impl<const X: char> AlgBase for $struct<X> { 
            fn symbol() -> String {
                String::default()
            }
        }

        impl<const X: char> FreeGenerator for $struct<X> {}

        impl<const X: char> One for $struct<X> {
            fn one() -> Self {
                Self::from(<$deg_type>::zero()) // x^0 = 1.
            }
        }

        impl<const X: char> Mul for $struct<X> {
            type Output = Self;
            fn mul(self, rhs: Self) -> Self::Output {
                Self(self.0 + rhs.0) // x^i * x^j = x^{i+j}.
            }
        }

        impl<const X: char> PolyGen for $struct<X> {
            type Degree = $deg_type;
            fn degree(&self) -> Self::Degree {
                self.0.clone()
            }
        }

        impl<const X: char, R> From<Vec<($deg_type, R)>> for PolyBase<$struct<X>, R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            fn from(data: Vec<($deg_type, R)>) -> Self {
                let data = data.into_iter().map(|(i, r)| 
                    ($struct::from(i), r)
                ).collect_vec();
                Self::from(data)
            }
        }
    };
}

macro_rules! impl_upoly_gen {
    ($struct:ident, $deg_type:ty) => {
        impl_poly_gen!($struct, $deg_type);

        impl<const X: char> Display for $struct<X> { 
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let d = self.degree();
                if d.is_zero() { 
                    write!(f, "1")
                } else { 
                    let e = if d.is_one() { 
                        String::from("")
                    } else { 
                        superscript(d as isize)
                    };
                    write!(f, "{X}{e}")
                }
            }
        }
    };
}

macro_rules! impl_mpoly_gen {
    ($struct:ident, $mdeg_type:ident<$deg_type:ty>) => {
        impl_poly_gen!($struct, $mdeg_type<$deg_type>);

        impl<const X: char> Display for $struct<X> { 
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let d = self.degree();
                if d.is_zero() { 
                    write!(f, "1")
                } else { 
                    for (i, d_i) in d.0 { 
                        let e = if d_i.is_one() { 
                            String::from("")
                        } else { 
                            superscript(d_i as isize)
                        };
                        write!(f, "{X}{}{}", subscript(i as isize), e)?;
                    }
                    write!(f, "")
                }
            }
        }

        impl<const X: char, R> From<Vec<(Vec<$deg_type>, R)>> for PolyBase<$struct<X>, R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            fn from(data: Vec<(Vec<$deg_type>, R)>) -> Self {
                let data = data.into_iter().map(|(v, r)| {
                    let mdeg = MDegree::from_vec(v);
                    ($struct::from(mdeg), r)
                }).collect_vec();
                Self::from(data)
            }
        }
    };
}

impl_upoly_gen!(UPolyGen,  usize);
impl_upoly_gen!(ULPolyGen, isize);
impl_mpoly_gen!(MPolyGen,  MDegree<usize>);
impl_mpoly_gen!(MLPolyGen, MDegree<isize>);

#[derive(Clone, PartialEq, Eq, Default, Debug, derive_more::Display)]
pub struct PolyBase<I, R>(LinComb<I, R>)
where 
    I: PolyGen, 
    R: Ring, for<'x> &'x R: RingOps<R>;

impl<I, R> From<Vec<(I, R)>> for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(data: Vec<(I, R)>) -> Self {
        Self(LinComb::from(data))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
 
    #[test]
    fn init_poly() { 
        let f = Poly::<'x', _>::from(vec![(0, 3), (1, 2), (2, 3)]);
        dbg!(f.to_string());
    }
 
    #[test]
    fn init_lpoly() { 
        let f = LPoly::<'x', _>::from(vec![(-1, 3), (0, 2), (2, 3)]);
        dbg!(f.to_string());
    }

    #[test]
    fn init_mpoly() { 
        let f = MPoly::<'x', i64>::from(vec![
            (vec![], 3),
            (vec![1, 0, 2], 2),
        ]);
        dbg!(f.to_string());
    }
}