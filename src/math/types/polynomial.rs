
use std::collections::BTreeMap;
use std::fmt::Display;
use std::iter::{Sum, Product};
use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Neg};
use itertools::Itertools;
use num_traits::{Zero, One};
use auto_impl_ops::auto_ops;

use crate::math::traits::{Ring, RingOps, AlgBase, AddMon, AddMonOps, AddGrpOps, MonOps, AddGrp, Mon};
use super::lin_comb::{LinComb, FreeGenerator};
use crate::utils::format::{subscript, superscript};

pub type Poly  <const X: char, R> = PolyBase<Mono<X, usize>, R>;          // univar
pub type LPoly <const X: char, R> = PolyBase<Mono<X, isize>, R>;          // univar, Laurent
pub type MPoly <const X: char, R> = PolyBase<Mono<X, MDegree<usize>>, R>; // multivar
pub type MLPoly<const X: char, R> = PolyBase<Mono<X, MDegree<isize>>, R>; // multivar, Laurent

pub trait PolyGen: 
    Mul<Output = Self> + 
    One + 
    PartialOrd + 
    Ord + 
    From<Self::Degree> +
    FreeGenerator
{
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

// Display for univar-type.
macro_rules! impl_display_u {
    ($struct:ident<X, $I:ty>) => {
        impl<const X: char> Display for $struct<X, $I> { 
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

impl_display_u!(Mono<X, usize>);
impl_display_u!(Mono<X, isize>);

// Display for multivar-type.
macro_rules! impl_display_m {
    ($struct:ident<X, $I:ty>) => {
        impl<const X: char> Display for $struct<X, $I> { 
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
    };
}

impl_display_m!(Mono<X, MDegree<usize>>);
impl_display_m!(Mono<X, MDegree<isize>>);

macro_rules! impl_poly_gen {
    ($struct:ident, $I:ty) => {
        impl<const X: char> AlgBase for $struct<X, $I> { 
            fn symbol() -> String {
                format!("{}", X)
            }
        }

        impl<const X: char> FreeGenerator for $struct<X, $I> {}

        impl<const X: char> PolyGen for $struct<X, $I> {
            type Degree = $I;
            fn degree(&self) -> Self::Degree {
                self.0.clone()
            }
        }
    };
}

impl_poly_gen!(Mono, usize);
impl_poly_gen!(Mono, isize);
impl_poly_gen!(Mono, MDegree<usize>);
impl_poly_gen!(Mono, MDegree<isize>);

#[derive(Clone, PartialEq, Eq, Default, Debug)]
pub struct PolyBase<I, R>
where 
    I: PolyGen, 
    R: Ring, for<'x> &'x R: RingOps<R>
{
    data: LinComb<I, R>,
    zero: (I, R)
}

impl<I, R> PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(data: LinComb<I, R>) -> Self { 
        Self { data, zero: (I::one(), R::zero()) }
    }

    pub fn from_deg(data: Vec<(I::Degree, R)>) -> Self {
        let data = data.into_iter().map(|(i, r)| 
            (I::from(i), r)
        ).collect_vec();
        Self::from(data)
    }
    
    pub fn as_lincomb(&self) -> &LinComb<I, R> { 
        &self.data
    }

    pub fn len(&self) -> usize { 
        self.data.len()
    }

    pub fn coeff(&self, i: &I) -> &R {
        self.data.coeff(i)
    }

    pub fn coeff_for(&self, i: I::Degree) -> &R {
        self.data.coeff(&I::from(i))
    }

    pub fn iter(&self) -> impl Iterator<Item = (&I, &R)> {
        self.data.iter()
    }

    pub fn into_iter(self) -> impl Iterator<Item = (I, R)> {
        self.data.into_iter()
    }

    pub fn reduce(&mut self) { 
        self.data.reduce()
    }

    pub fn reduced(&self) -> Self { 
        Self::new(self.data.reduced())
    }

    pub fn is_const(&self) -> bool { 
        self.iter().all(|(i, r)| 
            i.is_one() || !i.is_one() && r.is_zero()
        )
    }

    pub fn const_term(&self) -> &R { 
        self.coeff(&I::one())
    }

    pub fn lead_term(&self) -> (&I, &R) { 
        self.iter()
            .filter(|(_, r)| !r.is_zero())
            .max_by_key(|(i, _)| *i)
            .unwrap_or((&self.zero.0, &self.zero.1))
    }

    pub fn lead_coeff(&self) -> &R { 
        self.lead_term().1
    }
}

impl<I, J, R> PolyBase<I, R>
where I: PolyGen<Degree = MDegree<J>>, J: Zero, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_mdeg(data: Vec<(Vec<J>, R)>) -> Self {
        let data = data.into_iter().map(|(v, r)| {
            let mdeg = MDegree::from_vec(v);
            (I::from(mdeg), r)
        }).collect_vec();
        Self::from(data)
    }
}

impl<I, R> From<R> for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(a: R) -> Self {
        Self::new(LinComb::from((I::one(), a)))
    }
}

impl<I, R> From<Vec<(I, R)>> for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(data: Vec<(I, R)>) -> Self {
        Self::new(LinComb::from(data))
    }
}

impl<I, R> Display for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data.fmt(f)
    }
}

impl<I, R> Zero for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn zero() -> Self {
        Self::new(LinComb::zero())
    }

    fn is_zero(&self) -> bool {
        self.is_const() && self.const_term().is_zero()
    }
}

impl<I, R> One for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn one() -> Self {
        Self::new(LinComb::from((I::one(), R::one())))
    }

    fn is_one(&self) -> bool {
        self.is_const() && self.const_term().is_one()
    }
}

impl<I, R> Neg for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self::new(-self.data)
    }
}

impl<I, R> Neg for &PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = PolyBase<I, R>;
    fn neg(self) -> Self::Output {
        PolyBase::new(-&self.data)
    }
}

macro_rules! impl_assop {
    ($trait:ident, $method:ident) => {
        #[auto_ops]
        impl<I, R> $trait<&PolyBase<I, R>> for PolyBase<I, R>
        where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
            fn $method(&mut self, rhs: &PolyBase<I, R>) {
                self.data.$method(&rhs.data)
            }
        }
    };
}

impl_assop!(AddAssign, add_assign);
impl_assop!(SubAssign, sub_assign);
impl_assop!(MulAssign, mul_assign);

macro_rules! impl_accum {
    ($trait:ident, $method:ident, $accum_trait:ident, $accum_method:ident, $accum_init:ident) => {
        impl<I, R> $trait for PolyBase<I, R>
        where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
            fn $method<Iter: Iterator<Item = Self>>(iter: Iter) -> Self {
                let mut res = Self::$accum_init();
                for r in iter { Self::$accum_method(&mut res, r) }
                return res;
            }
        }

        impl<'a, I, R> $trait<&'a Self> for PolyBase<I, R>
        where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
            fn $method<Iter: Iterator<Item = &'a Self>>(iter: Iter) -> Self {
                let mut res = Self::$accum_init();
                for r in iter { Self::$accum_method(&mut res, r) }
                return res;
            }
        }
    };
}

impl_accum!(Sum, sum, AddAssign, add_assign, zero);
impl_accum!(Product, product, MulAssign, mul_assign, one);

macro_rules! impl_alg_op {
    ($trait:ident) => {
        impl<I, R> $trait<Self> for PolyBase<I, R>
        where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {}

        impl<I, R> $trait<PolyBase<I, R>> for &PolyBase<I, R>
        where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {}
    };
}

impl_alg_op!(AddMonOps);
impl_alg_op!(AddGrpOps);
impl_alg_op!(MonOps);
impl_alg_op!(RingOps);

impl<I, R> AlgBase for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn symbol() -> String {
        format!("{}[{}]", R::symbol(), I::symbol())
    }
}

impl<I, R> AddMon for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {}

impl<I, R> AddGrp for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {}

impl<I, R> Mon for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {}

// TODO must change impl for Laurent-type. 

impl<I, R> Ring for PolyBase<I, R>
where I: PolyGen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn inv(&self) -> Option<Self> {
        if self.is_const() { 
            let a = self.const_term();
            if let Some(ainv) = a.inv() { 
                return Some(Self::from(ainv))
            }
        }
        None
    }

    fn is_unit(&self) -> bool {
        self.is_const() && self.const_term().is_unit()
    }

    fn normalizing_unit(&self) -> Self {
        Self::from(self.lead_coeff().normalizing_unit())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
 
    #[test]
    fn init_poly() { 
        type P = Poly::<'x', i32>; 
        let f = P::from_deg(vec![(0, 3), (1, 2), (2, 3)]);
        assert_eq!(&f.to_string(), "3 + 2x + 3x²");
    }
 
    #[test]
    fn init_lpoly() { 
        type P = LPoly::<'x', i32>; 
        let f = P::from_deg(vec![(-1, 4), (0, 2), (2, 3)]);
        assert_eq!(&f.to_string(), "4x⁻¹ + 2 + 3x²");
    }

    #[test]
    fn init_mpoly() { 
        type P = MPoly::<'x', i32>; 
        let f = P::from_mdeg(vec![
            (vec![], 3),
            (vec![1], -1),
            (vec![3, 0, 2], 2),
        ]);
        assert_eq!(&f.to_string(), "3 - x₀ + 2x₀³x₂²");
    }

    #[test]
    fn init_mlpoly() { 
        type P = MLPoly::<'x', i32>; 
        let f = P::from_mdeg(vec![
            (vec![], 3),
            (vec![1], 1),
            (vec![-3, 1, 3], 2),
        ]);
        assert_eq!(&f.to_string(), "3 + 2x₀⁻³x₁x₂³ + x₀");
    }

    #[test]
    fn reduce() { 
        type P = Poly::<'x', i32>;
        let mut f = P::from_deg(vec![(0, 0), (1, 1), (2, 0)]);
        assert_eq!(f.len(), 3);

        f.reduce();
        assert_eq!(f.len(), 1);
    }

    #[test]
    fn coeff() { 
        type P = Poly::<'x', i32>;
        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        
        assert_eq!(f.coeff_for(0), &2);
        assert_eq!(f.coeff_for(1), &3);
        assert_eq!(f.coeff_for(2), &-4);
        assert_eq!(f.coeff_for(3), &0);
    }

    #[test]
    fn const_term() { 
        type P = Poly::<'x', i32>;
        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        assert_eq!(f.const_term(), &2);

        let f = P::from_deg(vec![(1, 3), (2, -4)]);
        assert_eq!(f.const_term(), &0);
    }

    #[test]
    fn lead_term() { 
        type P = Poly::<'x', i32>;

        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        let (x, a) = f.lead_term();
        assert_eq!(x.degree(), 2);
        assert_eq!(a, &-4);

        let f = P::zero();
        let (x, a) = f.lead_term();
        assert_eq!(x.degree(), 0);
        assert_eq!(a, &0);
    }

    #[test]
    fn add() { 
        type P = Poly::<'x', i32>;

        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        let g = P::from_deg(vec![(0, -3), (1, -3), (3, 5)]);
        assert_eq!(f + g, P::from_deg(vec![(0, -1), (1, 0), (2, -4), (3, 5)]));
    }

    #[test]
    fn neg() { 
        type P = Poly::<'x', i32>;

        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        assert_eq!(-f, P::from_deg(vec![(0, -2), (1, -3), (2, 4)]));
    }

    #[test]
    fn sub() { 
        type P = Poly::<'x', i32>;

        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        let g = P::from_deg(vec![(0, -3), (1, -3), (3, 5)]);
        assert_eq!(f - g, P::from_deg(vec![(0, 5), (1, 6), (2, -4), (3, -5)]));
    }

    #[test]
    fn mul() { 
        type P = Poly::<'x', i32>;

        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        let g = P::from_deg(vec![(0, -3), (1, -3), (3, 5)]);
        assert_eq!(f * g, P::from_deg(vec![(0, -6), (1, -15), (2, 3), (3, 22), (4, 15), (5, -20)]));
    }

    #[test]
    fn mul_scal() { 
        type P = Poly::<'x', i32>;

        let f = P::from_deg(vec![(0, 2), (1, 3), (2, -4)]);
        let g = P::from(3);
        assert_eq!(f * g, P::from_deg(vec![(0, 6), (1, 9), (2, -12)]));
    }
}