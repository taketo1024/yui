use std::collections::HashMap;
use std::fmt::{Display, Debug};
use std::hash::Hash;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, MulAssign};
use itertools::Itertools;
use num_traits::Zero;
use auto_impl_ops::auto_ops;

use crate::utils::collections::map;
use crate::utils::display::OrdForDisplay;

use super::super::traits::{AlgBase, AddMon, AddMonOps, AddGrp, AddGrpOps, Ring, RingOps, RMod, RModOps};
pub trait FreeGenerator: AlgBase + Hash + OrdForDisplay {}

#[derive(PartialEq, Eq, Clone, Debug, Default)]
pub struct LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{ 
    data: HashMap<X, R>,
    r_zero: R
}

impl<X, R> LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{ 
    pub fn new_raw(data: HashMap<X, R>) -> Self {
        Self { data, r_zero: R::zero() }
    }

    pub fn new(data: HashMap<X, R>) -> Self {
        let mut new = Self::new_raw(data);
        new.reduce();
        new
    }

    pub fn wrap(x: X) -> Self { 
        Self::new_raw(map!{ x => R::one() })
    }

    pub fn from_iter<I>(iter: I) -> Self 
    where I: Iterator<Item = (X, R)> {
        let data = iter.collect::<HashMap<_, _>>();
        Self::new(data)
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn keys(&self) -> impl Iterator<Item = &X> {
        self.data.keys()
    }

    pub fn coeff(&self, x: &X) -> &R { 
        self.data.get(x).unwrap_or(&self.r_zero)
    }

    pub fn iter(&self) -> impl Iterator<Item = (&X, &R)> {
        self.data.iter()
    }

    pub fn into_iter(self) -> impl Iterator<Item = (X, R)> {
        self.data.into_iter()
    }

    pub fn reduce(&mut self) { 
        self.data.retain(|_, r| !r.is_zero());
    }

    pub fn reduced(&self) -> Self { 
        let mut copy = self.clone();
        copy.reduce();
        copy
    }

    pub fn map_coeffs<F>(&self, f: F) -> Self 
    where F: Fn(&R) -> R { 
        let data = self.iter().map(|(x, r)| (x.clone(), f(r))).collect();
        Self::new(data)
    }

    pub fn map_coeffs_into<F>(self, f: F) -> Self 
    where F: Fn(R) -> R { 
        let data = self.into_iter().map(|(x, r)| (x, f(r))).collect();
        Self::new(data)
    }

    pub fn apply<F>(&self, f: F) -> Self 
    where F: Fn(&X) -> Vec<(X, R)> {
        let mut res = Self::zero();
        for (x, r) in self.iter() { 
            for (y, s) in f(x) { 
                res += (y, r * &s);
            }
        }
        res
    }
}

impl<X, R> From<(X, R)> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn from(pair: (X, R)) -> Self {
        Self::new(map!{ pair.0 => pair.1 })
    }
}

impl<X, R> From<Vec<(X, R)>> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn from(data: Vec<(X, R)>) -> Self {
        Self::from_iter(data.into_iter())
    }
}

impl<X, R> Display for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_zero() { 
            return write!(f, "0")
        }

        let mut initial = true;
        let sorted = self.iter().sorted_by(|(x, _), (y, _)| x.cmp_for_display(y) );
        for (x, r) in sorted {
            let r = r.to_string();
            let x = x.to_string();

            if initial { 
                if r == "1" { 
                    write!(f, "{x}")?;
                } else if r == "-1" { 
                    write!(f, "-{x}")?;
                } else if x == "1" {
                    write!(f, "{r}")?;
                } else { 
                    write!(f, "{r}{x}")?;
                };
                initial = false
            } else {
                let (sign, r) = if r.starts_with("-") { 
                    ("-", &r[1..]) 
                } else { 
                    ("+", r.as_str())
                };
                if r == "1" { 
                    write!(f, " {sign} {x}")?;
                } else if x == "1" { 
                    write!(f, " {sign} {r}")?;
                } else { 
                    write!(f, " {sign} {r}{x}")?;
                };
            }
        }
        
        Ok(())
    }
}

impl<X, R> Zero for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn zero() -> Self {
        Self::new_raw(HashMap::new())
    }

    fn is_zero(&self) -> bool {
        self.data.values().all(|r| r.is_zero())
    }
}

impl<X, R> Neg for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.map_coeffs_into(|r| -r)
    }
}

impl<X, R> Neg for &LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = LinComb<X, R>;

    fn neg(self) -> Self::Output {
        self.map_coeffs(|r| -r)
    }
}

impl<X, R> AddAssign<(X, R)> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn add_assign(&mut self, rhs: (X, R)) {
        if rhs.1.is_zero() { return }

        let data = &mut self.data;
        let (x, r) = rhs;
        if data.contains_key(&x) { 
            let v = data.get_mut(&x).unwrap();
            v.add_assign(r);
        } else { 
            data.insert(x, r);
        }
    }
}

impl<X, R> AddAssign<(&X, &R)> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn add_assign(&mut self, rhs: (&X, &R)) {
        if rhs.1.is_zero() { return }

        let data = &mut self.data;
        let (x, r) = rhs;
        if data.contains_key(x) { 
            let v = data.get_mut(x).unwrap();
            v.add_assign(r);
        } else { 
            data.insert(x.clone(), r.clone());
        }
    }
}

#[auto_ops]
impl<X, R> AddAssign<&LinComb<X, R>> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn add_assign(&mut self, rhs: &Self) {
        for e in rhs.data.iter() { 
            self.add_assign(e);
        }
        self.reduce()
    }
}

#[auto_ops]
impl<X, R> SubAssign<&LinComb<X, R>> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn sub_assign(&mut self, rhs: &Self) {
        for e in rhs.data.iter() { 
            let e = (e.0, &-e.1);
            self.add_assign(e);
        }
        self.reduce()
    }
}

#[auto_ops]
impl<X, R> MulAssign<&R> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn mul_assign(&mut self, rhs: &R) {
        let data = std::mem::take(&mut self.data);
        self.data = data.into_iter().map(|(x, r)| (x, &r * rhs)).collect();
        self.reduce()
    }
}

#[auto_ops]
impl<X, R> Mul for &LinComb<X, R>
where 
    X: FreeGenerator + Mul<Output = X>,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = LinComb<X, R>;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut res = Self::Output::zero();
        res.data.reserve(self.len() * rhs.len());

        for (x, r) in self.iter().filter(|(_, r)| !r.is_zero()) { 
            for (y, s) in rhs.iter().filter(|(_, r)| !r.is_zero()) { 
                let xy = x.clone() * y.clone();
                let rs = r * s;
                res += (xy, rs);
            }
        }
        
        res.reduce();
        res
    }
}

macro_rules! impl_accum {
    ($trait:ident, $method:ident, $accum_trait:ident, $accum_method:ident, $accum_init:ident) => {
        impl<X, R> $trait<Self> for LinComb<X, R>
        where X: FreeGenerator, R: Ring, for<'x> &'x R: RingOps<R> {
            fn $method<Iter: Iterator<Item = Self>>(iter: Iter) -> Self {
                let mut res = Self::$accum_init();
                for r in iter { Self::$accum_method(&mut res, r) }
                return res;
            }
        }

        impl<'a, X, R> $trait<&'a Self> for LinComb<X, R>
        where X: FreeGenerator, R: Ring, for<'x> &'x R: RingOps<R> {
            fn $method<Iter: Iterator<Item = &'a Self>>(iter: Iter) -> Self {
                let mut res = Self::$accum_init();
                for r in iter { Self::$accum_method(&mut res, r) }
                return res;
            }
        }
    }
}

impl_accum!(Sum, sum, AddAssign, add_assign, zero);

macro_rules! impl_alg_ops {
    ($trait:ident) => {
        impl<X, R> $trait<Self> for LinComb<X, R>
        where X: FreeGenerator, R: Ring, for<'x> &'x R: RingOps<R> {}

        impl<X, R> $trait<LinComb<X, R>> for &LinComb<X, R>
        where X: FreeGenerator, R: Ring, for<'x> &'x R: RingOps<R> {}
    };
}

impl_alg_ops!(AddMonOps);
impl_alg_ops!(AddGrpOps);

impl<X, R> AlgBase for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn symbol() -> String {
        format!("Free<{}; {}>", X::symbol(), R::symbol())
    }
}

impl<X, R> AddMon for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<X, R> AddGrp for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{}


impl<X, R> RModOps<R, Self> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<X, R> RModOps<R, LinComb<X, R>> for &LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<X, R> RMod for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type R = R;
}

#[cfg(test)]
mod tests {
    use std::fmt::Display;
    use num_traits::Zero;
    
    use crate::utils::collections::map;
    use crate::math::traits::AlgBase;
    use super::{FreeGenerator, LinComb};
 
    #[derive(Debug, Default, Hash, PartialEq, Eq, Clone, PartialOrd, Ord)]
    struct X(i32);

    impl AlgBase for X { 
        fn symbol() -> String {
            String::from("X")
        }
    }
    impl Display for X {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "X{}", self.0)
        }
    }
    impl FreeGenerator for X {}

    #[test]
    fn symbol() { 
        type L = LinComb<X, i32>;
        let symbol = L::symbol();
        assert_eq!(symbol, "Free<X; Z>");
    }

    #[test]
    fn fmt() { 
        type L = LinComb<X, i32>;

        let z = L::new(map!{ X(1) => 1 });
        assert_eq!(z.to_string(), "X1");

        let z = L::new(map!{ X(1) => -1 });
        assert_eq!(z.to_string(), "-X1");

        let z = L::new(map!{ X(1) => 2 });
        assert_eq!(z.to_string(), "2X1");

        let z = L::new(map!{ X(1) => 1, X(2) => 1 });
        assert_eq!(z.to_string(), "X1 + X2");

        let z = L::new(map!{ X(1) => -1, X(2) => -1 });
        assert_eq!(z.to_string(), "-X1 - X2");

        let z = L::new(map!{ X(1) => 2, X(2) => 3 });
        assert_eq!(z.to_string(), "2X1 + 3X2");

        let z = L::new(map!{ X(1) => -2, X(2) => -3 });
        assert_eq!(z.to_string(), "-2X1 - 3X2");
    }

    #[test]
    fn default() { 
        type L = LinComb<X, i32>;
        let z = L::default();
        assert!(z.data.is_empty());
    }

    #[test]
    fn eq() { 
        type L = LinComb<X, i32>;
        let z1 = L::new(map!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(map!{ X(2) => 2, X(1) => 1 });
        let z3 = L::new(map!{ X(1) => 1 });

        assert_eq!(z1, z2);
        assert_ne!(z1, z3);
    }

    #[test]
    fn zero() { 
        type L = LinComb<X, i32>;
        let z = L::zero();

        assert!(z.data.is_empty());
        assert!(z.is_zero());
    }

    #[test]
    fn reduce() { 
        type L = LinComb<X, i32>;

        let data = map!{ X(1) => 1, X(2) => 0, X(3) => 0 };
        let mut z = L::new_raw(data);
        
        assert_eq!(z.len(), 3);

        z.reduce();

        assert_eq!(z, L::new_raw(map!{ X(1) => 1 }));
        assert_eq!(z.len(), 1);
    }

    #[test]
    fn add() {
        type L = LinComb<X, i32>;
        let z1 = L::new(map!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(map!{ X(2) => 20, X(3) => 30 });
        let w = z1 + z2;

        assert_eq!(w, L::new(map!{ X(1) => 1, X(2) => 22, X(3) => 30 }));
    }

    #[test]
    fn add_ref() {
        type L = LinComb<X, i32>;
        let z1 = L::new(map!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(map!{ X(2) => 20, X(3) => 30 });
        let w = &z1 + &z2;

        assert_eq!(w, L::new(map!{ X(1) => 1, X(2) => 22, X(3) => 30 }));
    }

    #[test]
    fn add_assign() {
        type L = LinComb<X, i32>;
        let mut z1 = L::new(map!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(map!{ X(2) => 20, X(3) => 30 });
        z1 += z2;

        assert_eq!(z1, L::new(map!{ X(1) => 1, X(2) => 22, X(3) => 30 }));
    }

    #[test]
    fn add_assign_ref() {
        type L = LinComb<X, i32>;
        let mut z1 = L::new(map!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(map!{ X(2) => 20, X(3) => 30 });
        z1 += &z2;

        assert_eq!(z1, L::new(map!{ X(1) => 1, X(2) => 22, X(3) => 30 }));
    }

    #[test]
    fn sum() {
        type L = LinComb<X, i32>;
        let z1 = L::new(map!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(map!{ X(2) => 20, X(3) => 30 });
        let z3 = L::new(map!{ X(3) => 300, X(4) => 400 });
        let w: L = [z1, z2, z3].into_iter().sum();

        assert_eq!(w, L::new(map!{ X(1) => 1, X(2) => 22, X(3) => 330, X(4) => 400 }));
    }

    #[test]
    fn sum_ref() {
        type L = LinComb<X, i32>;
        let z1 = L::new(map!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(map!{ X(2) => 20, X(3) => 30 });
        let z3 = L::new(map!{ X(3) => 300, X(4) => 400 });
        let w: L = [&z1, &z2, &z3].into_iter().sum();

        assert_eq!(w, L::new(map!{ X(1) => 1, X(2) => 22, X(3) => 330, X(4) => 400 }));
    }

    #[test]
    fn neg() {
        type L = LinComb<X, i32>;
        let z = L::new(map!{ X(1) => 1, X(2) => 2 });
        assert_eq!(-z, L::new(map!{ X(1) => -1, X(2) => -2 }));
    }

    #[test]
    fn neg_ref() {
        type L = LinComb<X, i32>;
        let z = L::new(map!{ X(1) => 1, X(2) => 2 });
        assert_eq!(-(&z), L::new(map!{ X(1) => -1, X(2) => -2 }));
    }

    #[test]
    fn sub() {
        type L = LinComb<X, i32>;
        let z1 = L::new(map!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(map!{ X(2) => 20, X(3) => 30 });
        let w = z1 - z2;

        assert_eq!(w, L::new(map!{ X(1) => 1, X(2) => -18, X(3) => -30 }));
    }

    #[test]
    fn sub_ref() {
        type L = LinComb<X, i32>;
        let z1 = L::new(map!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(map!{ X(2) => 20, X(3) => 30 });
        let w = &z1 - &z2;

        assert_eq!(w, L::new(map!{ X(1) => 1, X(2) => -18, X(3) => -30 }));
    }

    #[test]
    fn sub_assign() {
        type L = LinComb<X, i32>;
        let mut z1 = L::new(map!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(map!{ X(2) => 20, X(3) => 30 });
        z1 -= z2;

        assert_eq!(z1, L::new(map!{ X(1) => 1, X(2) => -18, X(3) => -30 }));
    }

    #[test]
    fn sub_assign_ref() {
        type L = LinComb<X, i32>;
        let mut z1 = L::new(map!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(map!{ X(2) => 20, X(3) => 30 });
        z1 -= &z2;

        assert_eq!(z1, L::new(map!{ X(1) => 1, X(2) => -18, X(3) => -30 }));
    }

    #[test]
    fn mul() {
        type L = LinComb<X, i32>;
        let z = L::new(map!{ X(1) => 1, X(2) => 2 });
        let r = 2;
        let w = z * r;

        assert_eq!(w, L::new(map!{ X(1) => 2, X(2) => 4 }));
    }

    #[test]
    fn mul_ref() {
        type L = LinComb<X, i32>;
        let z = L::new(map!{ X(1) => 1, X(2) => 2 });
        let r = 2;
        let w = &z * &r;

        assert_eq!(w, L::new(map!{ X(1) => 2, X(2) => 4 }));
    }

    #[test]
    fn mul_assign() {
        type L = LinComb<X, i32>;
        let mut z = L::new(map!{ X(1) => 1, X(2) => 2 });
        let r = 2;
        z *= r;

        assert_eq!(z, L::new(map!{ X(1) => 2, X(2) => 4 }));
    }

    #[test]
    fn mul_assign_ref() {
        type L = LinComb<X, i32>;
        let mut z = L::new(map!{ X(1) => 1, X(2) => 2 });
        let r = 2;
        z *= &r;

        assert_eq!(z, L::new(map!{ X(1) => 2, X(2) => 4 }));
    }
}