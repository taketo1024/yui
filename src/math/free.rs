use std::collections::HashMap;
use std::fmt::{Display, Debug};
use std::hash::Hash;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, MulAssign};
use itertools::Itertools;
use num_traits::Zero;

use crate::utils::collections::hashmap;
use crate::utils::display::OrdForDisplay;

use super::traits::{Symbol, AlgBase, AddMon, AddMonOps, AddGrp, AddGrpOps, Ring, RingOps, RMod, RModOps};

pub trait FreeGenerator: Clone + PartialEq + Eq + Hash + Display + Debug + Send + Sync + Symbol + OrdForDisplay {}

#[derive(PartialEq, Eq, Clone, Debug)]
pub struct LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{ 
    data: HashMap<X, R>
}

impl<X, R> LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{ 
    pub fn new(data: HashMap<X, R>) -> Self {
        Self { data }
    }

    pub fn wrap(x: X) -> Self { 
        Self::new(hashmap!{ x => R::one() })
    }

    pub fn unwrap(self) -> X { 
        if self.data.len() == 1 {
            if let Some((x, r)) = self.data.into_iter().next() {
                if r.is_one() { 
                    return x
                }
            }
        } 
        panic!()
    }

    pub fn from_iter<I>(iter: I) -> Self 
    where I: Iterator<Item = (X, R)> {
        let data = iter.collect::<HashMap<_, _>>();
        Self::new(data)
    }

    pub fn iter(&self) -> std::collections::hash_map::Iter<'_, X, R> {
        self.data.iter()
    }

    pub fn into_iter(self) -> std::collections::hash_map::IntoIter<X, R> {
        self.data.into_iter()
    }

    pub fn drop_zeros(self) -> Self { 
        let data = self.into_iter().filter(|(_, r)| !r.is_zero()).collect();
        Self::new(data)
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

    pub fn len(&self) -> usize {
        self.data.len()
    }
}

impl<X, R> From<(X, R)> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn from(pair: (X, R)) -> Self {
        Self::new(hashmap!{ pair.0 => pair.1 })
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

impl<X, R> Default for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn default() -> Self {
        Self::new(HashMap::default())
    }
}

impl<X, R> Display for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut initial = true;
        let sorted = self.iter().sorted_by(|(x, _), (y, _)| x.cmp_for_display(y) );
        for (x, r) in sorted {
            let x = x.to_string();
            let r = r.to_string();

            if initial { 
                let r = if r == "1" { 
                    ""
                } else if r == "-1" { 
                    "-"
                } else {
                    &r
                };
                write!(f, "{r}{x}")?;
                initial = false
            } else {
                let sign = if r.starts_with("-") { "-" } else { "+" };
                let r = if r == "1" || r == "-1" { 
                    ""
                } else if r.starts_with("-") { 
                    &r[1..]
                } else { 
                    &r
                };
                write!(f, " {sign} {r}{x}")?;
            }
        }
        Ok(())
    }
}

impl<X, R> Symbol for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn symbol() -> String {
        format!("Free<{}; {}>", X::symbol(), R::symbol())
    }
}

impl<X, R> AlgBase for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<X, R> Zero for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn zero() -> Self {
        Self::default()
    }

    fn is_zero(&self) -> bool {
        self.data.values().all(|r| r.is_zero())
    }
}

impl<X, R> Add for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<'a, X, R> Add for &'a LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = LinComb<X, R>;

    fn add(self, rhs: Self) -> Self::Output {
        let mut res = self.clone();
        res += rhs;
        res
    }
}

impl<X, R> AddAssign<(X, R)> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn add_assign(&mut self, rhs: (X, R)) {
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

impl<'a, X, R> AddAssign<(&'a X, &'a R)> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn add_assign(&mut self, rhs: (&'a X, &'a R)) {
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

impl<X, R> AddAssign for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn add_assign(&mut self, rhs: Self) {
        for e in rhs.data.into_iter() { 
            self.add_assign(e);
        }
    }
}

impl<'a, X, R> AddAssign<&'a Self> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn add_assign(&mut self, rhs: &'a Self) {
        for e in rhs.data.iter() { 
            self.add_assign(e);
        }
    }
}

impl<X, R> Sum for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut res = Self::zero();
        for z in iter { 
            res += z
        }
        res
    }
}

impl<'a, X, R> Sum<&'a Self> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        let mut res = Self::zero();
        for z in iter { 
            res += z
        }
        res
    }
}

impl<X, R> AddMonOps<Self> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<'a, X, R> AddMonOps<LinComb<X, R>> for &'a LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<X, R> AddMon for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

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

impl<'a, X, R> Neg for &'a LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = LinComb<X, R>;

    fn neg(self) -> Self::Output {
        self.map_coeffs(|r| -r)
    }
}

impl<X, R> Sub for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = Self;

    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<'a, X, R> Sub for &'a LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = LinComb<X, R>;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut res = self.clone();
        res -= rhs;
        res
    }
}

impl<X, R> SubAssign<(X, R)> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn sub_assign(&mut self, rhs: (X, R)) {
        let data = &mut self.data;
        let (x, r) = rhs;

        if data.contains_key(&x) { 
            let v = data.get_mut(&x).unwrap();
            *v -= r;
        } else { 
            data.insert(x, -r);
        }
    }
}

impl<'a, X, R> SubAssign<(&'a X, &'a R)> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn sub_assign(&mut self, rhs: (&'a X, &'a R)) {
        let data = &mut self.data;
        let (x, r) = rhs;

        if data.contains_key(x) { 
            let v = data.get_mut(x).unwrap();
            v.sub_assign(r);
        } else { 
            data.insert(x.clone(), -r);
        }
    }
}

impl<X, R> SubAssign for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn sub_assign(&mut self, rhs: Self) {
        for e in rhs.data.into_iter() { 
            self.sub_assign(e);
        }
    }
}

impl<'a, X, R> SubAssign<&'a Self> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        for e in rhs.data.iter() { 
            self.sub_assign(e);
        }
    }
}

impl<X, R> AddGrpOps<Self> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<'a, X, R> AddGrpOps<LinComb<X, R>> for &'a LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<X, R> AddGrp for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<X, R> Mul<R> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = Self;

    fn mul(self, rhs: R) -> Self::Output {
        self.map_coeffs_into(|r| &r * &rhs)
    }
}

impl<'a, X, R> Mul<&'a R> for &'a LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = LinComb<X, R>;

    fn mul(self, rhs: &'a R) -> Self::Output {
        self.map_coeffs(|r| r * rhs)
    }
}

impl<X, R> MulAssign<R> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn mul_assign(&mut self, rhs: R) {
        self.mul_assign(&rhs);
    }
}

impl<'a, X, R> MulAssign<&'a R> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn mul_assign(&mut self, rhs: &'a R) {
        let data = std::mem::take(&mut self.data);
        self.data = data.into_iter().map(|(x, r)| (x, &r * rhs)).collect();
    }
}

impl<X, R> RModOps<R, R, Self> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<'a, X, R> RModOps<R, &'a R, LinComb<X, R>> for &'a LinComb<X, R>
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

impl<X, R> Mul for LinComb<X, R>
where 
    X: FreeGenerator, for<'x> &'x X: Mul<Output = X>,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = Self;
    
    fn mul(mut self, rhs: Self) -> Self {
        self *= rhs;
        self
    }
}

impl<'a, X, R> Mul for &'a LinComb<X, R>
where 
    X: FreeGenerator, for<'x> &'x X: Mul<Output = X>,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = LinComb<X, R>;
    
    fn mul(self, rhs: Self) -> LinComb<X, R> {
        let mut res = self.clone();
        res *= rhs;
        res
    }
}

impl<X, R> MulAssign for LinComb<X, R>
where 
    X: FreeGenerator, for<'x> &'x X: Mul<Output = X>,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn mul_assign(&mut self, rhs: Self) {
        self.mul_assign(&rhs)
    }
}

impl<'a, X, R> MulAssign<&'a LinComb<X, R>> for LinComb<X, R>
where 
    X: FreeGenerator, for<'x> &'x X: Mul<Output = X>,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn mul_assign(&mut self, rhs: &'a Self) {
        let mut data = HashMap::new();
        data.reserve(self.len() * rhs.len());

        for (x, r) in self.iter() { 
            for (y, s) in rhs.iter() { 
                data.insert(x * y, r * s);
            }
        }
        
        self.data = data;
    }
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, fmt::Display};

    use crate::math::traits::Symbol;
    use crate::utils::collections::hashmap;
    use num_traits::Zero;

    use super::{FreeGenerator, LinComb};
 
    #[derive(Debug, Hash, PartialEq, Eq, Clone, PartialOrd, Ord)]
    struct X(i32);

    impl Symbol for X { 
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

        let z = L::new(hashmap!{ X(1) => 1 });
        assert_eq!(z.to_string(), "X1");

        let z = L::new(hashmap!{ X(1) => -1 });
        assert_eq!(z.to_string(), "-X1");

        let z = L::new(hashmap!{ X(1) => 2 });
        assert_eq!(z.to_string(), "2X1");

        let z = L::new(hashmap!{ X(1) => 1, X(2) => 1 });
        assert_eq!(z.to_string(), "X1 + X2");

        let z = L::new(hashmap!{ X(1) => -1, X(2) => -1 });
        assert_eq!(z.to_string(), "-X1 - X2");

        let z = L::new(hashmap!{ X(1) => 2, X(2) => 3 });
        assert_eq!(z.to_string(), "2X1 + 3X2");

        let z = L::new(hashmap!{ X(1) => -2, X(2) => -3 });
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
        let z1 = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(hashmap!{ X(2) => 2, X(1) => 1 });
        let z3 = L::new(hashmap!{ X(1) => 1 });

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
    fn drop_zeros() { 
        type L = LinComb<X, i32>;
        let z = L::new(hashmap!{ X(1) => 1, X(2) => 0, X(3) => 0 });
        let z = z.drop_zeros();

        assert_eq!(z, L::new(hashmap!{ X(1) => 1 }));
    }

    #[test]
    fn add() {
        type L = LinComb<X, i32>;
        let z1 = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(hashmap!{ X(2) => 20, X(3) => 30 });
        let w = z1 + z2;

        assert_eq!(w, L::new(hashmap!{ X(1) => 1, X(2) => 22, X(3) => 30 }));
    }

    #[test]
    fn add_ref() {
        type L = LinComb<X, i32>;
        let z1 = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(hashmap!{ X(2) => 20, X(3) => 30 });
        let w = &z1 + &z2;

        assert_eq!(w, L::new(hashmap!{ X(1) => 1, X(2) => 22, X(3) => 30 }));
    }

    #[test]
    fn add_assign() {
        type L = LinComb<X, i32>;
        let mut z1 = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(hashmap!{ X(2) => 20, X(3) => 30 });
        z1 += z2;

        assert_eq!(z1, L::new(hashmap!{ X(1) => 1, X(2) => 22, X(3) => 30 }));
    }

    #[test]
    fn add_assign_ref() {
        type L = LinComb<X, i32>;
        let mut z1 = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(hashmap!{ X(2) => 20, X(3) => 30 });
        z1 += &z2;

        assert_eq!(z1, L::new(hashmap!{ X(1) => 1, X(2) => 22, X(3) => 30 }));
    }

    #[test]
    fn sum() {
        type L = LinComb<X, i32>;
        let z1 = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(hashmap!{ X(2) => 20, X(3) => 30 });
        let z3 = L::new(hashmap!{ X(3) => 300, X(4) => 400 });
        let w: L = [z1, z2, z3].into_iter().sum();

        assert_eq!(w, L::new(hashmap!{ X(1) => 1, X(2) => 22, X(3) => 330, X(4) => 400 }));
    }

    #[test]
    fn sum_ref() {
        type L = LinComb<X, i32>;
        let z1 = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(hashmap!{ X(2) => 20, X(3) => 30 });
        let z3 = L::new(hashmap!{ X(3) => 300, X(4) => 400 });
        let w: L = [&z1, &z2, &z3].into_iter().sum();

        assert_eq!(w, L::new(hashmap!{ X(1) => 1, X(2) => 22, X(3) => 330, X(4) => 400 }));
    }

    #[test]
    fn neg() {
        type L = LinComb<X, i32>;
        let z = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        assert_eq!(-z, L::new(hashmap!{ X(1) => -1, X(2) => -2 }));
    }

    #[test]
    fn neg_ref() {
        type L = LinComb<X, i32>;
        let z = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        assert_eq!(-(&z), L::new(hashmap!{ X(1) => -1, X(2) => -2 }));
    }

    #[test]
    fn sub() {
        type L = LinComb<X, i32>;
        let z1 = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(hashmap!{ X(2) => 20, X(3) => 30 });
        let w = z1 - z2;

        assert_eq!(w, L::new(hashmap!{ X(1) => 1, X(2) => -18, X(3) => -30 }));
    }

    #[test]
    fn sub_ref() {
        type L = LinComb<X, i32>;
        let z1 = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(hashmap!{ X(2) => 20, X(3) => 30 });
        let w = &z1 - &z2;

        assert_eq!(w, L::new(hashmap!{ X(1) => 1, X(2) => -18, X(3) => -30 }));
    }

    #[test]
    fn sub_assign() {
        type L = LinComb<X, i32>;
        let mut z1 = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(hashmap!{ X(2) => 20, X(3) => 30 });
        z1 -= z2;

        assert_eq!(z1, L::new(hashmap!{ X(1) => 1, X(2) => -18, X(3) => -30 }));
    }

    #[test]
    fn sub_assign_ref() {
        type L = LinComb<X, i32>;
        let mut z1 = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let z2 = L::new(hashmap!{ X(2) => 20, X(3) => 30 });
        z1 -= &z2;

        assert_eq!(z1, L::new(hashmap!{ X(1) => 1, X(2) => -18, X(3) => -30 }));
    }

    #[test]
    fn mul() {
        type L = LinComb<X, i32>;
        let z = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let r = 2;
        let w = z * r;

        assert_eq!(w, L::new(hashmap!{ X(1) => 2, X(2) => 4 }));
    }

    #[test]
    fn mul_ref() {
        type L = LinComb<X, i32>;
        let z = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let r = 2;
        let w = &z * &r;

        assert_eq!(w, L::new(hashmap!{ X(1) => 2, X(2) => 4 }));
    }

    #[test]
    fn mul_assign() {
        type L = LinComb<X, i32>;
        let mut z = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let r = 2;
        z *= r;

        assert_eq!(z, L::new(hashmap!{ X(1) => 2, X(2) => 4 }));
    }

    #[test]
    fn mul_assign_ref() {
        type L = LinComb<X, i32>;
        let mut z = L::new(hashmap!{ X(1) => 1, X(2) => 2 });
        let r = 2;
        z *= &r;

        assert_eq!(z, L::new(hashmap!{ X(1) => 2, X(2) => 4 }));
    }
}