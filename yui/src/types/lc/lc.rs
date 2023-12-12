use std::collections::HashMap;
use std::fmt::{Display, Debug};
use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, MulAssign};
use ahash::AHashMap;
use itertools::Itertools;
use num_traits::Zero;
use auto_impl_ops::auto_ops;
use crate::{Elem, AddMon, AddMonOps, AddGrp, AddGrpOps, Ring, RingOps, RMod, RModOps};

use super::gen::*;

#[derive(PartialEq, Eq, Clone, Default)]
#[derive(serde::Deserialize, serde::Serialize)]
#[serde(transparent)]
pub struct Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{ 
    data: AHashMap<X, R>,
    #[serde(skip)]
    r_zero: R
}

impl<X, R> Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{ 
    pub fn new() -> Self {
        let hasher = ahash::RandomState::with_seeds(0, 0, 0, 0);
        let data = AHashMap::with_hasher(hasher);
        let r_zero = R::zero();
        Self { data, r_zero }
    }

    pub fn clean(&mut self) { 
        self.data.retain(|_, r| !r.is_zero());
    }

    pub fn nterms(&self) -> usize {
        self.data.len()
    }

    pub fn any_term(&self) -> Option<(&X, &R)> { 
        self.iter().next()
    }

    pub fn gens(&self) -> impl Iterator<Item = &X> {
        self.data.keys()
    }

    pub fn is_gen(&self) -> bool { 
        self.nterms() == 1 && 
        self.iter().next().unwrap().1.is_one()
    }

    pub fn as_gen(&self) -> Option<X> { 
        if !self.is_gen() { 
            None?
        }
        self.iter().next().map(|(x, _)| x.clone())
    }

    pub fn coeff(&self, x: &X) -> &R { 
        self.data.get(x).unwrap_or(&self.r_zero)
    }

    pub fn iter(&self) -> impl Iterator<Item = (&X, &R)> {
        self.data.iter()
    }

    pub fn map<Y, S, F>(&self, f: F) -> Lc<Y, S>
    where 
        Y: Gen, 
        S: Ring, for<'x> &'x S: RingOps<S>,
        F: Fn(&X, &R) -> (Y, S) 
    { 
        self.iter().map(|(x, r)| f(x, r)).collect()
    }

    pub fn into_map<Y, S, F>(self, f: F) -> Lc<Y, S>
    where 
        Y: Gen, 
        S: Ring, for<'x> &'x S: RingOps<S>,
        F: Fn(X, R) -> (Y, S) 
    { 
        self.into_iter().map(|(x, r)| f(x, r)).collect()
    }

    pub fn map_coeffs<S, F>(&self, f: F) -> Lc<X, S>
    where 
        S: Ring, for<'x> &'x S: RingOps<S>, 
        F: Fn(&R) -> S 
    { 
        self.map(|x, r| (x.clone(), f(r)))
    }

    pub fn into_map_coeffs<S, F>(self, f: F) -> Lc<X, S>
    where 
        S: Ring, for<'x> &'x S: RingOps<S>, 
        F: Fn(R) -> S 
    { 
        self.into_map(|x, r| (x, f(r)))
    }

    pub fn map_gens<Y, F>(&self, f: F) -> Lc<Y, R>
    where 
        Y: Gen, 
        F: Fn(&X) -> Y 
    { 
        self.map(|x, r| (f(x), r.clone()))
    }

    pub fn into_map_gens<Y, F>(self, f: F) -> Lc<Y, R>
    where 
        Y: Gen, 
        F: Fn(X) -> Y 
    { 
        self.into_map(|x, r| (f(x), r))
    }

    pub fn filter_gens<F>(&self, f: F) -> Self
    where F: Fn(&X) -> bool { 
        self.iter().filter_map(|(x, a)| 
            if f(x) { 
                Some((x.clone(), a.clone()))
            } else { 
                None
            }
        ).collect()
    }

    pub fn into_filter_gens<F>(self, f: F) -> Self
    where F: Fn(&X) -> bool { 
        self.into_iter().filter(|(x, _)| f(&x)).collect()
    }

    pub fn apply<F>(&self, f: F) -> Self 
    where F: Fn(&X) -> Lc<X, R> {
        self.iter().flat_map(|(x, r)| { 
            f(x).into_iter().map(move |(y, s)| { 
                (y, r * &s)
            })
        }).collect()
    }

    pub fn fmt(&self, f: &mut std::fmt::Formatter<'_>, ascending: bool) -> std::fmt::Result {
        use crate::util::format::paren_expr;

        if self.data.is_empty() { 
            return write!(f, "0");
        }

        let mut elements = self.iter().sorted_by(|(x, _), (y, _)| 
            if ascending { 
                x.cmp_for_display(y) 
            } else {
                x.cmp_for_display(y).reverse()
            }
        );
        
        if let Some((x, r)) = elements.next() {
            let r = paren_expr(r);
            let x = x.to_string();

            if r == "1" { 
                write!(f, "{x}")?
            } else if r == "-1" { 
                write!(f, "-{x}")?
            } else if x == "1" {
                write!(f, "{r}")?
            } else { 
                write!(f, "{r}{x}")?
            };
        };

        for (x, r) in elements {
            let r = paren_expr(r);
            let x = x.to_string();

            let (op, r) = if let Some(r) = r.strip_prefix('-') { 
                ("-", r) 
            } else { 
                ("+", r.as_str())
            };

            if r == "1" { 
                write!(f, " {op} {x}")?
            } else if x == "1" { 
                write!(f, " {op} {r}")?
            } else { 
                write!(f, " {op} {r}{x}")?
            };
        }

        Ok(())
    }
}

impl<X, R> From<X> for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn from(x: X) -> Self {
        Self::from((x, R::one()))
    }    
}

impl<X, R> From<(X, R)> for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn from(value: (X, R)) -> Self {
        Self::from_iter([value])
    }
}

impl<X, R> From<HashMap<X, R>> for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn from(value: HashMap<X, R>) -> Self {
        Self::from_iter(value.into_iter())
    }
}

impl<X, R> FromIterator<(X, R)> for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn from_iter<T: IntoIterator<Item = (X, R)>>(iter: T) -> Self {
        let mut res = Self::new();
        for e in iter.into_iter() { 
            res.add_pair(e);
        }
        res.clean();
        res
    }
}

impl<X, R> IntoIterator for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Item = (X, R);
    type IntoIter = std::collections::hash_map::IntoIter<X, R>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

impl<X, R> Display for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt(f, true)
    }
}

impl<X, R> Debug for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt(f, true)
    }
}

impl<X, R> Zero for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn zero() -> Self {
        Self::new()
    }

    fn is_zero(&self) -> bool {
        self.data.is_empty()
    }
}

impl<X, R> Neg for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.into_map_coeffs(|r| -r)
    }
}

impl<X, R> Neg for &Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = Lc<X, R>;

    fn neg(self) -> Self::Output {
        self.map_coeffs(|r| -r)
    }
}

impl<X, R> Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    // must clean after call
    pub fn add_pair(&mut self, rhs: (X, R)) { 
        let (x, r) = rhs;
        if r.is_zero() { return }

        if self.data.contains_key(&x) { 
            let v = self.data.get_mut(&x).unwrap();
            v.add_assign(r);
        } else { 
            self.data.insert(x, r);
        }
    } 

    // must clean after call
    pub fn add_pair_ref(&mut self, rhs: (&X, &R)) { 
        let (x, r) = rhs;
        if r.is_zero() { return }

        if self.data.contains_key(x) { 
            let v = self.data.get_mut(x).unwrap();
            v.add_assign(r);
        } else { 
            self.data.insert(x.clone(), r.clone());
        }
    }
}

#[auto_ops]
impl<X, R> AddAssign<&Lc<X, R>> for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn add_assign(&mut self, rhs: &Self) {
        for e in rhs.data.iter() { 
            self.add_pair_ref(e);
        }
        self.clean()
    }
}

#[auto_ops]
impl<X, R> SubAssign<&Lc<X, R>> for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn sub_assign(&mut self, rhs: &Self) {
        for e in rhs.data.iter() { 
            self.add_pair_ref((e.0, &-e.1));
        }
        self.clean()
    }
}

#[auto_ops]
impl<X, R> MulAssign<&R> for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn mul_assign(&mut self, rhs: &R) {
        let data = std::mem::take(&mut self.data);
        self.data = data.into_iter().map(|(x, r)| (x, &r * rhs)).collect();
        self.clean()
    }
}

#[auto_ops]
impl<X, R> Mul for &Lc<X, R>
where 
    X: Gen + Mul<Output = X>,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = Lc<X, R>;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut res = Self::Output::zero();
        res.data.reserve(self.nterms() * rhs.nterms());

        for (x, r) in self.iter() { 
            for (y, s) in rhs.iter() { 
                let xy = x.clone() * y.clone();
                let rs = r * s;
                res.add_pair((xy, rs));
            }
        }
        
        res.clean();
        res
    }
}

macro_rules! impl_alg_ops {
    ($trait:ident) => {
        impl<X, R> $trait<Self> for Lc<X, R>
        where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {}

        impl<X, R> $trait<Lc<X, R>> for &Lc<X, R>
        where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {}
    };
}

impl_alg_ops!(AddMonOps);
impl_alg_ops!(AddGrpOps);

impl<X, R> Elem for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn math_symbol() -> String {
        format!("{}<{}>", R::math_symbol(), X::math_symbol())
    }
}

impl<X, R> AddMon for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<X, R> AddGrp for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{}


impl<X, R> RModOps<R, Self> for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<X, R> RModOps<R, Lc<X, R>> for &Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<X, R> RMod for Lc<X, R>
where
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type R = R;
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use crate::{Elem, AddMon};
    use crate::util::macros::hashmap;
    use crate::lc::{Free, Lc};
 
    type X = Free<i32>;
    fn e(i: i32) -> X { 
        X::from(i)
    }

    #[test]
    fn math_symbol() { 
        type L = Lc<X, i32>;
        let symbol = L::math_symbol();
        assert_eq!(symbol, "Z<Z>");
    }

    #[test]
    fn fmt() { 
        type L = Lc<X, i32>;

        let z = L::from(hashmap!{ e(1) => 1 });
        assert_eq!(z.to_string(), "<1>");

        let z = L::from(hashmap!{ e(1) => -1 });
        assert_eq!(z.to_string(), "-<1>");

        let z = L::from(hashmap!{ e(1) => 2 });
        assert_eq!(z.to_string(), "2<1>");

        let z = L::from(hashmap!{ e(1) => 1, e(2) => 1 });
        assert_eq!(z.to_string(), "<1> + <2>");

        let z = L::from(hashmap!{ e(1) => -1, e(2) => -1 });
        assert_eq!(z.to_string(), "-<1> - <2>");

        let z = L::from(hashmap!{ e(1) => 2, e(2) => 3 });
        assert_eq!(z.to_string(), "2<1> + 3<2>");

        let z = L::from(hashmap!{ e(1) => -2, e(2) => -3 });
        assert_eq!(z.to_string(), "-2<1> - 3<2>");
    }

    #[test]
    fn default() { 
        type L = Lc<X, i32>;
        let z = L::default();
        assert!(z.data.is_empty());
    }

    #[test]
    fn from_gen() { 
        type L = Lc<X, i32>;
        let x = e(0);
        let z = L::from(x);
        assert_eq!(z, L::from(hashmap!{ e(0) => 1 }));
    }

    #[test]
    fn from_pair() { 
        type L = Lc<X, i32>;
        let x = e(0);
        let z = L::from((x, 2));
        assert_eq!(z, L::from(hashmap!{ e(0) => 2 }));
    }

    #[test]
    fn from_iter() { 
        type L = Lc<X, i32>;
        let z = L::from_iter([(e(0), 1), (e(1), 0), (e(2), 2)]);

        assert!(!z.is_zero());
        assert_eq!(z.nterms(), 2);
        assert_eq!(z.coeff(&e(0)), &1);
        assert_eq!(z.coeff(&e(2)), &2);
    }

    #[test]
    fn into_gen() { 
        type L = Lc<X, i32>;
        let z = L::from(e(0));

        assert!(z.is_gen());
        assert_eq!(z.as_gen(), Some(e(0)));

        let z = L::from((e(0), 2));
        assert!(!z.is_gen());
        assert_eq!(z.as_gen(), None);

        let z = L::from_iter([(e(0), 1), (e(1), 1)]);
        assert!(!z.is_gen());
        assert_eq!(z.as_gen(), None);
    }

    #[test]
    fn eq() { 
        type L = Lc<X, i32>;
        let z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 2, e(1) => 1 });
        let z3 = L::from(hashmap!{ e(1) => 1 });

        assert_eq!(z1, z2);
        assert_ne!(z1, z3);
    }

    #[test]
    fn zero() { 
        type L = Lc<X, i32>;
        let z = L::zero();

        assert!(z.data.is_empty());
        assert!(z.is_zero());

        let z = L::from(hashmap!{ e(1) => 1 });

        assert!(!z.data.is_empty());
        assert!(!z.is_zero());
    }

    #[test]
    fn clean() { 
        type L = Lc<X, i32>;

        let mut z = L::from(hashmap!{ e(1) => 1, e(2) => 2, e(3) => 1 });
        z.add_pair((e(1), -1));
        z.add_pair((e(2), -1));
        z.add_pair((e(3), -1));
        
        assert_eq!(z.nterms(), 3);

        z.clean();

        assert_eq!(z, L::from(hashmap!{ e(2) => 1 }));
        assert_eq!(z.nterms(), 1);
    }

    #[test]
    fn add() {
        type L = Lc<X, i32>;
        let z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        let w = z1 + z2;

        assert_eq!(w, L::from(hashmap!{ e(1) => 1, e(2) => 22, e(3) => 30 }));
    }

    #[test]
    fn add_ref() {
        type L = Lc<X, i32>;
        let z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        let w = &z1 + &z2;

        assert_eq!(w, L::from(hashmap!{ e(1) => 1, e(2) => 22, e(3) => 30 }));
    }

    #[test]
    fn add_assign() {
        type L = Lc<X, i32>;
        let mut z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        z1 += z2;

        assert_eq!(z1, L::from(hashmap!{ e(1) => 1, e(2) => 22, e(3) => 30 }));
    }

    #[test]
    fn add_assign_ref() {
        type L = Lc<X, i32>;
        let mut z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        z1 += &z2;

        assert_eq!(z1, L::from(hashmap!{ e(1) => 1, e(2) => 22, e(3) => 30 }));
    }

    #[test]
    fn sum() {
        type L = Lc<X, i32>;
        let z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        let z3 = L::from(hashmap!{ e(3) => 300, e(4) => 400 });
        let w  = L::sum([z1, z2, z3]);

        assert_eq!(w, L::from(hashmap!{ e(1) => 1, e(2) => 22, e(3) => 330, e(4) => 400 }));
    }

    #[test]
    fn sum_ref() {
        type L = Lc<X, i32>;
        let z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        let z3 = L::from(hashmap!{ e(3) => 300, e(4) => 400 });
        let w  = L::sum([&z1, &z2, &z3]);

        assert_eq!(w, L::from(hashmap!{ e(1) => 1, e(2) => 22, e(3) => 330, e(4) => 400 }));
    }

    #[test]
    fn neg() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        assert_eq!(-z, L::from(hashmap!{ e(1) => -1, e(2) => -2 }));
    }

    #[test]
    fn neg_ref() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        assert_eq!(-(&z), L::from(hashmap!{ e(1) => -1, e(2) => -2 }));
    }

    #[test]
    fn sub() {
        type L = Lc<X, i32>;
        let z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        let w = z1 - z2;

        assert_eq!(w, L::from(hashmap!{ e(1) => 1, e(2) => -18, e(3) => -30 }));
    }

    #[test]
    fn sub_ref() {
        type L = Lc<X, i32>;
        let z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        let w = &z1 - &z2;

        assert_eq!(w, L::from(hashmap!{ e(1) => 1, e(2) => -18, e(3) => -30 }));
    }

    #[test]
    fn sub_assign() {
        type L = Lc<X, i32>;
        let mut z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        z1 -= z2;

        assert_eq!(z1, L::from(hashmap!{ e(1) => 1, e(2) => -18, e(3) => -30 }));
    }

    #[test]
    fn sub_assign_ref() {
        type L = Lc<X, i32>;
        let mut z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        z1 -= &z2;

        assert_eq!(z1, L::from(hashmap!{ e(1) => 1, e(2) => -18, e(3) => -30 }));
    }

    #[test]
    fn mul() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let r = 2;
        let w = z * r;

        assert_eq!(w, L::from(hashmap!{ e(1) => 2, e(2) => 4 }));
    }

    #[test]
    fn mul_ref() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let r = 2;
        let w = z * r;

        assert_eq!(w, L::from(hashmap!{ e(1) => 2, e(2) => 4 }));
    }

    #[test]
    fn mul_assign() {
        type L = Lc<X, i32>;
        let mut z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let r = 2;
        z *= r;

        assert_eq!(z, L::from(hashmap!{ e(1) => 2, e(2) => 4 }));
    }

    #[test]
    fn mul_assign_ref() {
        type L = Lc<X, i32>;
        let mut z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let r = 2;
        z *= &r;

        assert_eq!(z, L::from(hashmap!{ e(1) => 2, e(2) => 4 }));
    }

    #[test]
    fn map_coeffs() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let w = z.map_coeffs(|a| a * 10);

        assert_eq!(w, L::from(hashmap!{ e(1) => 10, e(2) => 20 }));
    }

    #[test]
    fn into_map_coeffs() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let w = z.into_map_coeffs(|a| a * 10);

        assert_eq!(w, L::from(hashmap!{ e(1) => 10, e(2) => 20 }));
    }

    #[test]
    fn map_gens() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let w = z.map_gens(|x| e(x.0 * 10));

        assert_eq!(w, L::from(hashmap!{ e(10) => 1, e(20) => 2 }));
    }

    #[test]
    fn into_map_gens() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let w = z.into_map_gens(|x| e(x.0 * 10));

        assert_eq!(w, L::from(hashmap!{ e(10) => 1, e(20) => 2 }));
    }

    #[test]
    fn filter_gens() { 
        type L = Lc<X, i32>;
        let z = L::from_iter( (1..10).map(|i| (e(i), i * 10)) );
        let w = z.filter_gens(|x| x.0 % 3 == 0 );
        assert_eq!(w, L::from(hashmap!{ e(3) => 30, e(6) => 60, e(9) => 90}))
    }

    #[test]
    fn serialize() { 
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let ser = serde_json::to_string(&z).unwrap();
        let deser = serde_json::from_str::<L>(&ser).unwrap();
        assert_eq!(z, deser);
    }
}