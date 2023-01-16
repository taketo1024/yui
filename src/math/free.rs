use std::collections::HashMap;
use std::fmt::{Display, Debug};
use std::hash::Hash;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, MulAssign};
use num_traits::Zero;

use super::traits::{Symbol, AlgBase, AddMon, AddMonOps, AddGrp, AddGrpOps, Ring, RingOps, RMod, RModOps};

pub trait FreeGenerator: Clone + PartialEq + Eq + Hash + Display + Debug + Send + Sync + Symbol {}

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
    fn new(data: HashMap<X, R>) -> Self {
        Self { data }
    }

    pub fn data(&self) -> &HashMap<X, R> {
        &self.data
    }

    pub fn drop_zeros(self) -> Self { 
        let data = self.data.into_iter().filter(|(_, r)| !r.is_zero()).collect();
        Self::new(data)
    }
}

impl<X, R> From<Vec<(X, R)>> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn from(data: Vec<(X, R)>) -> Self {
        let data = data.into_iter().collect::<HashMap<_, _>>();
        Self::new(data)
    }
}

impl<X, R> Default for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn default() -> Self {
        Self::new(HashMap::new())
    }
}

impl<X, R> Display for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        todo!()
    }
}

impl<X, R> Symbol for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn symbol() -> String {
        format!("LinComb<{}; {}>", X::symbol(), R::symbol())
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

    fn add(self, rhs: Self) -> Self::Output {
        let mut res = self;
        res += rhs;
        res
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

impl<X, R> AddAssign for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn add_assign(&mut self, rhs: Self) {
        let data = &mut self.data;
        for (x, r) in rhs.data.into_iter() { 
            if data.contains_key(&x) { 
                let v = data.get_mut(&x).unwrap();
                *v += r;
            } else { 
                data.insert(x, r);
            }
        }
    }
}

impl<'a, X, R> AddAssign<&'a Self> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn add_assign(&mut self, rhs: &'a Self) {
        let data = &mut self.data;
        for (x, r) in rhs.data.iter() { 
            if data.contains_key(&x) { 
                let v = data.get_mut(&x).unwrap();
                *v += r;
            } else { 
                data.insert(x.clone(), r.clone());
            }
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
        todo!()
    }
}

impl<'a, X, R> Neg for &'a LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = LinComb<X, R>;

    fn neg(self) -> Self::Output {
        todo!()
    }
}

impl<X, R> Sub for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        todo!()
    }
}

impl<'a, X, R> Sub for &'a LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = LinComb<X, R>;

    fn sub(self, rhs: Self) -> Self::Output {
        todo!()
    }
}

impl<X, R> SubAssign for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn sub_assign(&mut self, rhs: Self) {
        todo!()
    }
}

impl<'a, X, R> SubAssign<&'a Self> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        todo!()
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
        &self * &rhs
    }
}

impl<'a, X, R> Mul<&'a R> for &'a LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = LinComb<X, R>;

    fn mul(self, rhs: &'a R) -> Self::Output {
        todo!()
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

impl<X, R> MulAssign<R> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn mul_assign(&mut self, rhs: R) {
        self.mul_assign(&rhs)
    }
}

impl<'a, X, R> MulAssign<&'a R> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn mul_assign(&mut self, rhs: &'a R) {
        todo!()
    }
}

impl<X, R> RMod for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type R = R;
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use crate::math::traits::Symbol;
    use crate::utils::collections::hashmap;
    use derive_more::Display;
    use num_traits::Zero;

    use super::{FreeGenerator, LinComb};
 
    #[derive(Debug, Display, Hash, PartialEq, Eq, Clone)]
    struct X(i32);

    impl Symbol for X { 
        fn symbol() -> String {
            String::from("X")
        }
    }
    impl FreeGenerator for X {}

    #[test]
    fn symbol() { 
        type L = LinComb<X, i32>;
        let symbol = L::symbol();
        assert_eq!(symbol, "LinComb<X; Z>");
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

}