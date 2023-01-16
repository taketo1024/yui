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
    pub fn new(data: HashMap<X, R>) -> Self {
        Self { data }
    }

    pub fn data(&self) -> &HashMap<X, R> {
        &self.data
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

impl<X, R> Add for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl<'a, X, R> Add for &'a LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = LinComb<X, R>;

    fn add(self, rhs: Self) -> Self::Output {
        todo!()
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

impl<X, R> AddAssign for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn add_assign(&mut self, rhs: Self) {
        todo!()
    }
}

impl<'a, X, R> AddAssign<&'a Self> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn add_assign(&mut self, rhs: &'a Self) {
        todo!()
    }
}

impl<X, R> Sum for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        todo!()
    }
}

impl<'a, X, R> Sum<&'a Self> for LinComb<X, R>
where
    X: FreeGenerator,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        todo!()
    }
}

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