use std::ops::{Mul, MulAssign};
use num_traits::One;
use crate::Elem;

// Monoids (multiplicative)

pub trait MonOps<T = Self>: 
    Sized + 
    Mul<T, Output = T> + 
    for<'a> Mul<&'a T, Output = T> 
{}

pub trait Mon: 
    Elem + 
    One +
    MonOps + 
    MulAssign + 
    for<'a> MulAssign<&'a Self>
where
    for<'a> &'a Self: MonOps<Self>
{
    fn product<A, I>(itr: I) -> Self 
    where 
        Self: MulAssign<A>,
        I: IntoIterator<Item = A> 
    { 
        itr.into_iter().fold(Self::one(), |mut res, a| { 
            res *= a;
            res
        })
    }
}

#[cfg(test)]
mod tests { 
    use super::*;
    #[test]
    fn product() { 
        let a = i64::product([4,5,6]);
        assert_eq!(a, 120);
    }
}