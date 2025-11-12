use std::ops::{Add, AddAssign};
use num_traits::Zero;
use crate::Elem;

// Additive Monoids 
pub trait AddMonOps<T = Self>: 
    Sized + 
    Add<T, Output = T> +              // S + T -> T
    for<'a> Add<&'a T, Output = T>    // S + &T -> T
{}

pub trait AddMon: 
    Elem + 
    Zero +
    AddMonOps +                       // T + T -> T, T + &T -> T
    AddAssign +                       // T += T
    for<'a> AddAssign<&'a Self>       // T += &T
where 
    for<'a> &'a Self: AddMonOps<Self> // &T + T -> T, &T + &T -> T
{
    fn sum<A, I>(itr: I) -> Self 
    where 
        Self: AddAssign<A>,
        I: IntoIterator<Item = A> 
    { 
        itr.into_iter().fold(Self::zero(), |mut res, a| { 
            res += a;
            res
        })
    }
}

#[cfg(test)]
mod tests { 
    use super::*;
    #[test]
    fn sum() { 
        let a = i64::sum([4,5,6]);
        assert_eq!(a, 15);
    }
}