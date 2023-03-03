use std::iter::Product;
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
    MonOps + 
    MulAssign + 
    for<'a> MulAssign<&'a Self> + 
    Product<Self> + 
    for<'a> Product<&'a Self> +
    One
where
    for<'a> &'a Self: MonOps<Self>
{}