use std::iter::Sum;
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
    AddMonOps +                       // T + T -> T, T + &T -> T
    AddAssign +                       // T += T
    for<'a> AddAssign<&'a Self> +     // T += &T
    Sum<Self> +                       // [T] -> T
    for<'a> Sum<&'a Self> +           // [&T] -> T
    Zero
where 
    for<'a> &'a Self: AddMonOps<Self> // &T + T -> T, &T + &T -> T
{}