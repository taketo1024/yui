use std::ops::{Mul, Div};
use num_traits::One;
use crate::lc::Gen;

pub trait MonoOrd { 
    fn cmp_lex(&self, other: &Self) -> std::cmp::Ordering;
    fn cmp_grlex(&self, other: &Self) -> std::cmp::Ordering;
}

pub trait Mono: 
    From<Self::Deg> +
    One + 
    Mul<Output = Self> + 
    Div<Output = Self> + 
    MonoOrd + 
    Gen
{
    type Deg;

    fn deg(&self) -> Self::Deg;
    fn is_unit(&self) -> bool;
    fn inv(&self) -> Option<Self>;
    fn divides(&self, other: &Self) -> bool;
}