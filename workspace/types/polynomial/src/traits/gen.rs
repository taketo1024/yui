use std::ops::{Mul, Div};
use num_traits::One;
use yui_lin_comb::Gen;

pub trait PolyGen: 
    Mul<Output = Self> + 
    Div<Output = Self> + 
    One + 
    PartialOrd + 
    Ord + 
    From<Self::Deg> +
    Gen
{
    type Deg;

    fn deg(&self) -> Self::Deg;
    fn is_unit(&self) -> bool;
    fn inv(&self) -> Option<Self>;
    fn divides(&self, other: &Self) -> bool;
}