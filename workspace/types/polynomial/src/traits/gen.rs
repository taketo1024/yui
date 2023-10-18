use std::ops::Mul;
use num_traits::One;
use yui_lin_comb::Gen;

use super::MonoDeg;

pub trait MonoGen: 
    Mul<Output = Self> + 
    One + 
    PartialOrd + 
    Ord + 
    From<Self::Deg> +
    Gen
{
    type Deg: MonoDeg;

    fn deg(&self) -> Self::Deg;
    fn is_unit(&self) -> bool;
    fn inv(&self) -> Option<Self>;
}