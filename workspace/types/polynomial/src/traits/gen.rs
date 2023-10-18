use std::ops::Mul;
use num_traits::One;
use yui_lin_comb::Gen;

use super::MonoDeg;

pub trait MonoGen: 
    Mul<Output = Self> + 
    One + 
    PartialOrd + 
    Ord + 
    From<Self::Degree> +
    Gen
{
    type Degree: MonoDeg;
    fn degree(&self) -> Self::Degree;
    fn is_unit(&self) -> bool;
    fn inv(&self) -> Option<Self>;
}