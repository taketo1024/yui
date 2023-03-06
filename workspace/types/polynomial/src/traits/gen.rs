use std::ops::Mul;
use num_traits::One;
use yui_lin_comb::FreeGen;
use super::PolyDeg;

pub trait PolyGen: 
    Mul<Output = Self> + 
    One + 
    PartialOrd + 
    Ord + 
    From<Self::Degree> +
    FreeGen
{
    type Degree: PolyDeg;
    fn degree(&self) -> Self::Degree;
    fn is_unit(&self) -> bool;
    fn inv(&self) -> Option<Self>;
}