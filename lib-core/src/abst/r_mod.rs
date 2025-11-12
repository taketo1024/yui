use std::ops::{Mul, MulAssign};
use crate::{AddGrp, AddGrpOps, Ring, RingOps};

// R-Modules

pub trait RModOps<R, T>: 
    AddGrpOps<T> + 
    Mul<R, Output = T> + 
    for<'a> Mul<&'a R, Output = T>
where 
    R: Ring, for<'x> &'x R: RingOps<R>
{}

pub trait RMod:
    AddGrp + 
    RModOps<Self::R, Self> + 
    MulAssign<Self::R> +
    for<'a> MulAssign<&'a Self::R>
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>, 
    for<'a> &'a Self: RModOps<Self::R, Self>,
{
    type R;
}

