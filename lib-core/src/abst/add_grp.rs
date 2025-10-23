use std::ops::{Neg, Sub, SubAssign};
use crate::{AddMon, AddMonOps};

// Additive Groups 

pub trait AddGrpOps<T = Self>: 
    AddMonOps<T> + 
    Neg<Output = T> + 
    Sub<T, Output = T> +
    for<'a> Sub<&'a T, Output = T> 
{}

pub trait AddGrp: 
    AddMon + 
    AddGrpOps + 
    SubAssign + 
    for<'a> SubAssign<&'a Self>
where 
    for<'a> &'a Self: AddGrpOps<Self>
{}