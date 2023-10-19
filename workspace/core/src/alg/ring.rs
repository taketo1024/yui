use crate::{AddGrp, AddGrpOps, Mon, MonOps, Sign};

// Rings 

pub trait RingOps<T = Self>: 
    AddGrpOps<T> + 
    MonOps<T>
{}

pub trait Ring: 
    AddGrp + 
    Mon + 
    RingOps + 
    From<i32>
where
    for<'a> &'a Self: RingOps<Self>
{
    fn inv(&self) -> Option<Self>;
    fn is_unit(&self) -> bool;
    fn normalizing_unit(&self) -> Self;
    fn from_sign(s: Sign) -> Self { 
        Self::from(s.to_i32())
    }
}