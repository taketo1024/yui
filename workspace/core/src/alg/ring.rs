use crate::{Sign, AddGrp, AddGrpOps, Mon, MonOps};

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
    fn from_sign(e: Sign) -> Self {
        if e.is_positive() { Self::one() } else { -Self::one() }
    }
}

