use crate::{EucRing, EucRingOps};

// Fields

pub trait FieldOps<T = Self>: 
    EucRingOps<T>
{}

pub trait Field: 
    EucRing + 
    FieldOps
where 
    for<'a> &'a Self: FieldOps<Self> 
{}