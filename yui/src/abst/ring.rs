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

    fn is_pm_one(&self) -> bool { 
        self.is_one() || (-self).is_one()
    }
    fn from_sign(s: Sign) -> Self { 
        Self::from(s.to_i32())
    }
}

#[cfg(test)]
mod tests {
    use crate::Ring;
 
    #[test]
    fn is_pm_one() { 
        assert!(1.is_pm_one());
        assert!((-1).is_pm_one());
        assert!(!2.is_pm_one());
        assert!(!(-2).is_pm_one());
    }
}
