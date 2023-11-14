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
    fn from_sign(s: Sign) -> Self { 
        Self::from(s.to_i32())
    }

    fn inv(&self) -> Option<Self>;
    fn is_unit(&self) -> bool;
    fn normalizing_unit(&self) -> Self;

    fn normalized(&self) -> Self { 
        self.clone().into_normalized()
    }

    fn into_normalized(self) -> Self { 
        let u = self.normalizing_unit();
        if u.is_one() { 
            self
        } else { 
            self * u
        }
    }
    
    fn is_pm_one(&self) -> bool { 
        self.is_one() || (-self).is_one()
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

    #[test]
    fn normalized() {
        assert_eq!(3.normalized(), 3);
        assert_eq!((-3).normalized(), 3);
    } 

}