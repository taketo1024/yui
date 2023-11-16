use std::ops::Neg;
use num_traits::Signed;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Sign { 
    Neg, 
    Pos
}


impl Sign { 
    pub fn is_positive(&self) -> bool { 
        self == &Sign::Pos
    }

    pub fn is_negative(&self) -> bool { 
        !self.is_positive()
    }
}

impl From<i32> for Sign {
    fn from(value: i32) -> Self {
        match value { 
             1 => Sign::Pos,
            -1 => Sign::Neg,
             _ => panic!()
        }
    }
}

impl From<Sign> for i32 {
    fn from(value: Sign) -> Self {
        match value { 
            Sign::Pos =>  1,
            Sign::Neg => -1
        }
    }
}

impl Neg for Sign {
    type Output = Self;
    fn neg(self) -> Self {
        use Sign::*;
        match self { 
            Neg => Pos,
            Pos => Neg
        }
    }
}

pub trait GetSign { 
    fn sign(&self) -> Sign;
}

impl<T> GetSign for T where T: Signed {
    fn sign(&self) -> Sign {
        if self.is_positive() { 
            Sign::Pos 
        } else { 
            Sign::Neg 
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ord() {
        assert!(Sign::Neg < Sign::Pos)
    }
}