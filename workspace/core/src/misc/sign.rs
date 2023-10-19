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

    pub fn to_i32(&self) -> i32 { 
        if self.is_positive() { 1 } else { -1 }
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