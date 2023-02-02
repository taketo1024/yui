use num_traits::Signed;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sign { 
    Pos, Neg
}

impl Sign { 
    pub fn is_positive(&self) -> bool { 
        self == &Sign::Pos
    }

    pub fn is_negative(&self) -> bool { 
        self == &Sign::Neg
    }
}

impl<T> From<T> for Sign 
where T: Signed {
    fn from(i: T) -> Self {
        assert!(!i.is_zero());
        if i.is_positive() { Sign::Pos } else { Sign::Neg }
    }
}