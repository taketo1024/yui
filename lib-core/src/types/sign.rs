use std::ops::Neg;
use derive_more::{Display, Debug};
use is_even::IsEven;
use num_traits::Signed;

#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default, Display, Debug)]
#[cfg_attr(feature = "serde", derive(serde_repr::Serialize_repr, serde_repr::Deserialize_repr))]
#[repr(i8)]
pub enum Sign { 
    #[default]
    #[display("+")]
    #[debug("+")]
    Pos = 1,

    #[display("-")]
    #[debug("-")]
    Neg = -1
}

impl Sign { 
    pub fn is_positive(&self) -> bool { 
        self == &Sign::Pos
    }

    pub fn is_negative(&self) -> bool { 
        !self.is_positive()
    }

    pub fn from_parity<I: IsEven>(val: I) -> Self { 
        if val.is_even() { 
            Sign::Pos
        } else { 
            Sign::Neg
        }
    }
}

macro_rules! impl_int_conversion {
    ($t:tt) => {
        impl From<$t> for Sign {
            fn from(value: $t) -> Self {
                match value { 
                     1 => Sign::Pos,
                    -1 => Sign::Neg,
                     _ => panic!()
                }
            }
        }
        
        impl From<Sign> for $t {
            fn from(value: Sign) -> Self {
                match value { 
                    Sign::Pos =>  1,
                    Sign::Neg => -1
                }
            }
        }                
    };
}

impl_int_conversion!(i8);
impl_int_conversion!(i16);
impl_int_conversion!(i32);
impl_int_conversion!(i64);
impl_int_conversion!(isize);

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

    #[test]
    fn to_string() { 
        assert_eq!(&Sign::Neg.to_string(), "-");
        assert_eq!(&Sign::Pos.to_string(), "+");
    }

    #[cfg(feature = "serde")]
    #[test]
    fn serialize() { 
        let s = Sign::Pos;
        let ser = serde_json::to_string(&s).unwrap();
        assert_eq!(ser, "1");

        let des = serde_json::from_str(&ser).unwrap();
        assert_eq!(s, des);

        let s = Sign::Neg;
        let ser = serde_json::to_string(&s).unwrap();
        assert_eq!(ser, "-1");
        
        let des = serde_json::from_str(&ser).unwrap();
        assert_eq!(s, des);
    }
}