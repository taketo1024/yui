pub trait IntoDigits: Sized { 
    type Digit;
    fn into_digits(self) -> Vec<Self::Digit> { 
        self.into_rev_digits().into_iter().rev().collect()
    }

    fn into_rev_digits(self) -> Vec<Self::Digit> { 
        self.into_digits().into_iter().rev().collect()
    }
}

macro_rules! impl_into_digits {
    ($t: ty, $d: ty) => {
        impl IntoDigits for $t { 
            type Digit = $d;
            fn into_rev_digits(self) -> Vec<$d> {
                if self == 0 { return vec![0] }
        
                let mut num = self;
                (0..).map_while(|_| { 
                    if num > 0 { 
                        let d = (num % 10) as $d;
                        num /= 10;
                        Some(d)
                    } else {
                        None
                    }
                }).collect()
            }
        }                
    };
}

impl_into_digits!(usize, u8);

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn usize() { 
        let a = 123456789;
        assert_eq!(a.into_digits(), vec![1,2,3,4,5,6,7,8,9]);
        assert_eq!(a.into_rev_digits(), vec![9,8,7,6,5,4,3,2,1]);
    }
}