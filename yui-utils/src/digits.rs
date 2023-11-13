use num_traits::Zero;

pub trait Digits { 
    fn digits(&self) -> Vec<u8> { 
        self.digits_rev().into_iter().rev().collect()
    }

    fn digits_rev(&self) -> Vec<u8> { 
        self.digits().into_iter().rev().collect()
    }
}

impl Digits for usize { 
    fn digits_rev(&self) -> Vec<u8> {
        if self.is_zero() { return vec![0] }

        let mut num = *self;
        (0..).map_while(|_| { 
            if num > 0 { 
                let d = (num % 10) as u8;
                num /= 10;
                Some(d)
            } else {
                None
            }
        }).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn digits() { 
        assert_eq!(0.digits(), [0]);
        assert_eq!(123.digits(), [1, 2, 3]);
    }
}
