use core::fmt;
use std::fmt::Display;
use std::ops::Index;

use derive_more::Display;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Display, Default)]
pub enum Bit { 
    #[default]
    #[display(fmt="0")]
    Bit0, 

    #[display(fmt="1")]
    Bit1
}

impl Bit { 
    pub fn is_zero(&self) -> bool { 
        self == &Bit::Bit0
    }

    pub fn is_one(&self) -> bool { 
        self == &Bit::Bit1
    }
}

impl From<bool> for Bit {
    fn from(b: bool) -> Self {
        if b { 
            Bit::Bit1
        } else { 
            Bit::Bit0
        }
    }
}

macro_rules! impl_bit_from_int {
    ($t:ty) => {
        impl From<$t> for Bit {
            fn from(val: $t) -> Self {
                match val { 
                    0 => Bit::Bit0,
                    1 => Bit::Bit1,
                    _ => panic!()
                }
            }
        }
    };
}

impl_bit_from_int!(u8);
impl_bit_from_int!(u16);
impl_bit_from_int!(u32);
impl_bit_from_int!(u64);
impl_bit_from_int!(usize);
impl_bit_from_int!(i8);
impl_bit_from_int!(i16);
impl_bit_from_int!(i32);
impl_bit_from_int!(i64);
impl_bit_from_int!(isize);

#[derive(Clone, Copy, PartialEq, Eq, Ord, Hash, Debug, Default)]
pub struct BitSeq { 
    val: u64,
    len: usize
}

impl BitSeq { 
    pub const MAX_LEN: usize = 64;

    pub fn new(val: u64, len: usize) -> Self { 
        assert!(len <= Self::MAX_LEN);
        Self { val, len }
    }

    pub fn empty() -> Self { 
        Self::new(0, 0)
    }

    pub fn len(&self) -> usize { 
        self.len
    }

    pub fn is_empty(&self) -> bool { 
        self.len == 0
    } 

    pub fn weight(&self) -> usize { 
        let mut v = self.val;
        (0..self.len).filter(|_| {
            let b = (v & 1) == 1;
            v >>= 1;
            b
        }).count()
    }

    pub fn iter(&self) -> impl Iterator<Item = Bit> {
        let val = self.val;
        (0..self.len).rev().into_iter().map(move |i| {
            Bit::from((val >> i) & 1)
        })
    }

    pub fn set(&mut self, i: usize, b: Bit) {
        assert!(i < self.len);
        if self[i] == b { 
            return
        }

        let j = self.len - i - 1;

        if b.is_zero() { 
            self.val &= !(1 << j);
        } else { 
            self.val |= 1 << j;
        }
    }

    pub fn set_0(&mut self, i: usize) {
        self.set(i, Bit::Bit0)
    }

    pub fn set_1(&mut self, i: usize) {
        self.set(i, Bit::Bit1)
    }

    pub fn push(&mut self, b: Bit) {
        self.val = if b.is_zero() { 
            self.val << 1
        } else { 
            self.val << 1 | 1
        };
        self.len += 1;
    }

    pub fn append(&mut self, b: BitSeq) {
        assert!(self.len + b.len <= Self::MAX_LEN);
        self.val = self.val << b.len | b.val;
        self.len += b.len;
    }

    pub fn generate(len: usize) -> Vec<BitSeq> {
        assert!(len <= Self::MAX_LEN);
        (0..2_u64.pow(len as u32)).map(|v| Self::new(v, len)).collect()
    }
}

impl<T> FromIterator<T> for BitSeq
where Bit: From<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let mut val = 0;
        let mut len = 0;
        for b in iter.into_iter() {
            val <<= 1;
            len += 1;
            if Bit::from(b).is_one() {
                val |= 1;
            }
        }
        Self::new(val, len)
    }
}

impl Index<usize> for BitSeq {
    type Output = Bit;

    fn index(&self, i: usize) -> &Self::Output {
        assert!(i < self.len);
        let j = self.len - i - 1;
        match (self.val >> j) & 1 { 
            0 => &Bit::Bit0,
            1 => &Bit::Bit1,
            _ => panic!()
        }
    }
}

impl Display for BitSeq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.len { 
            self[i].fmt(f)?;
        }
        Ok(())
    }
}

impl PartialOrd for BitSeq {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let cmp = self.len.cmp(&other.len).then_with( || {
            self.val.cmp(&other.val)
        });
        Some(cmp)
    }
}

#[cfg(test)]
mod tests { 
    use Bit::*;
    use itertools::Itertools;
    use super::*;

    #[test]
    fn from_iter() { 
        let b = BitSeq::from_iter([0,1,1,0,1]);
        assert_eq!(b, BitSeq::new(0b01101, 5));
    }

    #[test]
    fn weight() { 
        let b = BitSeq::new(0b10110, 5);
        assert_eq!(b.weight(), 3);
    }

    #[test]
    fn index() { 
        let b = BitSeq::new(0b10110, 5);
        assert_eq!(b.len(), 5);
        assert_eq!(b[0], Bit1);
        assert_eq!(b[1], Bit0);
        assert_eq!(b[2], Bit1);
        assert_eq!(b[3], Bit1);
        assert_eq!(b[4], Bit0);
    }

    #[test]
    fn iter() { 
        let b = BitSeq::new(0b10110, 5);
        let v = b.iter().collect_vec();
        assert_eq!(v, vec![Bit1, Bit0, Bit1, Bit1, Bit0])
    }

    #[test]
    fn to_string() { 
        let b = BitSeq::new(0b10110, 5);
        let s = b.to_string();
        assert_eq!(s, "10110");
    }

    #[test]
    fn set() { 
        let mut b = BitSeq::new(0b10110, 5);

        b.set(0, Bit0);
        assert_eq!(b, BitSeq::new(0b00110, 5));

        b.set(2, Bit0);
        assert_eq!(b, BitSeq::new(0b00010, 5));

        b.set(3, Bit1);
        assert_eq!(b, BitSeq::new(0b00010, 5));

        b.set(4, Bit1);
        assert_eq!(b, BitSeq::new(0b00011, 5));
    }

    #[test]
    fn push() { 
        let mut b = BitSeq::new(0b10110, 5);

        b.push(Bit0);
        assert_eq!(b, BitSeq::new(0b101100, 6));

        b.push(Bit1);
        assert_eq!(b, BitSeq::new(0b1011001, 7));
    }

    #[test]
    fn append() { 
        let mut b0 = BitSeq::new(0b10110, 5);
        let b1 = BitSeq::new(0b0101, 4);

        b0.append(b1);

        assert_eq!(b0, BitSeq::new(0b101100101, 9));
    }

    #[test]
    fn generate() { 
        let v = BitSeq::generate(3);
        assert_eq!(v, vec![
            BitSeq::new(0b000, 3),
            BitSeq::new(0b001, 3),
            BitSeq::new(0b010, 3),
            BitSeq::new(0b011, 3),
            BitSeq::new(0b100, 3),
            BitSeq::new(0b101, 3),
            BitSeq::new(0b110, 3),
            BitSeq::new(0b111, 3),
        ]);
    }

    #[test]
    fn ord() { 
        let b0 = BitSeq::new(0b110, 3);
        let b1 = BitSeq::new(0b100, 3);
        let b2 = BitSeq::new(0b011, 3);

        assert!(b0 > b1);
        assert!(b1 > b2);
        assert!(b0 > b2);

        let b0 = BitSeq::new(0b0,  1);
        let b1 = BitSeq::new(0b00, 2);
        assert!(b0 < b1);
    }
}