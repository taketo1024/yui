use core::fmt;
use std::fmt::{Display, Debug};
use std::ops::Index;

use derive_more::{Display, DebugCustom};

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, DebugCustom, Display, Default)]
pub enum Bit { 
    #[default]
    #[display(fmt="0")]
    #[debug(fmt="0")]
    Bit0, 

    #[display(fmt="1")]
    #[debug(fmt="0")]
    Bit1
}

impl Bit { 
    pub fn is_zero(&self) -> bool { 
        self == &Bit::Bit0
    }

    pub fn is_one(&self) -> bool { 
        self == &Bit::Bit1
    }

    pub fn as_u64(&self) -> u64 { 
        if self.is_zero() { 0 } else { 1 }
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

#[derive(Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct BitSeq { 
    val: u64,
    len: usize
}

impl BitSeq { 
    pub const MAX_LEN: usize = 64;

    pub fn new(val: u64, len: usize) -> Self { 
        assert!(len <= Self::MAX_LEN);
        assert!(val < (1 << len));
        Self { val, len }
    }

    pub fn new_rev(val: u64, len: usize) -> Self { 
        let val = val.reverse_bits() >> (64 - len);
        Self::new(val, len)
    }

    pub fn empty() -> Self { 
        Self::new(0, 0)
    }

    pub fn zeros(len: usize) -> Self { 
        Self::new(0, len)
    }

    pub fn ones(len: usize) -> Self { 
        let val = (1 << len) - 1;
        Self::new(val, len)
    }

    pub fn len(&self) -> usize { 
        self.len
    }

    pub fn as_u64(&self) -> u64 { 
        self.val
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
        let mut val = self.val;
        
        (0..self.len).into_iter().map(move |_| {
            let b = val & 1;
            val >>= 1;
            Bit::from(b)
        })
    }

    pub fn set(&mut self, i: usize, b: Bit) {
        assert!(i < self.len);
        if b.is_zero() { 
            self.val &= !(1 << i);
        } else { 
            self.val |= 1 << i;
        }
    }

    pub fn set_0(&mut self, i: usize) {
        self.set(i, Bit::Bit0)
    }

    pub fn set_1(&mut self, i: usize) {
        self.set(i, Bit::Bit1)
    }

    pub fn push(&mut self, b: Bit) {
        if b.is_one() { 
            self.val |= 1 << self.len;
        }
        self.len += 1;
    }

    pub fn push_0(&mut self) { 
        self.push(Bit::Bit0)
    }

    pub fn push_1(&mut self) { 
        self.push(Bit::Bit1)
    }

    pub fn append(&mut self, b: BitSeq) {
        assert!(self.len + b.len <= Self::MAX_LEN);
        self.val = b.val << self.len | self.val;
        self.len += b.len;
    }

    pub fn remove(&mut self, i: usize) { 
        assert!(i < self.len);
        
        let a = self.val & !((1 << i + 1) - 1);
        let b = self.val & ((1 << i) - 1);

        self.val = a >> 1 | b;
        self.len -= 1;
    }

    pub fn insert(&mut self, i: usize, b: Bit) { 
        assert!(i <= self.len);
        assert!(self.len < Self::MAX_LEN);

        let mask = (1 << i) - 1;
        let a = self.val & !mask;
        let b = b.as_u64() << i;
        let c = self.val & mask;

        self.val = a << 1 | b | c;
        self.len += 1;
    }

    pub fn insert_0(&mut self, i: usize) { 
        self.insert(i, Bit::Bit0)
    }

    pub fn insert_1(&mut self, i: usize) { 
        self.insert(i, Bit::Bit1)
    }

    pub fn edit<F>(&self, f: F) -> Self
    where F: FnOnce(&mut BitSeq) {
        let mut copy = self.clone();
        f(&mut copy);
        copy
    }

    pub fn sub(&self, l: usize) -> Self { 
        assert!(l <= self.len);
        let val = self.val & ((1 << l) - 1);
        Self::new(val, l)
    }

    pub fn is_sub(&self, other: &Self) -> bool { 
        self.len <= other.len && 
        self.val == (other.val & ((1 << self.len) - 1))
    }

    pub fn generate(len: usize) -> Vec<BitSeq> {
        assert!(len <= Self::MAX_LEN);
        (0..2_u64.pow(len as u32)).map(|v| Self::new(v, len)).collect()
    }
}

impl<T> From<T> for BitSeq 
where Bit: From<T> {
    fn from(b: T) -> Self {
        let val = if Bit::from(b).is_zero() { 0 } else { 1 };
        Self::new(val, 1)
    }
}

impl<T, const N: usize> From<[T; N]> for BitSeq 
where Bit: From<T> {
    fn from(value: [T; N]) -> Self {
        Self::from_iter(value)
    }
}

impl<T> FromIterator<T> for BitSeq
where Bit: From<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let mut val = 0;
        let mut len = 0;
        for b in iter.into_iter() {
            if Bit::from(b).is_one() {
                val |= 1 << len;
            }
            len += 1;
        }
        Self::new(val, len)
    }
}

impl Index<usize> for BitSeq {
    type Output = Bit;

    fn index(&self, i: usize) -> &Self::Output {
        assert!(i < self.len);
        match (self.val >> i) & 1 { 
            0 => &Bit::Bit0,
            1 => &Bit::Bit1,
            _ => panic!()
        }
    }
}

impl Display for BitSeq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for b in self.iter() { 
            Display::fmt(&b, f)?;
        }
        Ok(())
    }
}

impl Debug for BitSeq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        Display::fmt(self, f)
    }
}


impl PartialOrd for BitSeq {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(&other))
    }
}

// TODO support lex-order (using generic parameter).

impl Ord for BitSeq {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.len().cmp(&other.len()).then_with( ||
            self.weight().cmp(&other.weight())
        ).then_with(|| 
            self.as_u64().cmp(&other.as_u64())
        )
    }
}

#[cfg(test)]
mod tests { 
    use Bit::*;
    use itertools::Itertools;
    use super::*;

    #[test]
    fn new() { 
        let b = BitSeq::new(0b10110, 5);
        assert_eq!(b.val, 22);
        assert_eq!(b.len, 5);
    }

    #[test]
    fn new_rev() { 
        let b = BitSeq::new_rev(0b01101, 5);
        assert_eq!(b.val, 22);
        assert_eq!(b.len, 5);
    }

    #[test]
    fn from_arr() { 
        let b = BitSeq::from([1,0,1,1,0]);
        assert_eq!(b, BitSeq::new(0b01101, 5));
    }

    #[test]
    fn from_iter() { 
        let b = BitSeq::from_iter([1,0,1,1,0]);
        assert_eq!(b, BitSeq::new(0b01101, 5));
    }

    #[test]
    fn weight() { 
        let b = BitSeq::new(0b10110, 5);
        assert_eq!(b.weight(), 3);
    }

    #[test]
    fn index() { 
        let b = BitSeq::new(0b01101, 5);
        assert_eq!(b.len(), 5);
        assert_eq!(b[0], Bit1);
        assert_eq!(b[1], Bit0);
        assert_eq!(b[2], Bit1);
        assert_eq!(b[3], Bit1);
        assert_eq!(b[4], Bit0);
    }

    #[test]
    fn iter() { 
        let b = BitSeq::new(0b01101, 5);
        let v = b.iter().collect_vec();
        assert_eq!(v, vec![Bit1, Bit0, Bit1, Bit1, Bit0])
    }

    #[test]
    fn to_string() { 
        let b = BitSeq::new(0b01101, 5);
        let s = b.to_string();
        assert_eq!(s, "10110");
    }

    #[test]
    fn set() { 
        let mut b = BitSeq::new(0b01101, 5);

        b.set(0, Bit1);
        assert_eq!(b, BitSeq::new(0b01101, 5));

        b.set(1, Bit0);
        assert_eq!(b, BitSeq::new(0b01101, 5));

        b.set(2, Bit0);
        assert_eq!(b, BitSeq::new(0b01001, 5));

        b.set(3, Bit0);
        assert_eq!(b, BitSeq::new(0b00001, 5));

        b.set(4, Bit1);
        assert_eq!(b, BitSeq::new(0b10001, 5));
    }

    #[test]
    fn remove() { 
        let mut b = BitSeq::new(0b100101, 6);

        b.remove(0);
        assert_eq!(b, BitSeq::new(0b10010, 5));

        b.remove(2);
        assert_eq!(b, BitSeq::new(0b1010, 4));

        b.remove(3);
        assert_eq!(b, BitSeq::new(0b010, 3));

        b.remove(2);
        assert_eq!(b, BitSeq::new(0b10, 2));
        
        b.remove(1);
        assert_eq!(b, BitSeq::new(0b0, 1));
        
        b.remove(0);
        assert_eq!(b, BitSeq::new(0b0, 0));
    }

    #[test]
    fn insert() { 
        let mut b = BitSeq::empty();

        b.insert(0, Bit1);
        assert_eq!(b, BitSeq::new(0b1, 1));

        b.insert(0, Bit0);
        assert_eq!(b, BitSeq::new(0b10, 2));

        b.insert(2, Bit0);
        assert_eq!(b, BitSeq::new(0b010, 3));

        b.insert(3, Bit1);
        assert_eq!(b, BitSeq::new(0b1010, 4));
    }

    #[test]
    fn push() { 
        let mut b = BitSeq::new(0b01101, 5);

        b.push(Bit0);
        assert_eq!(b, BitSeq::new(0b001101, 6));

        b.push(Bit1);
        assert_eq!(b, BitSeq::new(0b1001101, 7));
    }

    #[test]
    fn append() { 
        let mut b0 = BitSeq::new(0b10110, 5);
        let b1 = BitSeq::new(0b0101, 4);

        b0.append(b1);

        assert_eq!(b0, BitSeq::new(0b010110110, 9));
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
        // order priority: len > weight > val

        let b0 = BitSeq::new(0b0,  1);
        let b1 = BitSeq::new(0b00, 2);

        assert!(b0 < b1);
        
        let b0 = BitSeq::new(0b110, 3);
        let b1 = BitSeq::new(0b100, 3);
        let b2 = BitSeq::new(0b011, 3);

        assert!(b0 > b1);
        assert!(b1 < b2);
        assert!(b0 > b2);
    }

    #[test]
    fn sub() { 
        let b = BitSeq::new(0b10110, 5);

        assert_eq!(b.sub(0), BitSeq::empty());
        assert_eq!(b.sub(3), BitSeq::new(0b110, 3));
        assert_eq!(b.sub(5), b);
    }

    #[test]
    fn is_sub() { 
        let b0 = BitSeq::new(0b110,   3);
        let b1 = BitSeq::new(0b10110, 5);
        let b2 = BitSeq::new(0b11110, 5);

        assert!(b0.is_sub(&b1));
        assert!(b0.is_sub(&b2));
        assert!(!b1.is_sub(&b0));
        assert!(!b1.is_sub(&b2));
        assert!(!b2.is_sub(&b0));
        assert!(!b2.is_sub(&b1));
    }

    #[test]
    fn edit() {
        let b = BitSeq::new(0b10110, 5);
        let c = b.edit(|b| b.set_1(0));
        assert_eq!(c, BitSeq::new(0b10111, 5))
    }
}