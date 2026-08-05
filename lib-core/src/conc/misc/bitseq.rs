//! A single [`Bit`], and a packed sequence of bits [`BitSeq`] over a generic word `I`.

use core::fmt;
use std::fmt::{Display, Debug};
use std::hash::Hash;
use std::ops::{Add, AddAssign, BitAnd, BitAndAssign, BitOr, BitOrAssign, Index, Not, Shl, Shr, ShrAssign, Sub};
use std::str::FromStr;
use auto_impl_ops::auto_ops;

/// A single binary digit, `0` or `1`.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default, derive_more::Display, derive_more::Debug)]
#[cfg_attr(feature = "serde", derive(serde_repr::Serialize_repr, serde_repr::Deserialize_repr))]
#[repr(u8)]
pub enum Bit {
    #[default]
    #[display("0")]
    #[debug("0")]
    Bit0 = 0,

    #[display("1")]
    #[debug("1")]
    Bit1 = 1
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
impl_bit_from_int!(u128);
impl_bit_from_int!(usize);
impl_bit_from_int!(i8);
impl_bit_from_int!(i16);
impl_bit_from_int!(i32);
impl_bit_from_int!(i64);
impl_bit_from_int!(isize);

/// The unsigned word backing a [`BitSeq`]. Implemented for `u8` … `u128`.
pub trait BitRepr:
    Copy + Eq + Ord + Hash + Default + Debug + Send + Sync + 'static
    + BitAnd<Output = Self> + BitAndAssign
    + BitOr<Output = Self> + BitOrAssign
    + Not<Output = Self>
    + Shl<usize, Output = Self> + Shr<usize, Output = Self> + ShrAssign<usize>
    + Sub<Output = Self>
{
    const BITS: usize;
    const ZERO: Self;
    const ONE: Self;
    const MAX: Self;

    fn count_ones(self) -> u32;
    fn reverse_bits(self) -> Self;
    fn from_u128(v: u128) -> Self;
    fn to_u128(self) -> u128;
}

macro_rules! impl_bit_repr {
    ($($t:ty),* $(,)?) => {$(
        impl BitRepr for $t {
            const BITS: usize = <$t>::BITS as usize;
            const ZERO: Self = 0;
            const ONE: Self = 1;
            const MAX: Self = <$t>::MAX;

            fn count_ones(self) -> u32 {
                <$t>::count_ones(self)
            }

            fn reverse_bits(self) -> Self {
                <$t>::reverse_bits(self)
            }

            fn from_u128(v: u128) -> Self {
                v as $t
            }

            fn to_u128(self) -> u128 {
                self as u128
            }
        }
    )*};
}

impl_bit_repr!(u8, u16, u32, u64, u128);

/// A sequence of [`Bit`]s of length up to [`MAX_LEN`](BitSeq::MAX_LEN) = `I::BITS`,
/// packed into a single word `I` (default `u64`). The bit at index `i` is stored at
/// position `i` of `val`, i.e. the least-significant bit of `val` is `self[0]`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde_with::SerializeDisplay, serde_with::DeserializeFromStr))]
pub struct BitSeq<I: BitRepr = u64> {
    val: I,
    len: usize
}

impl<I: BitRepr> BitSeq<I> {
    pub const MAX_LEN: usize = I::BITS;

    pub fn new(val: I, len: usize) -> Self {
        assert!(len <= Self::MAX_LEN);
        assert!(len == Self::MAX_LEN || val < (I::ONE << len));
        Self { val, len }
    }

    /// Like [`new`](Self::new), but interprets `val` as bit-reversed so that
    /// the most-significant bit becomes `self[0]`.
    pub fn new_rev(val: I, len: usize) -> Self {
        if len == 0 {
            return Self::empty();
        }
        let val = val.reverse_bits() >> (Self::MAX_LEN - len);
        Self::new(val, len)
    }

    pub fn empty() -> Self {
        Self::new(I::ZERO, 0)
    }

    pub fn zeros(len: usize) -> Self {
        Self::new(I::ZERO, len)
    }

    pub fn ones(len: usize) -> Self {
        let val = if len == 0 { I::ZERO } else { I::MAX >> (Self::MAX_LEN - len) };
        Self::new(val, len)
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn as_u128(&self) -> u128 {
        self.val.to_u128()
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// The number of bits set to `1` (Hamming weight).
    pub fn weight(&self) -> usize {
        self.val.count_ones() as usize
    }

    pub fn iter(&self) -> impl Iterator<Item = Bit> + use<I> {
        let mut val = self.val;

        (0..self.len).map(move |_| {
            let b = val & I::ONE == I::ONE;
            val >>= 1;
            Bit::from(b)
        })
    }

    pub fn set(&mut self, i: usize, b: Bit) {
        assert!(i < self.len);
        if b.is_zero() {
            self.val &= !(I::ONE << i);
        } else {
            self.val |= I::ONE << i;
        }
    }

    pub fn set_0(&mut self, i: usize) {
        self.set(i, Bit::Bit0)
    }

    pub fn set_1(&mut self, i: usize) {
        self.set(i, Bit::Bit1)
    }

    pub fn push(&mut self, b: Bit) {
        assert!(self.len < Self::MAX_LEN);
        if b.is_one() {
            self.val |= I::ONE << self.len;
        }
        self.len += 1;
    }

    pub fn push_0(&mut self) {
        self.push(Bit::Bit0)
    }

    pub fn push_1(&mut self) {
        self.push(Bit::Bit1)
    }

    pub fn append(&mut self, b: BitSeq<I>) {
        assert!(self.len + b.len <= Self::MAX_LEN);

        // `self.len` may be `MAX_LEN`, where the shift below is undefined.
        if b.len == 0 {
            return
        }

        self.val |= b.val << self.len;
        self.len += b.len;
    }

    pub fn remove(&mut self, i: usize) {
        assert!(i < self.len);

        // shifted in two steps: `i + 1` may be `MAX_LEN`, where a shift is undefined.
        let hi = ((self.val >> i) >> 1) << i;
        let lo = self.val & ((I::ONE << i) - I::ONE);

        self.val = hi | lo;
        self.len -= 1;
    }

    pub fn insert(&mut self, i: usize, b: Bit) {
        assert!(i <= self.len);
        assert!(self.len < Self::MAX_LEN);

        let mask = (I::ONE << i) - I::ONE;
        let a = self.val & !mask;
        let b = if b.is_one() { I::ONE << i } else { I::ZERO };
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
    where F: FnOnce(&mut BitSeq<I>) {
        let mut copy = *self;
        f(&mut copy);
        copy
    }

    pub fn sub(&self, l: usize) -> Self {
        assert!(l <= self.len);
        let val = if l == Self::MAX_LEN { self.val } else { self.val & ((I::ONE << l) - I::ONE) };
        Self::new(val, l)
    }

    /// `true` if `self` is a prefix of `other`.
    pub fn is_sub(&self, other: &Self) -> bool {
        self.len <= other.len &&
        self.val == other.sub(self.len).val
    }

    /// Enumerate all `2^len` sequences of the given length, in ascending order of `val`.
    pub fn generate(len: usize) -> impl Iterator<Item = BitSeq<I>> {
        assert!(len < 128, "generate is only sensible for small lengths");
        assert!(len <= Self::MAX_LEN);
        (0 .. (1_u128 << len)).map(move |v| Self::new(I::from_u128(v), len))
    }
}

impl<I: BitRepr, T> From<T> for BitSeq<I>
where Bit: From<T> {
    fn from(b: T) -> Self {
        let val = if Bit::from(b).is_zero() { I::ZERO } else { I::ONE };
        Self::new(val, 1)
    }
}

impl<I: BitRepr, T, const N: usize> From<[T; N]> for BitSeq<I>
where Bit: From<T> {
    fn from(value: [T; N]) -> Self {
        Self::from_iter(value)
    }
}

impl<I: BitRepr, T> FromIterator<T> for BitSeq<I>
where Bit: From<T> {
    fn from_iter<Itr: IntoIterator<Item = T>>(iter: Itr) -> Self {
        let mut val = I::ZERO;
        let mut len = 0;
        for b in iter.into_iter() {
            if Bit::from(b).is_one() {
                val |= I::ONE << len;
            }
            len += 1;
        }
        Self::new(val, len)
    }
}

impl<I: BitRepr> FromStr for BitSeq<I> {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.chars().map(|c|
            match c {
                '0' => Ok(Bit::Bit0),
                '1' => Ok(Bit::Bit1),
                _   => Err("Invalid bit: {c}".into())
            }
        ).collect()
    }
}

impl<I: BitRepr> Index<usize> for BitSeq<I> {
    type Output = Bit;

    fn index(&self, i: usize) -> &Self::Output {
        assert!(i < self.len);
        if (self.val >> i) & I::ONE == I::ONE {
            &Bit::Bit1
        } else {
            &Bit::Bit0
        }
    }
}

impl<I: BitRepr> Display for BitSeq<I> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for b in self.iter() {
            Display::fmt(&b, f)?;
        }
        Ok(())
    }
}

impl<I: BitRepr> Debug for BitSeq<I> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        Display::fmt(self, f)
    }
}

impl<I: BitRepr> PartialOrd for BitSeq<I> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

// TODO support lex-order (using generic parameter).

/// Ordered by `len`, then [`weight`](BitSeq::weight), then raw `val`.
impl<I: BitRepr> Ord for BitSeq<I> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.len().cmp(&other.len()).then_with( ||
            self.weight().cmp(&other.weight())
        ).then_with(||
            self.val.cmp(&other.val)
        )
    }
}

#[auto_ops]
impl<I: BitRepr> AddAssign<&BitSeq<I>> for BitSeq<I> {
    fn add_assign(&mut self, rhs: &Self) {
        self.append(*rhs);
    }
}

#[auto_ops]
impl<I: BitRepr> AddAssign<Bit> for BitSeq<I> {
    fn add_assign(&mut self, b: Bit) {
        self.push(b);
    }
}

#[cfg(test)]
mod tests {
    use Bit::*;
    use itertools::Itertools;
    use super::*;

    type B = BitSeq;

    #[test]
    fn new() {
        let b = B::new(0b10110, 5);
        assert_eq!(b.val, 22);
        assert_eq!(b.len, 5);
    }

    #[test]
    fn new_rev() {
        let b = B::new_rev(0b01101, 5);
        assert_eq!(b.val, 22);
        assert_eq!(b.len, 5);
    }

    #[test]
    fn from_arr() {
        let b = B::from([1,0,1,1,0]);
        assert_eq!(b, B::new(0b01101, 5));
    }

    #[test]
    fn from_iter() {
        let b = B::from_iter([1,0,1,1,0]);
        assert_eq!(b, B::new(0b01101, 5));
    }

    #[test]
    fn weight() {
        let b = B::new(0b10110, 5);
        assert_eq!(b.weight(), 3);

        let b = B::new(0b0110101101, 10);
        assert_eq!(b.weight(), 6);
    }

    #[test]
    fn index() {
        let b = B::new(0b01101, 5);
        assert_eq!(b.len(), 5);
        assert_eq!(b[0], Bit1);
        assert_eq!(b[1], Bit0);
        assert_eq!(b[2], Bit1);
        assert_eq!(b[3], Bit1);
        assert_eq!(b[4], Bit0);
    }

    #[test]
    fn iter() {
        let b = B::new(0b01101, 5);
        let v = b.iter().collect_vec();
        assert_eq!(v, vec![Bit1, Bit0, Bit1, Bit1, Bit0])
    }

    #[test]
    fn to_string() {
        let b = B::new(0b01101, 5);
        let s = b.to_string();
        assert_eq!(s, "10110");
    }

    #[test]
    fn set() {
        let mut b = B::new(0b01101, 5);

        b.set(0, Bit1);
        assert_eq!(b, B::new(0b01101, 5));

        b.set(1, Bit0);
        assert_eq!(b, B::new(0b01101, 5));

        b.set(2, Bit0);
        assert_eq!(b, B::new(0b01001, 5));

        b.set(3, Bit0);
        assert_eq!(b, B::new(0b00001, 5));

        b.set(4, Bit1);
        assert_eq!(b, B::new(0b10001, 5));
    }

    #[test]
    fn remove() {
        let mut b = B::new(0b100101, 6);

        b.remove(0);
        assert_eq!(b, B::new(0b10010, 5));

        b.remove(2);
        assert_eq!(b, B::new(0b1010, 4));

        b.remove(3);
        assert_eq!(b, B::new(0b010, 3));

        b.remove(2);
        assert_eq!(b, B::new(0b10, 2));

        b.remove(1);
        assert_eq!(b, B::new(0b0, 1));

        b.remove(0);
        assert_eq!(b, B::new(0b0, 0));
    }

    #[test]
    fn remove_at_max_len() {
        // dropping the top bit of a full-length sequence.
        let n = B::MAX_LEN;

        let mut b = B::ones(n);
        b.remove(n - 1);
        assert_eq!(b, B::ones(n - 1));

        let mut b = B::zeros(n);
        b.set_1(n - 1);
        b.remove(n - 1);
        assert_eq!(b, B::zeros(n - 1));
    }

    #[test]
    fn insert() {
        let mut b = B::empty();

        b.insert(0, Bit1);
        assert_eq!(b, B::new(0b1, 1));

        b.insert(0, Bit0);
        assert_eq!(b, B::new(0b10, 2));

        b.insert(2, Bit0);
        assert_eq!(b, B::new(0b010, 3));

        b.insert(3, Bit1);
        assert_eq!(b, B::new(0b1010, 4));
    }

    #[test]
    fn push() {
        let mut b = B::new(0b01101, 5);

        b.push(Bit0);
        assert_eq!(b, B::new(0b001101, 6));

        b.push(Bit1);
        assert_eq!(b, B::new(0b1001101, 7));
    }

    #[test]
    fn push_by_add() {
        let mut b = B::new(0b01101, 5);

        b += Bit::Bit0;
        assert_eq!(b, B::new(0b001101, 6));

        b += Bit::Bit1;
        assert_eq!(b, B::new(0b1001101, 7));
    }

    #[test]
    fn append() {
        let mut b0 = B::new(0b10110, 5);
        let b1 = B::new(0b0101, 4);

        b0.append(b1);

        assert_eq!(b0, B::new(0b010110110, 9));
    }

    #[test]
    fn append_empty_to_full() {
        let full = B::ones(B::MAX_LEN);

        let mut b = full;
        b.append(B::empty());

        assert_eq!(b, full);
    }

    #[test]
    fn append_by_add() {
        let mut b0 = B::new(0b10110, 5);
        let b1 = B::new(0b0101, 4);

        b0 += b1;

        assert_eq!(b0, B::new(0b010110110, 9));
    }

    #[test]
    fn generate() {
        let v = B::generate(3).collect_vec();
        assert_eq!(v, vec![
            B::new(0b000, 3),
            B::new(0b001, 3),
            B::new(0b010, 3),
            B::new(0b011, 3),
            B::new(0b100, 3),
            B::new(0b101, 3),
            B::new(0b110, 3),
            B::new(0b111, 3),
        ]);
    }

    #[test]
    fn ord() {
        // order priority: len > weight > val

        let b0 = B::new(0b0,  1);
        let b1 = B::new(0b00, 2);

        assert!(b0 < b1);

        let b0 = B::new(0b110, 3);
        let b1 = B::new(0b100, 3);
        let b2 = B::new(0b011, 3);

        assert!(b0 > b1);
        assert!(b1 < b2);
        assert!(b0 > b2);
    }

    #[test]
    fn sub() {
        let b = B::new(0b10110, 5);

        assert_eq!(b.sub(0), B::empty());
        assert_eq!(b.sub(3), B::new(0b110, 3));
        assert_eq!(b.sub(5), b);
    }

    #[test]
    fn is_sub() {
        let b0 = B::new(0b110,   3);
        let b1 = B::new(0b10110, 5);
        let b2 = B::new(0b11110, 5);

        assert!(b0.is_sub(&b1));
        assert!(b0.is_sub(&b2));
        assert!(!b1.is_sub(&b0));
        assert!(!b1.is_sub(&b2));
        assert!(!b2.is_sub(&b0));
        assert!(!b2.is_sub(&b1));
    }

    #[test]
    fn edit() {
        let b = B::new(0b10110, 5);
        let c = b.edit(|b| b.set_1(0));
        assert_eq!(c, B::new(0b10111, 5))
    }

    #[test]
    fn u128_long() {
        type B128 = BitSeq<u128>;

        assert_eq!(B128::MAX_LEN, 128);

        let mut b = B128::zeros(100);
        assert_eq!(b.len(), 100);
        assert_eq!(b.weight(), 0);

        b.set_1(72);
        b.set_1(99);
        assert_eq!(b.weight(), 2);
        assert_eq!(b[72], Bit1);
        assert_eq!(b[71], Bit0);

        let ones = B128::ones(128);
        assert_eq!(ones.len(), 128);
        assert_eq!(ones.weight(), 128);

        let s = b.to_string();
        assert_eq!(s.len(), 100);
        assert_eq!(B128::from_str(&s).unwrap(), b);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn serialize() {
        let b = B::new(0b10110, 5);
        let ser = serde_json::to_string(&b).unwrap();
        assert_eq!(ser, "\"01101\"");

        let des = serde_json::from_str(&ser).unwrap();
        assert_eq!(b, des);
    }
}
