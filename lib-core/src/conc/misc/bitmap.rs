//! Compact bitmap over a small element type. `BitMap<E, S>` is a single
//! word `S` (the storage) whose bit `e.into()` is set iff `e` is present.
//! `S = u128` covers element indices `0..128`, and [`U256`](super::u256::U256)
//! doubles that; extend by impl-ing [`BitStorage`] for a wider type.

use std::marker::PhantomData;
use std::ops::{BitAnd, BitOr, BitOrAssign, Shl, Sub};

/// Single-word storage backing a [`BitMap`]. Implemented for `u8`..`u128` and
/// [`U256`](super::u256::U256).
pub trait BitStorage:
    Copy + Eq + Default
    + BitAnd<Output = Self> + BitOr<Output = Self> + BitOrAssign
    + Shl<u32, Output = Self> + Sub<Output = Self>
{
    /// Number of bits this storage can hold (capacity of the bitmap).
    const WIDTH: u32;

    fn one() -> Self;
    fn is_zero(self) -> bool;
    fn count_ones(self) -> u32;
    fn trailing_zeros(self) -> u32;
}

macro_rules! impl_bit_storage {
    ($($t:ty),+ $(,)?) => {
        $(
            impl BitStorage for $t {
                const WIDTH: u32 = <$t>::BITS;
                fn one() -> Self { 1 }
                fn is_zero(self) -> bool { self == 0 }
                fn count_ones(self) -> u32 { <$t>::count_ones(self) }
                fn trailing_zeros(self) -> u32 { <$t>::trailing_zeros(self) }
            }
        )+
    };
}

impl_bit_storage!(u8, u16, u32, u64, u128);

#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug, Default)]
pub struct BitMap<E, S> {
    bits: S,
    _phantom: PhantomData<E>,
}

impl<E, S> BitMap<E, S>
where
    E: Copy + Into<u32> + TryFrom<u32>,
    S: BitStorage,
{
    pub fn new() -> Self {
        Self { bits: S::default(), _phantom: PhantomData }
    }

    pub fn insert(&mut self, e: E) {
        let bit = e.into();
        // note: not a `debug_assert` — in release the shift is masked and the write wraps silently.
        assert!(bit < S::WIDTH, "BitMap index {bit} ≥ storage width {}", S::WIDTH);
        self.bits |= S::one() << bit;
    }

    pub fn contains(&self, e: E) -> bool {
        let bit = e.into();
        // hot read path; every present bit already went through the checked `insert`.
        debug_assert!(bit < S::WIDTH, "BitMap index {bit} ≥ storage width {}", S::WIDTH);
        !(self.bits & (S::one() << bit)).is_zero()
    }

    pub fn len(&self) -> usize {
        self.bits.count_ones() as usize
    }

    pub fn is_empty(&self) -> bool {
        self.bits.is_zero()
    }

    /// Iterate elements in ascending order. The smallest element (if any) is
    /// `self.iter().next()`.
    pub fn iter(&self) -> Iter<E, S> {
        Iter { bits: self.bits, _phantom: PhantomData }
    }
}

impl<E, S> FromIterator<E> for BitMap<E, S>
where
    E: Copy + Into<u32> + TryFrom<u32>,
    S: BitStorage,
{
    fn from_iter<I: IntoIterator<Item = E>>(iter: I) -> Self {
        let mut m = Self::new();
        for e in iter { m.insert(e); }
        m
    }
}

impl<E, S: BitStorage> BitOrAssign for BitMap<E, S> {
    fn bitor_assign(&mut self, rhs: Self) {
        self.bits |= rhs.bits;
    }
}

impl<E, S: BitStorage> BitOr for BitMap<E, S> {
    type Output = Self;
    fn bitor(self, rhs: Self) -> Self {
        Self { bits: self.bits | rhs.bits, _phantom: PhantomData }
    }
}

pub struct Iter<E, S> {
    bits: S,
    _phantom: PhantomData<E>,
}

impl<E, S> Iterator for Iter<E, S>
where
    E: TryFrom<u32>,
    S: BitStorage,
{
    type Item = E;
    fn next(&mut self) -> Option<E> {
        if self.bits.is_zero() { return None; }
        let bit = self.bits.trailing_zeros();
        self.bits = self.bits & (self.bits - S::one());
        E::try_from(bit).ok()
    }
}

impl<E, S> ExactSizeIterator for Iter<E, S>
where
    E: TryFrom<u32>,
    S: BitStorage,
{
    fn len(&self) -> usize {
        self.bits.count_ones() as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    type B = BitMap<u8, u128>;

    #[test]
    fn empty() {
        let m = B::new();
        assert!(m.is_empty());
        assert_eq!(m.len(), 0);
        assert_eq!(m.iter().next(), None);
        assert!(!m.contains(0));
        assert_eq!(m.iter().collect::<Vec<u8>>(), Vec::<u8>::new());
    }

    #[test]
    fn insert_contains() {
        let mut m = B::new();
        m.insert(3);
        m.insert(7);
        m.insert(3);
        assert_eq!(m.len(), 2);
        assert!(m.contains(3));
        assert!(m.contains(7));
        assert!(!m.contains(5));
    }

    #[test]
    fn from_iter_iter_ascending() {
        let m: B = [10u8, 3, 7, 0].into_iter().collect();
        assert_eq!(m.len(), 4);
        assert_eq!(m.iter().next(), Some(0));
        assert_eq!(m.iter().collect::<Vec<u8>>(), vec![0, 3, 7, 10]);
    }

    #[test]
    fn bitor_union() {
        let a: B = [1, 5].into_iter().collect();
        let b: B = [5, 9].into_iter().collect();
        let u = a | b;
        assert_eq!(u.iter().collect::<Vec<u8>>(), vec![1, 5, 9]);
    }

    #[test]
    fn small_storage_u8() {
        // Max-8 BitMap: indices 0..8.
        type Small = BitMap<u8, u8>;
        let mut m = Small::new();
        m.insert(0);
        m.insert(7);
        assert_eq!(m.len(), 2);
        assert_eq!(m.iter().collect::<Vec<u8>>(), vec![0, 7]);
    }

    #[test]
    #[should_panic(expected = "BitMap index")]
    fn small_storage_out_of_range_panics() {
        let mut m: BitMap<u8, u8> = BitMap::new();
        m.insert(8);  // u8 only has bits 0..8
    }
}
