//! [`U256`]: a 256-bit word, used as a wide [`BitStorage`] for [`BitMap`](super::bitmap::BitMap).

use std::ops::{BitAnd, BitOr, BitOrAssign, Shl, Sub};

use super::bitmap::BitStorage;

/// 256-bit storage: two little-endian `u128` limbs `[low, high]`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug, Default)]
pub struct U256([u128; 2]);

impl BitAnd for U256 {
    type Output = Self;
    fn bitand(self, rhs: Self) -> Self {
        Self([self.0[0] & rhs.0[0], self.0[1] & rhs.0[1]])
    }
}

impl BitOr for U256 {
    type Output = Self;
    fn bitor(self, rhs: Self) -> Self {
        Self([self.0[0] | rhs.0[0], self.0[1] | rhs.0[1]])
    }
}

impl BitOrAssign for U256 {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0[0] |= rhs.0[0];
        self.0[1] |= rhs.0[1];
    }
}

impl Shl<u32> for U256 {
    type Output = Self;
    fn shl(self, rhs: u32) -> Self {
        let [lo, hi] = self.0;
        if rhs == 0 {
            self
        } else if rhs < 128 {
            Self([lo << rhs, (hi << rhs) | (lo >> (128 - rhs))])
        } else {
            Self([0, lo << (rhs - 128)])
        }
    }
}

impl Sub for U256 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let (lo, borrow) = self.0[0].overflowing_sub(rhs.0[0]);
        let hi = self.0[1].wrapping_sub(rhs.0[1]).wrapping_sub(borrow as u128);
        Self([lo, hi])
    }
}

impl BitStorage for U256 {
    const WIDTH: u32 = 256;

    fn one() -> Self {
        Self([1, 0])
    }

    fn is_zero(self) -> bool {
        self.0 == [0, 0]
    }

    fn count_ones(self) -> u32 {
        self.0[0].count_ones() + self.0[1].count_ones()
    }

    fn trailing_zeros(self) -> u32 {
        if self.0[0] != 0 {
            self.0[0].trailing_zeros()
        } else if self.0[1] != 0 {
            128 + self.0[1].trailing_zeros()
        } else {
            256
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::bitmap::BitMap;

    #[test]
    fn shl_across_limbs() {
        let one = U256::one();
        let x = one << 200;
        assert_eq!(x.trailing_zeros(), 200);
        assert_eq!(x.count_ones(), 1);

        let y = one << 127;
        assert_eq!((y << 1).trailing_zeros(), 128);
    }

    #[test]
    fn sub_borrows() {
        let x = U256::one() << 128;   // lowest bit of the high limb
        let y = x - U256::one();      // all 128 low bits set
        assert_eq!(y.count_ones(), 128);
        assert_eq!(y.trailing_zeros(), 0);
    }

    #[test]
    fn bitmap_over_128() {
        let mut m: BitMap<u16, U256> = BitMap::new();
        m.insert(3);
        m.insert(130);
        m.insert(255);
        assert_eq!(m.len(), 3);
        assert!(m.contains(130));
        assert!(!m.contains(129));
        assert_eq!(m.iter().collect::<Vec<u16>>(), vec![3, 130, 255]);
    }
}
