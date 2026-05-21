use derive_more::{Display, Debug};
use num_traits::Zero;
use yui_core::{GetSign, Sign};

#[derive(Clone, Copy, PartialEq, Eq, Display, Debug)]
#[display("{}", _0)]
#[  debug("{}", _0)]
pub struct BraidGen(i32);

impl BraidGen {
    pub fn new(index: usize, sign: Sign) -> Self {
        assert!(!index.is_zero());
        if sign.is_positive() {
            Self(index as i32 )
        } else {
            Self(-(index as i32))
        }
    }

    pub fn index(&self) -> usize {
        self.0.unsigned_abs() as usize
    }

    pub fn sign(&self) -> Sign {
        self.0.sign()
    }

    pub fn inv(&self) -> Self {
        Self(-self.0)
    }
}

impl From<i32> for BraidGen {
    fn from(value: i32) -> Self {
        assert!(!value.is_zero());
        Self(value)
    }
}

pub(super) fn from_raw(value: i32) -> BraidGen {
    BraidGen(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_pos() {
        let g = BraidGen::new(1, Sign::Pos);
        assert_eq!(g.index(), 1);
        assert!(g.sign().is_positive());
    }

    #[test]
    fn new_neg() {
        let g = BraidGen::new(2, Sign::Neg);
        assert_eq!(g.index(), 2);
        assert!(g.sign().is_negative());
    }

    #[test]
    #[should_panic]
    fn new_zero_panics() {
        BraidGen::new(0, Sign::Pos);
    }

    #[test]
    fn from_i32() {
        let g = BraidGen::from(3);
        assert_eq!(g.index(), 3);
        assert!(g.sign().is_positive());

        let g = BraidGen::from(-3);
        assert_eq!(g.index(), 3);
        assert!(g.sign().is_negative());
    }

    #[test]
    #[should_panic]
    fn from_zero_panics() {
        let _ = BraidGen::from(0);
    }

    #[test]
    fn inv() {
        let g = BraidGen::from(2);
        let h = g.inv();
        assert_eq!(h.index(), 2);
        assert!(h.sign().is_negative());
        assert_eq!(h.inv(), g);
    }

    #[test]
    fn to_string() {
        assert_eq!(BraidGen::from(3).to_string(), "3");
        assert_eq!(BraidGen::from(-3).to_string(), "-3");
    }
}
