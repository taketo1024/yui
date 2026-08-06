use derive_more::{Display, Debug};
use num_traits::Zero;
use yui_core::num::{GetSign, Sign};

#[derive(Clone, Copy, PartialEq, Eq, Display, Debug)]
#[display("{}", _0)]
#[  debug("{}", _0)]
pub struct BraidGen(i8);

impl BraidGen {
    pub fn new(val: i8) -> Self {
        Self::from(val)
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

macro_rules! impl_from_int {
    ($($t:ty),* $(,)?) => {
        $(
            impl From<$t> for BraidGen {
                fn from(value: $t) -> Self {
                    assert!(!value.is_zero());
                    let v = i8::try_from(value).expect("BraidGen value must fit in i8");
                    Self(v)
                }
            }
        )*
    };
}

impl_from_int!(i8, i16, i32, i64);

pub(super) fn from_raw(value: i8) -> BraidGen {
    BraidGen(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_pos() {
        let g = BraidGen::new(1);
        assert_eq!(g.index(), 1);
        assert!(g.sign().is_positive());
    }

    #[test]
    fn new_neg() {
        let g = BraidGen::new(-2);
        assert_eq!(g.index(), 2);
        assert!(g.sign().is_negative());
    }

    #[test]
    #[should_panic]
    fn new_zero_panics() {
        BraidGen::new(0);
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
