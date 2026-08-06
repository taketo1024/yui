//! [`RangeExt`]: shift the endpoints of a range, and [`empty_range`].

use std::ops::{Add, Neg, Range, RangeInclusive, Sub};

/// The empty `RangeInclusive<isize>`, for a grading support with no elements.
/// `start > end`, so it yields nothing and `contains` is always false.
#[allow(clippy::reversed_empty_ranges)] // deliberately reversed: that is what makes it empty.
pub fn empty_range() -> RangeInclusive<isize> {
    0 ..= -1
}

/// Shift the endpoints of a [`Range`] or [`RangeInclusive`] by independent
/// left/right offsets.
pub trait RangeExt
where Self::Idx: Copy, Self: Sized {
    type Idx; 
    fn mv(&self, l: Self::Idx, r: Self::Idx) -> Self;
    fn shift(&self, a: Self::Idx) -> Self { 
        self.mv(a, a)
    }
}

impl<Idx> RangeExt for Range<Idx>
where Idx: Copy + Add<Output = Idx> + Sub<Output = Idx> + Neg<Output = Idx> {
    type Idx = Idx;

    fn mv(&self, l: Self::Idx, r: Self::Idx) -> Self {
        (self.start + l) .. (self.end + r)
    }
}

impl<Idx> RangeExt for RangeInclusive<Idx>
where Idx: Copy + Add<Output = Idx> + Sub<Output = Idx> + Neg<Output = Idx> {
    type Idx = Idx;
    
    fn mv(&self, l: Self::Idx, r: Self::Idx) -> Self {
        (*self.start() + l) ..= (*self.end() + r)
    }
}

#[cfg(test)]
mod tests {
    use crate::ext::{empty_range, RangeExt};

    #[test]
    fn empty() {
        let r = empty_range();
        assert!(r.is_empty());
        assert_eq!(r.clone().count(), 0);
        assert!(!r.contains(&0));
    }

    #[test]
    fn range() {
        let r = -1 .. 3;
        assert_eq!(r.mv(2, 3), 1 .. 6);
        assert_eq!(r.shift(2), 1 .. 5);
    }

    #[test]
    fn range_incl() { 
        let r = -1 ..= 3;
        assert_eq!(r.mv(2, 3), 1 ..= 6);
        assert_eq!(r.shift(2), 1 ..= 5);
    }
}