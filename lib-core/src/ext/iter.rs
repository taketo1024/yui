use std::ops::RangeInclusive;
use itertools::{Itertools, MinMaxResult};

/// Extension methods on [`Iterator`].
pub trait IteratorExt: Iterator + Sized {
    /// The inclusive range `min..=max` of the items, or `None` if the iterator
    /// is empty.
    fn range(self) -> Option<RangeInclusive<Self::Item>>
    where Self::Item: Ord + Copy {
        match self.minmax() {
            MinMaxResult::NoElements       => None,
            MinMaxResult::OneElement(x)    => Some(x ..= x),
            MinMaxResult::MinMax(min, max) => Some(min ..= max),
        }
    }
}

impl<I: Iterator> IteratorExt for I {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn range_empty() {
        let r: Option<RangeInclusive<i32>> = std::iter::empty().range();
        assert_eq!(r, None);
    }

    #[test]
    fn range_single() {
        assert_eq!([7].into_iter().range(), Some(7 ..= 7));
    }

    #[test]
    fn range_sorted() {
        assert_eq!([1, 2, 3, 4, 5].into_iter().range(), Some(1 ..= 5));
    }

    #[test]
    fn range_unsorted() {
        assert_eq!([3, -2, 5, 1, -7, 4].into_iter().range(), Some(-7 ..= 5));
    }

    #[test]
    fn range_duplicates() {
        assert_eq!([2, 2, 2].into_iter().range(), Some(2 ..= 2));
    }
}
