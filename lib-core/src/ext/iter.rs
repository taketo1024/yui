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
