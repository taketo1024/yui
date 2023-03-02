use std::{ops::Index};
use yui_core::{Ring, RingOps};
use crate::{AdditiveIndex, AdditiveIndexRange, RModStr};

pub trait RModGrid: Index<Self::Index>
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>,
    Self::Index: AdditiveIndex,
    Self::IndexRange: AdditiveIndexRange<Item = Self::Index>
{
    type R;
    type Index;
    type IndexRange;

    fn in_range(&self, k: Self::Index) -> bool;
    fn range(&self) -> Self::IndexRange;

    fn is_free(&self) -> bool { 
        self.range().all(|i| self[i].is_free())
    }

    fn is_zero(&self) -> bool { 
        self.range().all(|i| self[i].is_zero())
    }
}