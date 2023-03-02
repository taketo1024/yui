use std::{ops::Index};
use yui_core::{Ring, RingOps};
use crate::{GridIdx, GridItr, RModStr};

pub trait RModGrid: Index<Self::Idx>
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>,
    Self::Idx: GridIdx,
    Self::IdxIter: GridItr<Item = Self::Idx>
{
    type R;
    type Idx;
    type IdxIter;

    fn contains_idx(&self, k: Self::Idx) -> bool;
    fn indices(&self) -> Self::IdxIter;

    fn is_free(&self) -> bool { 
        self.indices().all(|i| self[i].is_free())
    }

    fn is_zero(&self) -> bool { 
        self.indices().all(|i| self[i].is_zero())
    }
}