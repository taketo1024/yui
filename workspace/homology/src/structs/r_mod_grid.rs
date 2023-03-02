use std::collections::HashMap;
use std::ops::{Index};
use yui_core::{Ring, RingOps};
use crate::{AdditiveIndex, AdditiveIndexRange, RModStr, RModGrid};

pub struct GenericRModGrid<S, I>
where 
    S: RModStr,
    S::R: Ring, for<'x> &'x S::R: RingOps<S::R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    range: I,
    grid: HashMap<I::Item, S>,
    zero: S
}

impl<S, I> GenericRModGrid<S, I>
where 
    S: RModStr,
    S::R: Ring, for<'x> &'x S::R: RingOps<S::R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    pub fn new<F>(range: I, mut f: F) -> Self
    where F: FnMut(I::Item) -> Option<S> {
        let grid = range.clone().filter_map(|i| {
            if let Some(s) = f(i) { 
                Some((i, s))
            } else { 
                None
            }
        }).collect();
        let zero = S::zero();
        Self { range, grid, zero }
    }
}

impl<S, I> Index<I::Item> for GenericRModGrid<S, I>
where
    S: RModStr,
    S::R: Ring, for<'x> &'x S::R: RingOps<S::R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    type Output = S;

    fn index(&self, index: I::Item) -> &Self::Output {
        if let Some(s) = self.grid.get(&index) { 
            s
        } else { 
            &self.zero
        }
    }
}

impl<S, I> RModGrid for GenericRModGrid<S, I>
where 
    S: RModStr,
    S::R: Ring, for<'x> &'x S::R: RingOps<S::R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    type R = S::R;
    type Index = I::Item;
    type IndexRange = I;

    fn in_range(&self, k: Self::Index) -> bool {
        self.grid.contains_key(&k)
    }

    fn range(&self) -> Self::IndexRange {
        self.range.clone()
    }
}