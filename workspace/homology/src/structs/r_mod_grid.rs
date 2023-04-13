use std::collections::HashMap;
use std::ops::{Index};
use yui_core::{Ring, RingOps};
use crate::{GridIdx, GridItr, RModStr, Grid};

pub struct GenericRModGrid<S, I>
where 
    S: RModStr,
    S::R: Ring, for<'x> &'x S::R: RingOps<S::R>,
    I: GridItr,
    I::Item: GridIdx
{
    range: I,
    grid: HashMap<I::Item, S>,
    zero: S
}

impl<S, I> GenericRModGrid<S, I>
where 
    S: RModStr,
    S::R: Ring, for<'x> &'x S::R: RingOps<S::R>,
    I: GridItr,
    I::Item: GridIdx
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
    I: GridItr,
    I::Item: GridIdx
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

impl<S, I> Grid for GenericRModGrid<S, I>
where 
    S: RModStr,
    S::R: Ring, for<'x> &'x S::R: RingOps<S::R>,
    I: GridItr,
    I::Item: GridIdx
{
    type Idx = I::Item;
    type IdxIter = I;

    fn contains_idx(&self, k: Self::Idx) -> bool {
        self.grid.contains_key(&k)
    }

    fn indices(&self) -> Self::IdxIter {
        self.range.clone()
    }
}