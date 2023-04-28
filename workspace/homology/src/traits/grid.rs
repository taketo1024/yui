use std::fmt::Display;
use std::hash::Hash;
use std::ops::{Add, Neg, Sub, RangeInclusive};

pub trait GridIdx:
    Clone
    + Copy
    + PartialEq
    + Eq
    + Hash
    + Display
    + Add<Output = Self>
    + Sub<Output = Self>
    + Neg<Output = Self>
{}

impl<T> GridIdx for T where
    T: Clone
        + Copy
        + PartialEq
        + Eq
        + Hash
        + Display
        + Add<Output = Self>
        + Sub<Output = Self>
        + Neg<Output = Self>
{}

pub trait GridItr: Iterator + Clone
where
    Self::Item: GridIdx,
{}

impl<T> GridItr for T
where
    T: Iterator + Clone,
    T::Item: GridIdx,
{}

pub trait Grid: Sized
where 
    Self::Idx: GridIdx,
    Self::IdxIter: GridItr<Item = Self::Idx>,
{
    type Idx;
    type IdxIter;
    type Output;

    fn indices(&self) -> Self::IdxIter;
    fn get(&self, i: Self::Idx) -> Option<&Self::Output>;

    fn contains_idx(&self, i: Self::Idx) -> bool {
        self.get(i).is_some()
    }

    fn iter(&self) -> GridIter<'_, Self> { 
        GridIter { grid: self, iter: self.indices() }
    }
}

pub struct GridIter<'a, G>
where G: Grid { 
    grid: &'a G,
    iter: G::IdxIter,
}

impl<'a, G> Iterator for GridIter<'a, G>
where G: Grid {
    type Item = (G::Idx, &'a G::Output);

    fn next(&mut self) -> Option<Self::Item> {
        let Some(i) = self.iter.next() else { return None };
        let v = self.grid.get(i).unwrap();
        Some((i, v))
    }
}

pub trait Shift { 
    fn shift(self, i: isize) -> Self;
}

impl Shift for RangeInclusive<isize> { 
    fn shift(self, i: isize) -> Self {
        RangeInclusive::new(
            self.start() + i, 
            self.end() + i
        )
    }
}

