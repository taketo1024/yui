use std::ops::Index;
use std::fmt::Display;

use ahash::AHashMap;
use itertools::Itertools;
use crate::{GridDeg, isize2, isize3};

pub type Grid1<E> = Grid<isize,  E>;
pub type Grid2<E> = Grid<isize2, E>;
pub type Grid3<E> = Grid<isize3, E>;
pub type GridIter<I> = std::vec::IntoIter<I>;

pub trait GridTrait<I>
where I: GridDeg { 
    type Itr: Iterator<Item = I>;
    type Output;

    fn support(&self) -> Self::Itr;
    fn is_supported(&self, i: I) -> bool;
    fn get(&self, i: I) -> &Self::Output;
}

#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Grid<I, E>
where I: GridDeg { 
    support: Vec<I>,
    data: AHashMap<I, E>,
    #[cfg_attr(feature = "serde", serde(skip))]
    default: E
}

impl<I, E> Grid<I, E>
where I: GridDeg { 
    fn new(support: Vec<I>, data: AHashMap<I, E>, default: E) -> Self { 
        Self { support, data, default }
    }

    pub fn generate<It, F>(support: It, e_map: F) -> Self
    where 
        It: IntoIterator<Item = I>, 
        F: FnMut(I) -> E,
        E: Default
    {
        Self::generate_with_default(support, e_map, E::default())
    }

    pub fn generate_with_default<It, F>(support: It, mut e_map: F, default: E) -> Self
    where 
        It: IntoIterator<Item = I>, 
        F: FnMut(I) -> E
    {
        let support = support.into_iter().collect_vec();
        let data = support.iter().map(|&i| (i, e_map(i))).collect();
        Self::new(support, data, default)
    }

    pub fn get_default(&self) -> &E { 
        &self.default
    }

    pub fn insert(&mut self, i: I, e: E) {
        self.data.insert(i, e);
    }

    pub fn remove(&mut self, i: I) -> Option<E> { 
        self.data.remove(&i)
    }

    pub fn get_mut(&mut self, i: I) -> Option<&mut E> { 
        self.data.get_mut(&i)
    }

    pub fn iter(&self) -> impl Iterator<Item = (I, &E)> {
        self.support.iter().map(|&i| (i, self.get(i)))
    }

    pub fn map<E2, F>(&self, mut f: F) -> Grid<I, E2>
    where F: FnMut(&E) -> E2 
    {
        let d = f(self.get_default());
        Grid::generate_with_default(
            self.support(), 
            |i| f(self.get(i)),
            d
        )
    }
}

impl<I, E> Default for Grid<I, E>
where I: GridDeg, E: Default {
    fn default() -> Self {
        Self::new(Vec::default(), AHashMap::default(), E::default())
    }
}

impl<I, E> IntoIterator for Grid<I, E>
where I: GridDeg {
    type Item = (I, E);
    type IntoIter = std::vec::IntoIter<(I, E)>;

    fn into_iter(mut self) -> Self::IntoIter {
        let support = self.support();
        support.flat_map(|i| { 
            self.data.remove(&i).map(|e| (i, e))
        }).collect_vec().into_iter()
    }
}


impl<I, E> GridTrait<I> for Grid<I, E>
where I: GridDeg { 
    type Itr = GridIter<I>;
    type Output = E;

    fn support(&self) -> Self::Itr {
        self.support.clone().into_iter()
    }

    fn is_supported(&self, i: I) -> bool {
        self.data.contains_key(&i)
    }

    fn get(&self, i: I) -> &E {
        self.data.get(&i).unwrap_or(&self.default)
    }
}

impl<I, E> Index<I> for Grid<I, E>
where I: GridDeg {
    type Output = E;
    fn index(&self, i: I) -> &Self::Output {
        self.get(i)
    }
}

impl<I, E> FromIterator<(I, E)> for Grid<I, E>
where I: GridDeg, E: Default {
    fn from_iter<T: IntoIterator<Item = (I, E)>>(iter: T) -> Self {
        let init = (vec![], AHashMap::new());
        let (support, data) = iter.into_iter().fold(init, |(mut support, mut data), (i, e)| {
            support.push(i);
            data.insert(i, e);
            (support, data)
        });
        Self::new(support, data, E::default())
    }
}

macro_rules! impl_index {
    ($t:ident, $t2:ident, $t3:ident) => {
        impl<E> Index<($t, $t)> for Grid<$t2, E> {
            type Output = E;
            fn index(&self, i: ($t, $t)) -> &Self::Output {
                self.get(i.into())
            }
        }
        
        impl<E> Index<($t, $t, $t)> for Grid<$t3, E> {
            type Output = E;
            fn index(&self, i: ($t, $t, $t)) -> &Self::Output {
                self.get(i.into())
            }
        }
    };
}

impl_index!(isize, isize2, isize3);

pub trait DisplayAt<I>: GridTrait<I>
where I: GridDeg { 
    fn display_at(&self, i: I) -> String;
}

impl<I, E> DisplayAt<I> for Grid<I, E>
where I: GridDeg, E: Display {
    fn display_at(&self, i: I) -> String {
        self.get(i).to_string()
    }
}

pub trait DisplaySeq<I>: DisplayAt<I>
where I: GridDeg {
    fn display_seq(&self, label: &str) -> String {
        use yui::util::format::table;
        table(label, [""].iter(), self.support(), |_, &i| {
            self.display_at(i)
        })
    }

    fn print_seq(&self, label: &str) {
        println!("{}", self.display_seq(label))
    }
}

macro_rules! impl_print_seq {
    ($t:ident) => {
        impl<E> DisplaySeq<$t> for Grid<$t, E> where E: Display {}
    };
}

impl_print_seq!(isize);

pub trait Grid2Trait<I, I2> {
    fn rows(&self) -> Vec<I>;
    fn cols(&self) -> Vec<I>;
}

impl<G> Grid2Trait<isize, isize2> for G
where G: GridTrait<isize2> {
    fn rows(&self) -> Vec<isize> {
        self.support().map(|t| t.1).unique().sorted().collect_vec()
    }

    fn cols(&self) -> Vec<isize> {
        self.support().map(|t| t.0).unique().sorted().collect_vec()
    }
}

pub trait DisplayTable<I, I2>: DisplayAt<I2> + Grid2Trait<I, I2>
where I: GridDeg, I2: GridDeg + From<(I, I)> {
    fn display_table(&self, label0: &str, label1: &str) -> String {
        use yui::util::format::table;

        let head = format!("{}\\{}", label1, label0);
        let rows = self.rows().into_iter().rev();
        let cols = self.cols();

        table(head, rows, cols, |&j, &i|
            self.display_at((i, j).into())
        )
    }

    fn print_table(&self, label0: &str, label1: &str) {
        println!("{}", self.display_table(label0, label1))
    }
}

macro_rules! impl_print_table {
    ($i:ident, $i2:ident) => {
        impl<E> DisplayTable<$i, $i2> for Grid<$i2, E> where E: Display {}
    };
}

impl_print_table!(isize, isize2);

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn grid() { 
        let g = Grid1::generate(0..=3, |i| i * 10);

        assert!( g.is_supported( 1));
        assert!(!g.is_supported(-1));
        assert_eq!(g.get( 1), &10);
        assert_eq!(g.get(-1), &0); // default

        let _seq = g.display_seq("i");
        // println!("{_seq}");
    }

    #[test]
    fn grid2() { 
        use cartesian::cartesian;
        let g = Grid2::generate(cartesian!(0..=3, 0..=2).map(|(i, j)| isize2(i, j)), |i| i.0 * 10 + i.1);

        assert!( g.is_supported(isize2(1, 2)));
        assert!(!g.is_supported(isize2(3, 3)));
        assert_eq!(g.get(isize2(1, 2)), &12);
        assert_eq!(g.get(isize2(3, 3)), &0);

        let _table = g.display_table("i", "j");
        // println!("{_table}");
    }
}