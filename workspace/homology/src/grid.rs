use std::collections::HashMap;
use std::ops::Index;
use std::fmt::Display;

use itertools::Itertools;
use yui_core::{isize2, usize2, Deg, isize3};

pub trait GridTrait<I>
where I: Deg { 
    type Itr: Iterator<Item = I>;
    type E;

    fn support(&self) -> Self::Itr;
    fn is_supported(&self, i: I) -> bool;
    fn get(&self, i: I) -> &Self::E;
}

pub type Grid<E>  = GridBase<isize,  E>;
pub type Grid2<E> = GridBase<isize2, E>;
pub type Grid3<E> = GridBase<isize3, E>;

pub type GridIter<I> = std::vec::IntoIter<I>;

pub struct GridBase<I, E>
where I: Deg { 
    support: Vec<I>,
    data: HashMap<I, E>,
    dflt: E
}

impl<I, E> GridBase<I, E>
where I: Deg, E: Default { 
    pub fn new<It, F>(support: It, mut e_map: F) -> Self
    where 
        It: Iterator<Item = I>, 
        F: FnMut(I) -> E
    {
        let support = support.collect_vec();
        let data = support.iter().map(|&i| (i, e_map(i))).collect();
        Self::new_raw(support, data)
    }

    pub fn new_raw(support: Vec<I>, data: HashMap<I, E>) -> Self { 
        Self { support, data, dflt: E::default() }
    }

    pub fn insert(&mut self, i: I, e: E) {
        self.data.insert(i, e);
    }
}

impl<I, E> GridTrait<I> for GridBase<I, E>
where I: Deg { 
    type Itr = GridIter<I>;
    type E = E;

    fn support(&self) -> Self::Itr {
        self.support.clone().into_iter()
    }

    fn is_supported(&self, i: I) -> bool {
        self.data.contains_key(&i)
    }

    fn get(&self, i: I) -> &E {
        self.data.get(&i).unwrap_or(&self.dflt)
    }
}

impl<I, E> Index<I> for GridBase<I, E>
where I: Deg {
    type Output = E;
    fn index(&self, i: I) -> &Self::Output {
        self.get(i)
    }
}

pub trait DisplayForGrid {
    fn display_for_grid(&self) -> String;
}

impl<T> DisplayForGrid for T
where T: Display {
    fn display_for_grid(&self) -> String {
        self.to_string()
    }
}

pub trait DisplaySeq<I> {
    fn display_seq(&self) -> String;
    fn print_seq(&self) {
        println!("{}", self.display_seq())
    }
}

macro_rules! impl_print_seq {
    ($t:ident) => {
        impl<T> DisplaySeq<$t> for T
        where T: GridTrait<$t>, T::E: DisplayForGrid {
            fn display_seq(&self) -> String {
                use yui_utils::table;
                let str = table("i", [""].iter(), self.support(), |_, &i| {
                    self.get(i).display_for_grid()
                });
                str
            }
        }                
    };
}

impl_print_seq!(isize);
impl_print_seq!(usize);

pub trait DisplayTable<I> {
    fn display_table(&self) -> String;
    fn print_table(&self) {
        println!("{}", self.display_table())
    }
}

macro_rules! impl_print_table {
    ($t:ident) => {
        impl<T> DisplayTable<$t> for T
        where T: GridTrait<$t>, T::E: DisplayForGrid {
            fn display_table(&self) -> String {
                use yui_utils::table;
        
                let cols = self.support().map(|$t(i, _)| i).unique().sorted();
                let rows = self.support().map(|$t(_, j)| j).unique().sorted().rev();
        
                let str = table("j\\i", rows, cols, |&j, &i| {
                    self.get($t(i, j)).display_for_grid()
                });
        
                str
            }
        }
    };
}

impl_print_table!(isize2);
impl_print_table!(usize2);

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn grid() { 
        let g = Grid::new(0..=3, |i| i * 10);

        assert_eq!(g.is_supported( 1), true);
        assert_eq!(g.is_supported(-1), false);
        assert_eq!(g.get( 1), &10);
        assert_eq!(g.get(-1), &0); // default

        let _seq = g.display_seq();
        // println!("{_seq}");
    }

    #[test]
    fn grid2() { 
        use cartesian::cartesian;
        let g = Grid2::new(cartesian!(0..=3, 0..=2).map(|(i, j)| isize2(i, j)), |i| i.0 * 10 + i.1);

        assert_eq!(g.is_supported(isize2(1, 2)), true);
        assert_eq!(g.is_supported(isize2(3, 3)), false);
        assert_eq!(g.get(isize2(1, 2)), &12);
        assert_eq!(g.get(isize2(3, 3)), &0);

        let _table = g.display_table();
        // println!("{_table}");
    }
}