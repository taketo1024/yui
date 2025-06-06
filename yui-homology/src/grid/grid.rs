use std::ops::{Index, RangeInclusive};
use std::fmt::Display;

use ahash::AHashMap;
use itertools::Itertools;
use crate::{GridDeg, isize2, usize2, isize3, usize3};

pub trait GridTrait<I>
where I: GridDeg { 
    type Support: Iterator<Item = I>;
    type Item;

    fn support(&self) -> Self::Support;
    fn is_supported(&self, i: I) -> bool;
    fn get(&self, i: I) -> &Self::Item;
    fn get_default(&self) -> &Self::Item;
}

pub type Grid1<E> = Grid<isize,  E>;
pub type Grid2<E> = Grid<isize2, E>;
pub type Grid3<E> = Grid<isize3, E>;

pub type GridIter<I> = std::vec::IntoIter<I>;

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

impl<E> Grid1<E> {
    pub fn truncated(&self, range: RangeInclusive<isize>) -> Self
    where E: Clone { 
        let support = self.support().filter(|i| range.contains(i));
        Self::generate_with_default(support, |i| self[i].clone(), self.default.clone())
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
    type Support = GridIter<I>;
    type Item = E;

    fn support(&self) -> Self::Support {
        self.support.clone().into_iter()
    }

    fn is_supported(&self, i: I) -> bool {
        self.data.contains_key(&i)
    }

    fn get(&self, i: I) -> &E {
        self.data.get(&i).unwrap_or(&self.default)
    }

    fn get_default(&self) -> &E { 
        &self.default
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
impl_index!(usize, usize2, usize3);

pub trait DisplaySeq<I> {
    fn display_seq(&self, label: &str) -> String;
    fn print_seq(&self, label: &str) {
        println!("{}", self.display_seq(label))
    }
}

macro_rules! impl_print_seq {
    ($t:ident) => {
        impl<G> DisplaySeq<$t> for G
        where G: GridTrait<$t>, G::Item: Display {
            fn display_seq(&self, label: &str) -> String { 
                use yui::util::format::table;

                table(label, [""].iter(), self.support(), |_, &i| {
                    self.get(i).to_string()
                })
            }
        }                
    };
}

impl_print_seq!(isize);
impl_print_seq!(usize);

pub trait DisplayTable<I> {
    fn display_table(&self, label0: &str, label1: &str) -> String;
    fn print_table(&self, label0: &str, label1: &str) {
        println!("{}", self.display_table(label0, label1))
    }
}

macro_rules! impl_print_table {
    ($t:ident) => {
        impl<G> DisplayTable<$t> for G
        where G: GridTrait<$t>, G::Item: Display {
            fn display_table(&self, label0: &str, label1: &str) -> String {
                use yui::util::format::table;

                let def_str = self.get_default().to_string();
                let head = format!("{}\\{}", label1, label0);
                let cols = self.support().map(|$t(i, _)| i).unique().sorted();
                let rows = self.support().map(|$t(_, j)| j).unique().sorted().rev();
        
                table(head, rows, cols, |&j, &i| {
                    let str = self.get($t(i, j)).to_string();
                    if str == def_str { 
                        ".".to_string()
                    } else { 
                        str
                    }
                })
            }
        }
    };
}

impl_print_table!(isize2);
impl_print_table!(usize2);

#[cfg(feature = "tex")]
pub mod tex { 
    use super::*;
    use yui::tex::{TeX, tex_table};

    pub trait TeXTable<I> {
        fn tex_table(&self, caption: &str, head: &str) -> String;
    }

    macro_rules! impl_tex_table {
        ($t:ident) => {
            impl<G> TeXTable<$t> for G
            where G: GridTrait<$t>, G::Item: TeX {
                fn tex_table(&self, caption: &str, head: &str) -> String {
                    let def_str = self.get_default().tex_string();
                    let cols = self.support().map(|$t(i, _)| i).unique().sorted();
                    let rows = self.support().map(|$t(_, j)| j).unique().sorted().rev();
            
                    tex_table(caption, head, rows, cols, |&j, &i| {
                        let str = self.get($t(i, j)).tex_string();
                        if str == def_str { 
                            ".".to_string()
                        } else { 
                            str
                        }
                    }, true, false)
                }
            }
        };
    }
    
    impl_tex_table!(isize2);
}

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
        let g = Grid2::generate(
            cartesian!(0..=3, 0..=2).map(|(i, j)| isize2(i, j)), 
            |i| i.0 * 10 + i.1
        );

        assert!( g.is_supported(isize2(1, 2)));
        assert!(!g.is_supported(isize2(3, 3)));
        assert_eq!(g.get(isize2(1, 2)), &12);
        assert_eq!(g.get(isize2(3, 3)), &0);

        let _table = g.display_table("i", "j");
        // println!("{_table}");
    }

    #[test]
    fn table_rmod() { 
        use yui::FF2;
        use crate::GenericSummand;
        use cartesian::cartesian;
        
        let g = Grid2::generate(
            cartesian!(0..=3, 0..=2).map(|(i, j)| isize2(i, j)), 
            |i| GenericSummand::<_, FF2>::generate(i.0, (i.0 * 10 + i.1) as usize, vec![], None)
        );
        let _table = g.display_table("i", "j");
        // println!("{_table}");
    }

    #[cfg(feature = "tex")]
    #[test]
    fn tex_table() { 
        use super::tex::TeXTable;
        use cartesian::cartesian;
        let g = Grid2::generate(
            cartesian!(0..=3, 0..=2).map(|(i, j)| isize2(i, j)), 
            |i| (i.0 * 10 + i.1) as i32
        );
        let _table = g.tex_table("Caption", "i, j");
        // println!("{_table}");
    }

    #[cfg(feature = "tex")]
    #[test]
    fn tex_table_rmod() { 
        use super::tex::TeXTable;
        use crate::GenericSummand;
        use cartesian::cartesian;
        
        let g = Grid2::generate(
            cartesian!(0..=3, 0..=2).map(|(i, j)| isize2(i, j)), 
            |i| GenericSummand::<_, i32>::generate(i.0, (i.0 * 10 + i.1) as usize, vec![], None)
        );
        let _table = g.tex_table("Caption", "i, j");
        // println!("{_table}");
    }
}