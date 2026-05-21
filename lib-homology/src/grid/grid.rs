use std::fmt::Display;
use std::ops::{Index, RangeInclusive};

use ahash::AHashMap;
use itertools::Itertools;
use crate::utils::{ToSeqString, ToTableString};
use crate::{GridDeg, isize2, usize2, isize3, usize3};

pub type Grid1<E> = Grid<isize,  E>;
pub type Grid2<E> = Grid<isize2, E>;
pub type Grid3<E> = Grid<isize3, E>;

pub type GridIter<'a, I, V> = std::collections::hash_map::Keys<'a, I, V>;

#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Grid<I, E>
where I: GridDeg {
    data: AHashMap<I, E>,
    #[cfg_attr(feature = "serde", serde(skip))]
    default: E,
}

impl<I, E> Grid<I, E>
where I: GridDeg {
    fn new(data: impl IntoIterator<Item = (I, E)>, default: E) -> Self {
        Self { data: data.into_iter().collect(), default }
    }

    pub fn generate<It, F>(support: It, mut e_map: F) -> Self
    where
        It: IntoIterator<Item = I>,
        F: FnMut(I) -> E,
        E: Default,
    {
        Self::new(support.into_iter().map(|i| (i, e_map(i))), E::default())
    }

    pub fn generate_filtered<It, F>(support: It, mut e_map: F) -> Self
    where
        It: IntoIterator<Item = I>,
        F: FnMut(I) -> Option<E>,
        E: Default,
    {
        Self::new(support.into_iter().filter_map(|i| e_map(i).map(|e| (i, e))), E::default())
    }

    pub fn support(&self) -> GridIter<'_, I, E> {
        self.data.keys()
    }

    pub fn is_supported(&self, i: I) -> bool {
        self.data.contains_key(&i)
    }

    pub fn get(&self, i: I) -> Option<&E> {
        self.data.get(&i)
    }

    pub fn get_default(&self) -> &E {
        &self.default
    }

    pub fn iter(&self) -> impl Iterator<Item = (I, &E)> {
        self.data.iter().map(|(&i, e)| (i, e))
    }
}

impl<E> Grid1<E> {
    pub fn truncated(&self, range: RangeInclusive<isize>) -> Self
    where E: Clone {
        Self::new(
            self.data.iter().filter_map(|(&i, e)| range.contains(&i).then_some((i, e.clone()))),
            self.default.clone(),
        )
    }
}

impl<I, E> Default for Grid<I, E>
where I: GridDeg, E: Default {
    fn default() -> Self {
        Self::new(AHashMap::default(), E::default())
    }
}

impl<I, E> IntoIterator for Grid<I, E>
where I: GridDeg {
    type Item = (I, E);
    type IntoIter = std::collections::hash_map::IntoIter<I, E>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

impl<I, E> Index<I> for Grid<I, E>
where I: GridDeg {
    type Output = E;
    fn index(&self, i: I) -> &Self::Output {
        self.data.get(&i).unwrap_or(&self.default)
    }
}

macro_rules! impl_index {
    ($t:ident, $t2:ident, $t3:ident) => {
        impl<E> Index<($t, $t)> for Grid<$t2, E> {
            type Output = E;
            fn index(&self, i: ($t, $t)) -> &Self::Output {
                &self[$t2::from(i)]
            }
        }

        impl<E> Index<($t, $t, $t)> for Grid<$t3, E> {
            type Output = E;
            fn index(&self, i: ($t, $t, $t)) -> &Self::Output {
                &self[$t3::from(i)]
            }
        }
    };
}

impl_index!(isize, isize2, isize3);
impl_index!(usize, usize2, usize3);

impl<E: Display> ToSeqString<isize> for Grid<isize, E> {
    fn label(&self) -> String {
        "i".to_string()
    }

    fn indices(&self) -> Vec<isize> {
        self.support().sorted().cloned().collect()
    }

    fn entry_at(&self, i: &isize) -> String {
        self.get(*i).map(|e| e.to_string()).unwrap_or_else(|| ".".to_string())
    }
}

impl<E: Display> ToTableString<isize> for Grid<isize2, E> {
    fn labels(&self) -> (String, String) {
        ("i".to_string(), "j".to_string())
    }

    fn indices(&self) -> (Vec<isize>, Vec<isize>) {
        let is = self.support().map(|&isize2(i, _)| i).unique().sorted().collect();
        let js = self.support().map(|&isize2(_, j)| j).unique().sorted().collect();
        (is, js)
    }

    fn entry_at(&self, i: &isize, j: &isize) -> String {
        self.get(isize2(*i, *j)).map(|e| e.to_string()).unwrap_or_else(|| ".".to_string())
    }
}

#[cfg(feature = "tex")]
mod tex_impl {
    use super::*;
    use yui_core::{TeX, tex_table};
    use crate::utils::tex::TeXTable;

    macro_rules! impl_tex_table {
        ($t:ident) => {
            impl<E> TeXTable<$t> for Grid<$t, E>
            where E: TeX {
                fn tex_table(&self, caption: &str, head: &str) -> String {
                    let cols = self.support().map(|&$t(i, _)| i).unique().sorted();
                    let rows = self.support().map(|&$t(_, j)| j).unique().sorted().rev();

                    tex_table(caption, head, rows, cols, |&j, &i| {
                        self.get($t(i, j))
                            .map(|e| e.tex_string())
                            .unwrap_or_else(|| ".".to_string())
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
        assert_eq!(g.get( 1), Some(&10));
        assert_eq!(g.get(-1), None);
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
        assert_eq!(g.get(isize2(1, 2)), Some(&12));
        assert_eq!(g.get(isize2(3, 3)), None);
    }
}