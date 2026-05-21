use std::fmt::Display;
use std::hash::Hash;
use std::ops::Index;

use ahash::AHashMap;
use delegate::delegate;
use itertools::Itertools;

use crate::utils::{ToSeqString, ToTableString};
use crate::{isize2, isize3};

pub type Grid1<V> = Grid<isize,  V>;
pub type Grid2<V> = Grid<isize2, V>;
pub type Grid3<V> = Grid<isize3, V>;

/// A sparse `K`-indexed table of values, backed by a `HashMap` plus a stored
/// `V::default()`. `Index<K>` returns a reference to the stored value if the
/// key is present, otherwise to the default — so `grid[k]` never panics.
#[derive(Clone)]
pub struct Grid<K, V>
where K: Hash + Eq + Copy, V: Default {
    data: AHashMap<K, V>,
    default: V,
}

// Manual Default impl: derive(Default) would also require `K: Default`,
// which is unnecessary since `AHashMap::default()` doesn't need it.
impl<K, V> Default for Grid<K, V>
where K: Hash + Eq + Copy, V: Default {
    fn default() -> Self {
        Self { data: AHashMap::default(), default: V::default() }
    }
}

impl<K, V> Grid<K, V>
where K: Hash + Eq + Copy, V: Default {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn get_default(&self) -> &V {
        &self.default
    }

    delegate! {
        to self.data {
            #[expr(self.data.get(&k))]
            pub fn get(&self, k: K) -> Option<&V>;
            #[expr(self.data.contains_key(&k))]
            pub fn contains_key(&self, k: K) -> bool;

            pub fn keys(&self) -> impl Iterator<Item = &K> + '_;
            pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> + '_;
            pub fn insert(&mut self, k: K, v: V) -> Option<V>;
            pub fn len(&self) -> usize;
            pub fn is_empty(&self) -> bool;
        }
    }
}

impl<K, V> Index<K> for Grid<K, V>
where K: Hash + Eq + Copy, V: Default {
    type Output = V;
    fn index(&self, k: K) -> &V {
        self.data.get(&k).unwrap_or(&self.default)
    }
}

impl<V> Index<(isize, isize)> for Grid<isize2, V>
where V: Default {
    type Output = V;
    fn index(&self, k: (isize, isize)) -> &V {
        &self[isize2::from(k)]
    }
}

impl<V> Index<(isize, isize, isize)> for Grid<isize3, V>
where V: Default {
    type Output = V;
    fn index(&self, k: (isize, isize, isize)) -> &V {
        &self[isize3::from(k)]
    }
}

impl<K, V> FromIterator<(K, V)> for Grid<K, V>
where K: Hash + Eq + Copy, V: Default {
    fn from_iter<It: IntoIterator<Item = (K, V)>>(iter: It) -> Self {
        Self { data: iter.into_iter().collect(), default: V::default() }
    }
}

impl<K, V> IntoIterator for Grid<K, V>
where K: Hash + Eq + Copy, V: Default {
    type Item = (K, V);
    type IntoIter = std::collections::hash_map::IntoIter<K, V>;
    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

impl<V: Display + Default> ToSeqString<isize> for Grid<isize, V> {
    fn label(&self) -> String {
        "i".to_string()
    }

    fn indices(&self) -> Vec<isize> {
        self.keys().copied().sorted().collect()
    }

    fn entry_at(&self, i: &isize) -> String {
        self.get(*i).map(|v| v.to_string()).unwrap_or_else(|| ".".to_string())
    }
}

impl<V: Display + Default> ToTableString<isize> for Grid<isize2, V> {
    fn labels(&self) -> (String, String) {
        ("i".to_string(), "j".to_string())
    }

    fn indices(&self) -> (Vec<isize>, Vec<isize>) {
        let is = self.keys().map(|&isize2(i, _)| i).unique().sorted().collect();
        let js = self.keys().map(|&isize2(_, j)| j).unique().sorted().collect();
        (is, js)
    }

    fn entry_at(&self, i: &isize, j: &isize) -> String {
        self.get(isize2(*i, *j)).map(|v| v.to_string()).unwrap_or_else(|| ".".to_string())
    }
}

#[cfg(feature = "tex")]
impl<V> crate::utils::tex::TeXTable<isize2> for Grid<isize2, V>
where V: yui_core::TeX + Default {
    fn tex_table(&self, caption: &str, head: &str) -> String {
        let cols = self.keys().map(|&isize2(i, _)| i).unique().sorted();
        let rows = self.keys().map(|&isize2(_, j)| j).unique().sorted().rev();

        yui_core::tex_table(caption, head, rows, cols, |&j, &i| {
            self.get(isize2(i, j))
                .map(|e| e.tex_string())
                .unwrap_or_else(|| ".".to_string())
        }, true, false)
    }
}
