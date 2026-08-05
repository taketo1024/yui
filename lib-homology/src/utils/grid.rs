//! [`Grid<K, V>`]: a sparse `K`-indexed table backed by a hash map plus a
//! stored default value, used as the storage primitive for [`GrMod`].

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

impl<V> crate::utils::tex::ToTexSeq<isize> for Grid<isize, V>
where V: Display + yui_core::TeX + Default {
    fn tex_entry_at(&self, i: &isize) -> String {
        self.get(*i).map(|v| v.tex_string()).unwrap_or_else(|| ".".to_string())
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

impl<V> crate::utils::tex::ToTexTable<isize> for Grid<isize2, V>
where V: Display + yui_core::TeX + Default {
    fn tex_entry_at(&self, i: &isize, j: &isize) -> String {
        self.get(isize2(*i, *j)).map(|e| e.tex_string()).unwrap_or_else(|| ".".to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_falls_back_to_the_default() {
        // the point of `Grid`: a missing key reads as `V::default()` rather than panicking.
        let g = Grid1::<i32>::from_iter([(0, 3), (2, 5)]);

        assert_eq!(g[0], 3);
        assert_eq!(g[2], 5);
        assert_eq!(g[1], 0);
        assert_eq!(g[-100], 0);
        assert_eq!(g.get_default(), &0);

        // reading an absent key must not insert it.
        assert_eq!(g.len(), 2);
        assert!(!g.contains_key(1));
    }

    #[test]
    fn empty_grid() {
        let g = Grid1::<i32>::new();
        assert!(g.is_empty());
        assert_eq!(g.len(), 0);
        assert_eq!(g[0], 0);
        assert_eq!(g.get(0), None);
        assert_eq!(g.keys().count(), 0);
    }

    #[test]
    fn get_distinguishes_absent_from_default_valued() {
        // `index` cannot tell the two apart, `get` must.
        let g = Grid1::<i32>::from_iter([(0, 0)]);
        assert_eq!(g[0], 0);
        assert_eq!(g[1], 0);
        assert_eq!(g.get(0), Some(&0));
        assert_eq!(g.get(1), None);
    }

    #[test]
    fn insert_returns_the_previous_value() {
        let mut g = Grid1::<i32>::new();
        assert_eq!(g.insert(0, 3), None);
        assert_eq!(g.insert(0, 5), Some(3));
        assert_eq!(g[0], 5);
        assert_eq!(g.len(), 1);
    }

    #[test]
    fn tuple_index_matches_the_isize2_key() {
        let g = Grid2::<i32>::from_iter([(isize2(1, -2), 7)]);
        assert_eq!(g[isize2(1, -2)], 7);
        assert_eq!(g[(1, -2)], 7);
        assert_eq!(g[(0, 0)], 0);
    }

    #[test]
    fn tuple_index_matches_the_isize3_key() {
        let g = Grid3::<i32>::from_iter([(isize3(1, -2, 3), 7)]);
        assert_eq!(g[isize3(1, -2, 3)], 7);
        assert_eq!(g[(1, -2, 3)], 7);
        assert_eq!(g[(0, 0, 0)], 0);
    }

    #[test]
    fn iter_round_trips() {
        let g = Grid1::<i32>::from_iter([(0, 3), (2, 5)]);
        let mut es = g.iter().map(|(&k, &v)| (k, v)).collect::<Vec<_>>();
        es.sort();
        assert_eq!(es, vec![(0, 3), (2, 5)]);

        let mut es = g.clone().into_iter().collect::<Vec<_>>();
        es.sort();
        assert_eq!(es, vec![(0, 3), (2, 5)]);
    }

    #[test]
    fn seq_string_lists_only_the_stored_keys() {
        // `indices()` is the key set, so a gap is skipped rather than shown as `.` — index 1
        // does not appear at all.
        let g = Grid1::<i32>::from_iter([(0, 3), (2, 5)]);
        fn cells(s: &str) -> Vec<Vec<String>> {
            s.lines()
                .map(|l| l.split_whitespace().map(str::to_string).collect::<Vec<_>>())
                .filter(|l| !l.is_empty())
                .collect()
        }

        assert_eq!(cells(&g.to_seq_string()), cells("
            i  0  2
               3  5
        "));
    }
}
