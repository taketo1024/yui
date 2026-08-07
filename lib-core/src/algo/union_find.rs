//! Disjoint-set structures: [`UnionFind`] over `0..n` (path compression and
//! union by rank, from `petgraph`), and [`KeyedUnionFind`] over hashable keys.

use std::collections::HashMap;
use std::hash::Hash;

use indexmap::IndexSet;
use petgraph::unionfind::UnionFind as PgUnionFind;

/// A disjoint-set data structure over `0..n`. Backed by [`petgraph::unionfind::UnionFind`],
/// which uses path compression and union by rank for `O(α(n))` amortized operations.
#[derive(Clone)]
pub struct UnionFind {
    inner: PgUnionFind<usize>,
}

impl UnionFind {
    pub fn new(n: usize) -> Self {
        Self { inner: PgUnionFind::new(n) }
    }

    pub fn extend(&mut self, l: usize) {
        for _ in 0..l {
            self.inner.new_set();
        }
    }

    pub fn size(&self) -> usize {
        self.inner.len()
    }

    pub fn root(&self, i: usize) -> usize {
        self.inner.find(i)
    }

    pub fn is_same(&self, i: usize, j: usize) -> bool {
        self.inner.equiv(i, j)
    }

    pub fn union(&mut self, i: usize, j: usize) {
        self.inner.union(i, j);
    }

    pub fn into_disjoint(self) -> Vec<Vec<usize>> {
        let labeling = self.inner.into_labeling();
        let n = labeling.len();
        labeling.iter().enumerate().fold(
            (Vec::<Vec<usize>>::new(), vec![None::<usize>; n]),
            |(mut groups, mut idx_of), (i, &r)| {
                let idx = *idx_of[r].get_or_insert_with(|| {
                    groups.push(vec![]);
                    groups.len() - 1
                });
                groups[idx].push(i);
                (groups, idx_of)
            },
        ).0
    }
}

pub struct KeyedUnionFind<X> where X: Eq + Hash {
    inner: UnionFind,
    keys: IndexSet<X>,
}

impl<X> KeyedUnionFind<X> where X: Eq + Hash {
    pub fn new() -> Self {
        Self { inner: UnionFind::new(0), keys: IndexSet::new() }
    }

    /// Insert a key, returning its index. Returns the existing index if `x` is already present.
    pub fn insert(&mut self, x: X) -> usize {
        let (idx, was_new) = self.keys.insert_full(x);
        if was_new {
            self.inner.extend(1);
        }
        idx
    }

    fn index_of(&self, x: &X) -> usize {
        self.keys.get_index_of(x).expect("key is not in the union-find")
    }

    fn element_at(&self, i: usize) -> &X {
        &self.keys[i]
    }

    pub fn size(&self) -> usize {
        self.inner.size()
    }

    pub fn contains(&self, x: &X) -> bool {
        self.keys.contains(x)
    }

    pub fn root(&self, x: &X) -> &X {
        let i = self.index_of(x);
        let j = self.inner.root(i);
        self.element_at(j)
    }

    pub fn is_same(&self, x: &X, y: &X) -> bool {
        self.root(x) == self.root(y)
    }

    pub fn union(&mut self, x: &X, y: &X) {
        let i = self.index_of(x);
        let j = self.index_of(y);
        self.inner.union(i, j);
    }

    pub fn into_disjoint(self) -> Vec<Vec<X>> {
        let Self { inner, keys } = self;
        let group = inner.into_disjoint();

        let mut map: HashMap<usize, X> = keys.into_iter().enumerate().collect();

        group.iter().map(|l|
            l.iter().map(|i| map.remove(i).unwrap()).collect()
        ).collect()
    }
}

impl<X> FromIterator<X> for KeyedUnionFind<X>
where X: Hash + Eq {
    fn from_iter<T: IntoIterator<Item = X>>(keys: T) -> Self {
        let keys: IndexSet<X> = keys.into_iter().collect();
        let inner = UnionFind::new(keys.len());
        Self { inner, keys }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use super::*;

    fn group_as_sets(g: Vec<Vec<usize>>) -> HashSet<Vec<usize>> {
        g.into_iter().map(|mut v| { v.sort(); v }).collect()
    }

    #[test]
    fn test() {
        let mut u = UnionFind::new(4);

        assert_eq!(u.size(), 4);
        assert!(!u.is_same(0, 1));
        assert!(!u.is_same(1, 2));
        assert!(!u.is_same(2, 3));
        assert_eq!(group_as_sets(u.clone().into_disjoint()), HashSet::from([vec![0], vec![1], vec![2], vec![3]]));

        u.union(0, 1);

        assert!( u.is_same(0, 1));
        assert!(!u.is_same(1, 2));
        assert!(!u.is_same(2, 3));
        assert_eq!(group_as_sets(u.clone().into_disjoint()), HashSet::from([vec![0, 1], vec![2], vec![3]]));

        u.union(2, 3);

        assert!( u.is_same(0, 1));
        assert!(!u.is_same(1, 2));
        assert!( u.is_same(2, 3));
        assert_eq!(group_as_sets(u.clone().into_disjoint()), HashSet::from([vec![0, 1], vec![2, 3]]));

        u.union(1, 3);

        assert!( u.is_same(0, 1));
        assert!( u.is_same(1, 2));
        assert!( u.is_same(2, 3));
        assert_eq!(group_as_sets(u.clone().into_disjoint()), HashSet::from([vec![0, 1, 2, 3]]));
    }

    fn keyed_with_unions(unions: &[(&'static str, &'static str)]) -> KeyedUnionFind<&'static str> {
        let mut u = KeyedUnionFind::from_iter(["a", "b", "c", "d"]);
        for (x, y) in unions { u.union(x, y); }
        u
    }

    fn sorted_disjoint(u: KeyedUnionFind<&'static str>) -> HashSet<Vec<&'static str>> {
        u.into_disjoint().into_iter().map(|mut v| { v.sort(); v }).collect()
    }

    #[test]
    fn test_hash_no_unions() {
        let u = keyed_with_unions(&[]);
        assert_eq!(u.size(), 4);
        assert!(!u.is_same(&"a", &"b"));
        assert!(!u.is_same(&"b", &"c"));
        assert!(!u.is_same(&"c", &"d"));
        assert_eq!(sorted_disjoint(u), HashSet::from([vec!["a"], vec!["b"], vec!["c"], vec!["d"]]));
    }

    #[test]
    fn test_hash_one_union() {
        let u = keyed_with_unions(&[("a", "b")]);
        assert!( u.is_same(&"a", &"b"));
        assert!(!u.is_same(&"b", &"c"));
        assert!(!u.is_same(&"c", &"d"));
        assert_eq!(sorted_disjoint(u), HashSet::from([vec!["a", "b"], vec!["c"], vec!["d"]]));
    }

    #[test]
    fn test_hash_two_unions() {
        let u = keyed_with_unions(&[("a", "b"), ("c", "d")]);
        assert!( u.is_same(&"a", &"b"));
        assert!(!u.is_same(&"b", &"c"));
        assert!( u.is_same(&"c", &"d"));
        assert_eq!(sorted_disjoint(u), HashSet::from([vec!["a", "b"], vec!["c", "d"]]));
    }

    #[test]
    fn test_hash_all_unioned() {
        let u = keyed_with_unions(&[("a", "b"), ("c", "d"), ("b", "d")]);
        assert!(u.is_same(&"a", &"b"));
        assert!(u.is_same(&"b", &"c"));
        assert!(u.is_same(&"c", &"d"));
        assert_eq!(sorted_disjoint(u), HashSet::from([vec!["a", "b", "c", "d"]]));
    }
}
