use std::collections::HashMap;
use std::hash::Hash;
use std::rc::Rc;

use itertools::Itertools;
use petgraph::unionfind::UnionFind as PgUnionFind;

/// A disjoint-set data structure over `0..n`. Backed by [`petgraph::unionfind::UnionFind`],
/// which uses path compression and union by rank for `O(α(n))` amortized operations.
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

    pub fn group(&self) -> Vec<Vec<usize>> {
        let n = self.size();
        (0..n).into_group_map_by(|&i| self.root(i))
            .into_iter()
            .sorted_by_key(|&(i, _)| i)
            .map(|(_, l)| l)
            .collect()
    }
}

pub struct KeyedUnionFind<X> where X: Eq + Hash {
    inner: UnionFind,
    keys: Vec<Rc<X>>,
    dict: HashMap<Rc<X>, usize>
}

impl<X> KeyedUnionFind<X> where X: Eq + Hash {
    pub fn new() -> Self {
        Self { inner: UnionFind::new(0), keys: vec![], dict: HashMap::new() }
    }

    pub fn insert(&mut self, x: X) -> usize {
        let i = self.size();
        let x = Rc::new(x);

        self.inner.extend(1);
        self.keys.push(Rc::clone(&x));
        self.dict.insert(x, i);

        i
    }

    fn index_of(&self, x: &X) -> usize {
        self.dict[x]
    }

    fn element_at(&self, i: usize) -> &X {
        &self.keys[i]
    }

    pub fn size(&self) -> usize {
        self.inner.size()
    }

    pub fn contains(&self, x: &X) -> bool {
        self.dict.contains_key(x)
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

    pub fn group(&self) -> Vec<Vec<&X>> {
        self.inner.group().iter().map(|l|
            l.iter().map(|&i|
                self.element_at(i)
            ).collect()
        ).collect()
    }

    pub fn into_group(mut self) -> Vec<Vec<X>> {
        let group = self.inner.group();
        let keys = std::mem::take(&mut self.keys);

        std::mem::drop(self);

        let mut map = keys.into_iter().enumerate().collect::<HashMap<_, _>>();

        group.iter().map(|l|
            l.iter().map(|i| {
                let x = map.remove(i).unwrap();

                Rc::into_inner(x).unwrap()
            }).collect()
        ).collect()
    }
}

impl<X> FromIterator<X> for KeyedUnionFind<X>
where X: Hash + Eq {
    fn from_iter<T: IntoIterator<Item = X>>(keys: T) -> Self {
        let keys = keys.into_iter().map(|e| Rc::new(e)).collect_vec();
        let dict = keys.iter().enumerate().map(|(i, e)| (Rc::clone(e), i)).collect();
        let n = keys.len();
        let inner = UnionFind::new(n);
        Self { inner, keys, dict }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use super::*;

    fn group_as_sets(g: Vec<Vec<usize>>) -> HashSet<Vec<usize>> {
        g.into_iter().map(|mut v| { v.sort(); v }).collect()
    }

    fn group_as_str_sets<'a>(g: Vec<Vec<&'a &'a str>>) -> HashSet<Vec<&'a str>> {
        g.into_iter().map(|v| {
            let mut v: Vec<&str> = v.into_iter().copied().collect();
            v.sort();
            v
        }).collect()
    }

    #[test]
    fn test() {
        let mut u = UnionFind::new(4);

        assert_eq!(u.size(), 4);
        assert!(!u.is_same(0, 1));
        assert!(!u.is_same(1, 2));
        assert!(!u.is_same(2, 3));
        assert_eq!(group_as_sets(u.group()), HashSet::from([vec![0], vec![1], vec![2], vec![3]]));

        u.union(0, 1);

        assert!( u.is_same(0, 1));
        assert!(!u.is_same(1, 2));
        assert!(!u.is_same(2, 3));
        assert_eq!(group_as_sets(u.group()), HashSet::from([vec![0, 1], vec![2], vec![3]]));

        u.union(2, 3);

        assert!( u.is_same(0, 1));
        assert!(!u.is_same(1, 2));
        assert!( u.is_same(2, 3));
        assert_eq!(group_as_sets(u.group()), HashSet::from([vec![0, 1], vec![2, 3]]));

        u.union(1, 3);

        assert!( u.is_same(0, 1));
        assert!( u.is_same(1, 2));
        assert!( u.is_same(2, 3));
        assert_eq!(group_as_sets(u.group()), HashSet::from([vec![0, 1, 2, 3]]));
    }

    #[test]
    fn test_hash() {
        let mut u = KeyedUnionFind::from_iter(["a", "b", "c", "d"]);

        assert_eq!(u.size(), 4);
        assert!(!u.is_same(&"a", &"b"));
        assert!(!u.is_same(&"b", &"c"));
        assert!(!u.is_same(&"c", &"d"));
        assert_eq!(group_as_str_sets(u.group()), HashSet::from([vec!["a"], vec!["b"], vec!["c"], vec!["d"]]));

        u.union(&"a", &"b");

        assert!( u.is_same(&"a", &"b"));
        assert!(!u.is_same(&"b", &"c"));
        assert!(!u.is_same(&"c", &"d"));
        assert_eq!(group_as_str_sets(u.group()), HashSet::from([vec!["a", "b"], vec!["c"], vec!["d"]]));

        u.union(&"c", &"d");

        assert!( u.is_same(&"a", &"b"));
        assert!(!u.is_same(&"b", &"c"));
        assert!( u.is_same(&"c", &"d"));
        assert_eq!(group_as_str_sets(u.group()), HashSet::from([vec!["a", "b"], vec!["c", "d"]]));

        u.union(&"b", &"d");

        assert!(u.is_same(&"a", &"b"));
        assert!(u.is_same(&"b", &"c"));
        assert!(u.is_same(&"c", &"d"));
        assert_eq!(group_as_str_sets(u.group()), HashSet::from([vec!["a", "b", "c", "d"]]));

        let into_groups: Vec<Vec<&str>> = u.into_group()
            .into_iter()
            .map(|mut v| { v.sort(); v })
            .collect();
        assert_eq!(into_groups, vec![vec!["a", "b", "c", "d"]]);
    }
}
