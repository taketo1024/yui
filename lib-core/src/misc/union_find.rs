use std::collections::HashMap;
use std::hash::Hash;
use std::rc::Rc;

use itertools::Itertools;

pub struct UnionFind { 
    p: Vec<usize>
}

impl UnionFind { 
    pub fn new(n: usize) -> Self { 
        Self { p: (0..n).collect() }
    }

    pub fn extend(&mut self, l: usize) { 
        let n = self.p.len();
        self.p.extend(n .. n + l);
    }

    pub fn size(&self) -> usize { 
        self.p.len()
    }

    pub fn root(&self, i: usize) -> usize { 
        let p = self.p[i];
        if p == i { 
            i
        } else { 
            self.root(p)
        }
    }
    
    pub fn is_same(&self, i: usize, j: usize) -> bool { 
        self.root(i) == self.root(j)
    }

    pub fn union(&mut self, i: usize, j: usize) {
        use std::cmp::Ordering::*;
        let ri = self.root(i);
        let rj = self.root(j);

        match usize::cmp(&ri, &rj) {
            Less    => self.p[rj] = ri,
            Equal   => (),
            Greater => self.p[ri] = rj,
        }
    }

    pub fn group(&self) -> Vec<Vec<usize>> { 
        let n = self.size();
        (0..n).into_group_map_by(|&i| self.root(i)).into_iter().sorted_by_key(|&(i, _)| i).map(|(_, l)| l).collect()
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
                let x = Rc::into_inner(x).unwrap();
                x
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
    use itertools::Itertools;

    use super::*;
 
    #[test]
    fn test() { 
        let mut u = UnionFind::new(4);
        
        assert_eq!(u.size(), 4);
        assert_eq!(&u.p, &vec![0,1,2,3]);
        assert_eq!((0..4).map(|i| u.root(i)).collect_vec(), vec![0,1,2,3]);

        assert!(!u.is_same(0, 1));
        assert!(!u.is_same(1, 2));
        assert!(!u.is_same(2, 3));
        assert_eq!(u.group(), vec![vec![0], vec![1], vec![2], vec![3]]);

        u.union(0, 1);

        assert_eq!(&u.p, &vec![0,0,2,3]);
        assert_eq!((0..4).map(|i| u.root(i)).collect_vec(), vec![0,0,2,3]);

        assert!( u.is_same(0, 1));
        assert!(!u.is_same(1, 2));
        assert!(!u.is_same(2, 3));
        assert_eq!(u.group(), vec![vec![0, 1], vec![2], vec![3]]);

        u.union(2, 3);

        assert_eq!(&u.p, &vec![0,0,2,2]);
        assert_eq!((0..4).map(|i| u.root(i)).collect_vec(), vec![0,0,2,2]);

        assert!( u.is_same(0, 1));
        assert!(!u.is_same(1, 2));
        assert!( u.is_same(2, 3));
        assert_eq!(u.group(), vec![vec![0, 1], vec![2, 3]]);

        u.union(1, 3);

        assert_eq!(&u.p, &vec![0,0,0,2]);
        assert_eq!((0..4).map(|i| u.root(i)).collect_vec(), vec![0,0,0,0]);

        assert!( u.is_same(0, 1));
        assert!( u.is_same(1, 2));
        assert!( u.is_same(2, 3));
        assert_eq!(u.group(), vec![vec![0, 1, 2, 3]]);
    }

    #[test]
    fn test_hash() { 
        let mut u = KeyedUnionFind::from_iter(["a", "b", "c", "d"]);
        
        assert_eq!(u.size(), 4);
        assert_eq!(u.root(&"a"), &"a");
        assert_eq!(u.root(&"b"), &"b");
        assert_eq!(u.root(&"c"), &"c");
        assert_eq!(u.root(&"d"), &"d");

        assert!(!u.is_same(&"a", &"b"));
        assert!(!u.is_same(&"b", &"c"));
        assert!(!u.is_same(&"c", &"d"));
        assert_eq!(u.group(), vec![vec![&"a"], vec![&"b"], vec![&"c"], vec![&"d"]]);

        u.union(&"a", &"b");

        assert_eq!(u.root(&"a"), &"a");
        assert_eq!(u.root(&"b"), &"a");
        assert_eq!(u.root(&"c"), &"c");
        assert_eq!(u.root(&"d"), &"d");

        assert!( u.is_same(&"a", &"b"));
        assert!(!u.is_same(&"b", &"c"));
        assert!(!u.is_same(&"c", &"d"));
        assert_eq!(u.group(), vec![vec![&"a", &"b"], vec![&"c"], vec![&"d"]]);

        u.union(&"c", &"d");

        assert_eq!(u.root(&"a"), &"a");
        assert_eq!(u.root(&"b"), &"a");
        assert_eq!(u.root(&"c"), &"c");
        assert_eq!(u.root(&"d"), &"c");

        assert!( u.is_same(&"a", &"b"));
        assert!(!u.is_same(&"b", &"c"));
        assert!( u.is_same(&"c", &"d"));
        assert_eq!(u.group(), vec![vec![&"a", &"b"], vec![&"c", &"d"]]);

        u.union(&"b", &"d");

        assert_eq!(u.root(&"a"), &"a");
        assert_eq!(u.root(&"b"), &"a");
        assert_eq!(u.root(&"c"), &"a");
        assert_eq!(u.root(&"d"), &"a");

        assert!(u.is_same(&"a", &"b"));
        assert!(u.is_same(&"b", &"c"));
        assert!(u.is_same(&"c", &"d"));

        assert_eq!(u.group(), vec![vec![&"a", &"b", &"c", &"d"]]);
        assert_eq!(u.into_group(), vec![vec!["a", "b", "c", "d"]]);
    }
}