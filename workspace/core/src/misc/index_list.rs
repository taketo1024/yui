use std::fmt::Debug;
use std::hash::Hash;
use std::ops::Index;

use bimap::BiHashMap;
use ahash::RandomState as ARandomState;
use itertools::Itertools;

#[derive(Default, Clone)]
pub struct IndexList<E>
where E: Eq + Hash {
    data: BiHashMap<usize, E, ARandomState, ARandomState>
}

impl<E> IndexList<E>
where E: Eq + Hash {
    pub fn new() -> Self {
        Self { data: BiHashMap::with_hashers(ARandomState::new(), ARandomState::new()) }
    }

    pub fn len(&self) -> usize { 
        self.data.len()
    }

    pub fn contains(&self, x: &E) -> bool { 
        self.data.contains_right(x)
    }

    pub fn index_of(&self, x: &E) -> Option<usize> {
        self.data.get_by_right(x).cloned()
    }

    pub fn iter(&self) -> impl Iterator<Item = &E> { 
        let n = self.len();
        (0..n).map(|i| 
            self.data.get_by_left(&i).unwrap()
        )
    }
}

impl<E> Debug for IndexList<E>
where E: Eq + Hash + Debug {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.iter().collect_vec())
    }
}

impl<E> FromIterator<E> for IndexList<E>
where E: Eq + Hash {
    fn from_iter<T: IntoIterator<Item = E>>(iter: T) -> Self {
        let data = iter.into_iter().enumerate().collect();
        Self { data }
    }
}

impl<E> IntoIterator for IndexList<E>
where E: Eq + Hash {
    type Item = E;
    type IntoIter = std::vec::IntoIter<E>;

    fn into_iter(self) -> Self::IntoIter {
        let mut data = self.data;
        let n = data.len();
        (0..n).map(move |i| 
            data.remove_by_left(&i).unwrap().1
        ).collect_vec().into_iter()
    }
}

impl<E> Index<usize> for IndexList<E>
where E: Eq + Hash {
    type Output = E;

    fn index(&self, index: usize) -> &Self::Output {
        self.data.get_by_left(&index).unwrap()
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn display() {
        let list = IndexList::from_iter([3,7,4,5,1]);
        assert_eq!(format!("{:?}", list), "[3, 7, 4, 5, 1]");
    }
    
    #[test]
    fn index_of() {
        let list = IndexList::from_iter([3,7,4,5,1]);
        assert_eq!(list[0], 3);
        assert_eq!(list[4], 1);
        assert_eq!(list.index_of(&3), Some(0));
        assert_eq!(list.index_of(&4), Some(2));
        assert_eq!(list.index_of(&2), None);
    }

    #[test]
    fn into_iter() {
        let list = IndexList::from_iter([3,7,4,5,1]);
        let into = list.into_iter().collect_vec();
        assert_eq!(into, vec![3,7,4,5,1]);
    }
}