use std::hash::Hash;
use std::ops::Index;

use bimap::BiHashMap;

#[derive(Default, Clone)]
pub struct IndexList<E>
where E: Eq + Hash {
    data: BiHashMap<usize, E>
}

impl<E> IndexList<E>
where E: Eq + Hash {
    pub fn new() -> Self {
        Self { data: BiHashMap::new() }
    }

    pub fn len(&self) -> usize { 
        self.data.len()
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

    pub fn into_iter(self) -> impl Iterator<Item = E> { 
        let mut data = self.data;
        let n = data.len();
        (0..n).map(move |i| 
            data.remove_by_left(&i).unwrap().1
        )
    }
}

impl<E> FromIterator<E> for IndexList<E>
where E: Eq + Hash {
    fn from_iter<T: IntoIterator<Item = E>>(iter: T) -> Self {
        let data = iter.into_iter().enumerate().collect();
        Self { data }
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
    fn test() {
        let list = IndexList::from_iter([3,7,4,5,1]);
        assert_eq!(list[0], 3);
        assert_eq!(list[4], 1);
        assert_eq!(list.index_of(&3), Some(0));
        assert_eq!(list.index_of(&4), Some(2));
        assert_eq!(list.index_of(&2), None);
    }
}