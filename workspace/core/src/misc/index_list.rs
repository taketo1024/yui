use std::hash::Hash;
use std::ops::Index;

use bimap::BiHashMap;

#[derive(Default)]
pub struct IndexList<E>
where E: Eq + Hash {
    data: BiHashMap<usize, E>
}

impl<E> IndexList<E>
where E: Eq + Hash {
    pub fn new<Itr>(gens: Itr) -> Self
    where Itr: Iterator<Item = E> { 
        let data = gens.enumerate().collect();
        Self { data }
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

impl<E> Index<usize> for IndexList<E>
where E: Eq + Hash {
    type Output = E;

    fn index(&self, index: usize) -> &Self::Output {
        self.data.get_by_left(&index).unwrap()
    }
}