use itertools::Itertools;

pub struct UnionFind { 
    roots: Vec<usize>
}

impl UnionFind { 
    pub fn new(n: usize) -> Self { 
        let roots = (0..n).collect();
        Self { roots }
    }

    pub fn size(&self) -> usize { 
        self.roots.len()
    }

    pub fn is_same(&self, i: usize, j: usize) -> bool { 
        self.roots[i] == self.roots[j]
    }

    pub fn union(&mut self, i: usize, j: usize) {
        let ri = self.roots[i];
        let rj = self.roots[j];

        match usize::cmp(&ri, &rj) {
            std::cmp::Ordering::Less    => self.roots[j] = ri,
            std::cmp::Ordering::Equal   => (),
            std::cmp::Ordering::Greater => self.roots[i] = rj,
        }
    }

    pub fn group(&self) -> Vec<Vec<usize>> { 
        let n = self.size();
        (0..n).into_group_map_by(|&i| self.roots[i]).into_values().collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::UnionFind;
 
    #[test]
    fn test() { 
        let mut u = UnionFind::new(4);
        
        assert_eq!(u.size(), 4);
        assert!(!u.is_same(0, 1));
        assert!(!u.is_same(0, 2));
        assert!(!u.is_same(1, 2));
        assert!(!u.is_same(0, 3));

        u.union(0, 1);

        assert!( u.is_same(0, 1));
        assert!(!u.is_same(0, 2));
        assert!(!u.is_same(1, 2));
        assert!(!u.is_same(0, 3));

        u.union(1, 2);

        assert!( u.is_same(0, 1));
        assert!( u.is_same(0, 2));
        assert!( u.is_same(1, 2));
        assert!(!u.is_same(0, 3));
    }
}