use itertools::Itertools;

pub struct UnionFind { 
    p: Vec<usize>
}

impl UnionFind { 
    pub fn new(n: usize) -> Self { 
        Self { p: (0..n).collect() }
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

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::UnionFind;
 
    #[test]
    fn test() { 
        let mut u = UnionFind::new(4);
        
        assert_eq!(u.size(), 4);
        assert_eq!(&u.p, &vec![0,1,2,3]);
        assert_eq!((0..4).map(|i| u.root(i)).collect_vec(), vec![0,1,2,3]);

        assert!(!u.is_same(0, 1));
        assert!(!u.is_same(1, 2));
        assert!(!u.is_same(2, 3));

        u.union(0, 1);

        assert_eq!(&u.p, &vec![0,0,2,3]);
        assert_eq!((0..4).map(|i| u.root(i)).collect_vec(), vec![0,0,2,3]);

        assert!( u.is_same(0, 1));
        assert!(!u.is_same(1, 2));
        assert!(!u.is_same(2, 3));

        u.union(2, 3);

        assert_eq!(&u.p, &vec![0,0,2,2]);
        assert_eq!((0..4).map(|i| u.root(i)).collect_vec(), vec![0,0,2,2]);

        assert!( u.is_same(0, 1));
        assert!(!u.is_same(1, 2));
        assert!( u.is_same(2, 3));

        u.union(1, 3);

        assert_eq!(&u.p, &vec![0,0,0,2]);
        assert_eq!((0..4).map(|i| u.root(i)).collect_vec(), vec![0,0,0,0]);

        assert!( u.is_same(0, 1));
        assert!( u.is_same(1, 2));
        assert!( u.is_same(2, 3));

    }
}