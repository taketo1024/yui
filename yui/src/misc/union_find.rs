use itertools::Itertools;

pub struct UnionFind { 
    roots: Vec<usize>
}

impl UnionFind { 
    pub fn new(n: usize) -> Self { 
        let roots = (0..n).collect();
        Self { roots }
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
        let n = self.roots.len();
        (0..n).into_group_map_by(|&i| self.roots[i]).into_values().collect()
    }
}