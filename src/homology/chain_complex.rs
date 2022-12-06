use std::ops::RangeInclusive;
use std::collections::HashMap;
use std::hash::Hash;
use sprs::{CsMat, CsVec, TriMat};
use num_traits::One;
use crate::math::traits::{Ring, RingOps};
use crate::matrix::sparse::*;

pub trait ChainGenerator: PartialEq + Eq + Hash {}

impl<T> ChainGenerator for T 
where T: PartialEq + Eq + Hash {}

pub trait ChainComplex 
where 
    Self::R: Ring + CsMatElem, 
    for<'x> &'x Self::R: RingOps<Self::R> 
{
    type R;
    type Generator: ChainGenerator;

    fn hdeg_range(&self) -> RangeInclusive<isize>;
    fn generators(&self, k: isize) -> Vec<&Self::Generator>;

    fn rank(&self, k: isize) -> usize { 
        if self.hdeg_range().contains(&k) { 
            self.generators(k).len()
        } else {
            0
        }
    }

    fn d_degree(&self) -> isize { 1 } 
    fn d_matrix(&self, k: isize) -> CsMat<Self::R>;
    fn differentiate(&self, k: isize, x:&Self::Generator) -> Vec<(&Self::Generator, Self::R)>;

    fn make_d_matrix(&self, k: isize) -> CsMat<Self::R> {
        let source = self.generators(k);
        let target = self.generators(k + self.d_degree());
        let (m, n) = (target.len(), source.len());

        let t_ind = target.into_iter()
                .enumerate()
                .map(|(i, y)| (y, i))
                .collect::<HashMap<_, _>>();

        let mut trip = TriMat::new((m, n));

        for (j, x) in source.iter().enumerate() {
            let ys = self.differentiate(k, x);
            for (y, a) in ys {
                if !t_ind.contains_key(y) { continue }
                let i = t_ind[y];
                trip.add_triplet(i, j, a);
            }
        }

        trip.to_csc()
    }

    fn vectorize_x(&self, k: isize, x:&Self::Generator) -> CsVec<Self::R> {
        self.vectorize(k, hashmap![x => Self::R::one()])
    }

    fn vectorize(&self, k: isize, z:HashMap<&Self::Generator, Self::R>) -> CsVec<Self::R> {
        let gens = self.generators(k);
        let n = gens.len();

        let mut v_ind: Vec<usize> = vec![];
        let mut v_val: Vec<Self::R> = vec![];

        for (x, a) in z { 
            let Some(i) = gens.iter().position(|&z| x == z) else { 
                continue 
            };
            v_ind.push(i);
            v_val.push(a);
        }

        CsVec::new(n, v_ind, v_val)
    }

    fn check_d_at(&self, k: isize) { 
        let d1 = self.d_matrix(k);
        let d2 = self.d_matrix(k + self.d_degree());
        let res = &d2 * &d1;
        assert!( res.is_zero() );
    }

    fn check_d_all(&self) {
        for k in self.hdeg_range() { 
            let k2 = k + self.d_degree();
            if !self.hdeg_range().contains(&k2) { continue }

            self.check_d_at(k);
        }
    }
}

pub struct SimpleChainComplex<R> 
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> { 
    generators: Vec<Vec<usize>>,
    d_matrices: Vec<CsMat<R>>,
    d_degree: isize
}

impl<R> SimpleChainComplex<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> { 
    pub fn new(d_matrices: Vec<CsMat<R>>, d_degree: isize) -> Self {
        assert!(d_degree == 1 || d_degree == -1);

        let mut count = 0;
        let mut gens = |n: usize| { 
            let g = (count .. count + n).into_iter().collect();
            count += n;
            g
        };

        let mut generators: Vec<_> = d_matrices.iter().map(|d| {
            let n = if d_degree > 0 { d.cols() } else { d.rows() };
            gens(n)
        }).collect();

        if let Some(d) = d_matrices.last() { 
            let n = if d_degree > 0 { d.rows() } else { d.cols() };
            generators.push(gens(n));
        }

        SimpleChainComplex { generators, d_matrices, d_degree }
    }
}

impl<R> ChainComplex for SimpleChainComplex<R> 
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> { 
    type R = R;
    type Generator = usize;

    fn hdeg_range(&self) -> RangeInclusive<isize> {
        let n = self.generators.len() as isize;
        0..=(n - 1)
    }

    fn d_degree(&self) -> isize {
        self.d_degree
    }

    fn rank(&self, k: isize) -> usize {
        if self.hdeg_range().contains(&k) { 
            let k = k as usize;
            self.generators[k].len()
        } else {
            0
        }
   }

    fn generators(&self, k: isize) -> Vec<&Self::Generator> {
        if self.hdeg_range().contains(&k) { 
            let k = k as usize;
            self.generators[k].iter().collect()
        } else {
            vec![]
        }
    }

    fn d_matrix(&self, k: isize) -> CsMat<R> {
        let n = *self.hdeg_range().end();
        let (range, idx) = if self.d_degree > 0 { 
            (0 ..= n - 1, k as usize)
        } else {
            (1 ..= n, (k - 1) as usize)
        };
        if range.contains(&k) { 
            self.d_matrices[idx].clone()
        } else {
            let m = self.rank(k + self.d_degree);
            let n = self.rank(k);
            CsMat::zero((m, n))
        }
    }

    fn differentiate(&self, k: isize, x: &Self::Generator) -> Vec<(&Self::Generator, R)> {
        assert!(self.hdeg_range().contains(&k));
        let v = self.vectorize_x(k, x);
        let d = self.d_matrix(k);
        let w = &d * &v;

        let gens = self.generators(k + self.d_degree());

        let res = w.iter().map(|(i, a)| { 
            (gens[i], a.clone())
        }).collect();

        res
    }
}

#[cfg(test)]
pub mod tests { 
    use super::{SimpleChainComplex, ChainComplex};
    use num_traits::Zero;
    use sprs::{CsMat, TriMat};

    #[test]
    fn empty_complex() {
        let c = sample_empty();
        assert_eq!(c.hdeg_range(), 0..=-1);
        assert_eq!(c.d_degree(), -1);
        assert_eq!(c.rank(0), 0);
    }

    #[test]
    fn complex_d3() {
        let c = sample_d3();

        assert_eq!(c.hdeg_range(), 0..=3);
        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c.rank(0), 4);
        assert_eq!(c.rank(1), 6);
        assert_eq!(c.rank(2), 4);
        assert_eq!(c.rank(3), 1);

        c.check_d_all();
    }

    #[test]
    fn complex_t2() {
        let c = sample_t2();

        assert_eq!(c.hdeg_range(), 0..=2);
        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c.rank(0), 9);
        assert_eq!(c.rank(1), 27);
        assert_eq!(c.rank(2), 18);

        c.check_d_all();
    }


    #[test]
    fn complex_rp2() {
        let c = sample_rp2();

        assert_eq!(c.hdeg_range(), 0..=2);
        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c.rank(0), 6);
        assert_eq!(c.rank(1), 15);
        assert_eq!(c.rank(2), 10);

        c.check_d_all();
    }

    // below : test data // 

    pub fn mat(size:(usize, usize), data: Vec<i32>) -> CsMat<i32> {
        let n = size.1;
        let mut trip = TriMat::new(size);
        for (k, a) in data.into_iter().enumerate() {
            if a.is_zero() { continue }
            let (i, j) = (k / n, k % n);
            trip.add_triplet(i, j, a);
        }
        trip.to_csc()
    }

    pub fn sample_empty() -> SimpleChainComplex<i32> {
        SimpleChainComplex::new(vec![], -1)
    }

    pub fn sample_d3() -> SimpleChainComplex<i32> {
        SimpleChainComplex::new(
            vec![
                mat((4, 6), vec![-1, -1, 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 1]),
                mat((6, 4), vec![1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1] ),
                mat((4, 1), vec![-1, 1, -1, 1])
            ],
            -1
        )
    }

    pub fn sample_t2() -> SimpleChainComplex<i32> {
        SimpleChainComplex::new(
            vec![
                mat((9, 27), vec![-1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, -1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1]),
                mat((27, 18), vec![1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1])
            ],
            -1
        )
    }

    pub fn sample_rp2() -> SimpleChainComplex<i32> { 
        SimpleChainComplex::new(
            vec![
                mat((6, 15), vec![-1, -1, 0, 0, 0, 0, 0, -1, -1, 0, -1, 0, 0, 0, 0, 1, 0, -1, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 1, 1, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1] ),
                mat((15, 10), vec![1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1] ),
            ],
            -1
        )
    }
}