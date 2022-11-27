use std::ops::RangeInclusive;
use std::collections::HashMap;
use std::hash::Hash;
use sprs::{CsMat, CsVec, TriMat};
use num_traits::{One, Zero};
use crate::matrix::CsMatElem;

pub trait ChainGenerator: PartialEq + Eq + Hash {}

pub trait ChainComplex {
    type R: CsMatElem;
    type Generator: ChainGenerator;

    fn hdeg_range(&self) -> RangeInclusive<isize>;
    fn in_hdeg_range(&self, k: isize) -> bool { self.hdeg_range().contains(&k) }
    fn generators(&self, k: isize) -> &Vec<Self::Generator>;
    fn rank(&self, k: isize) -> u32 { if self.in_hdeg_range(k) { self.generators(k).len() as u32 } else { 0 } } 

    fn d_degree(&self) -> isize { 1 } 
    fn d_matrix(&self, k: isize) -> &CsMat<Self::R>;
    fn differentiate(&self, k: isize, x:&Self::Generator) -> Vec<(&Self::Generator, Self::R)>;

    fn make_d_matrix(&self, k: isize) -> CsMat<Self::R> {
        let source = self.generators(k);
        let target = self.generators(k + self.d_degree());
        let (m, n) = (target.len(), source.len());

        let t_ind = target.iter()
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
            let Some(i) = gens.iter().position(|z| x == z) else { continue };
            v_ind.push(i);
            v_val.push(a);
        }

        CsVec::new(n, v_ind, v_val)
    }

    fn check_d_at(&self, k: isize) { 
        let d1 = self.d_matrix(k);
        let d2 = self.d_matrix(k + self.d_degree());
        let res = d2 * d1;
        assert!( res.data().iter().all(|a| a.is_zero()) );
    }

    fn check_d_all(&self) {
        for k in self.hdeg_range() { 
            let k2 = k + self.d_degree();
            if !self.hdeg_range().contains(&k2) { continue }

            self.check_d_at(k);
        }
    }
}

pub struct SimpleChainComplex<X, R> 
where 
    X: ChainGenerator, 
    R: CsMatElem
{ 
    generators: Vec<Vec<X>>,
    d_matrices: Vec<CsMat<R>>,
    d_degree: isize
}

impl<R, X> ChainComplex for SimpleChainComplex<X, R> 
where 
    X: ChainGenerator, 
    R: CsMatElem
{ 
    type R = R;
    type Generator = X;

    fn hdeg_range(&self) -> RangeInclusive<isize> {
        let n = self.generators.len() as isize;
        0..=(n - 1)
    }

    fn d_degree(&self) -> isize {
        self.d_degree
    }

    fn generators(&self, k: isize) -> &Vec<X> {
        assert!(self.hdeg_range().contains(&k));
        &self.generators[k as usize]
    }

    fn d_matrix(&self, k: isize) -> &CsMat<R> {
        assert!(self.hdeg_range().contains(&k));
        &self.d_matrices[k as usize]
    }

    fn differentiate(&self, k: isize, x: &X) -> Vec<(&X, R)> {
        assert!(self.hdeg_range().contains(&k));
        let v = self.vectorize_x(k, x);
        let d = self.d_matrix(k);
        let w = d * &v;

        let gens = self.generators(k + self.d_degree());

        let res = w.iter().map(|(i, a)| { 
            (&gens[i], a.clone())
        }).collect();

        res
    }
}

#[cfg(test)]
mod tests { 
    use super::{SimpleChainComplex, ChainGenerator, ChainComplex};
    use std::ops::Range;
    use num_traits::Zero;
    use sprs::{CsMat, TriMat};

    #[derive(PartialEq, Eq, Hash)]
    struct X { val: i32 }
    impl ChainGenerator for X {}

    fn gens(range: Range<i32>) -> Vec<X> { 
        range.map(|i| X{ val: i }).collect() 
    }

    fn mat(size:(usize, usize), data: Vec<i32>) -> CsMat<i32> {
        let n = size.1;
        let mut trip = TriMat::new(size);
        for (k, a) in data.into_iter().enumerate() {
            if a.is_zero() { continue }
            let (i, j) = (k / n, k % n);
            trip.add_triplet(i, j, a);
        }
        trip.to_csc()
    }

    #[test]
    fn empty_complex() {
        let c: SimpleChainComplex<X, i32> = SimpleChainComplex { 
            generators: vec![
                gens(0..0),
                gens(0..0),
                gens(0..0),
            ],
            d_matrices: vec![
                mat((0, 0), vec![]),
                mat((0, 0), vec![]),
                mat((0, 0), vec![]),
            ],
            d_degree: -1
        };

        assert_eq!(c.hdeg_range(), 0..=2);
        assert_eq!(c.d_degree(), -1);
    }

    #[test]
    fn complex_d3() {
        let c: SimpleChainComplex<X, i32> = SimpleChainComplex { 
            generators: vec![
                gens(0..4),
                gens(4..10),
                gens(10..14),
                gens(14..15)
            ],
            d_matrices: vec![
                mat((0, 4), vec![]),
                mat((4, 6), vec![-1, -1, 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 1]),
                mat((6, 4), vec![1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1] ),
                mat((4, 1), vec![-1, 1, -1, 1])
            ],
            d_degree: -1
        };

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
        let c: SimpleChainComplex<X, i32> = SimpleChainComplex { 
            generators: vec![
                gens(0..9),
                gens(9..36),
                gens(36..54)
            ],
            d_matrices: vec![
                mat((0, 9), vec![]),
                mat((9, 27), vec![-1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, -1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1]),
                mat((27, 18), vec![1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1])
            ],
            d_degree: -1
        };

        assert_eq!(c.hdeg_range(), 0..=2);
        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c.rank(0), 9);
        assert_eq!(c.rank(1), 27);
        assert_eq!(c.rank(2), 18);

        c.check_d_all();
    }


    #[test]
    fn complex_rp2() {
        let c: SimpleChainComplex<X, i32> = SimpleChainComplex { 
            generators: vec![
                gens(0..6),
                gens(6..21),
                gens(21..31)
            ],
            d_matrices: vec![
                mat((0, 6), vec![]),
                mat((6, 15), vec![-1, -1, 0, 0, 0, 0, 0, -1, -1, 0, -1, 0, 0, 0, 0, 1, 0, -1, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 1, 1, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1] ),
                mat((15, 10), vec![1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1] ),
            ],
            d_degree: -1
        };

        assert_eq!(c.hdeg_range(), 0..=2);
        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c.rank(0), 6);
        assert_eq!(c.rank(1), 15);
        assert_eq!(c.rank(2), 10);

        c.check_d_all();
    }
}