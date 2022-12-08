use std::ops::RangeInclusive;
use sprs::CsMat;

use crate::math::traits::{Ring, RingOps};
use crate::math::matrix::CsMatElem;
use super::complex::{ChainComplex, ChainComplexSparseD};

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

    pub fn d_matrix(&self, k: isize) -> CsMat<R> {
        let n = *self.range().end();
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
}

impl<R> ChainComplex for SimpleChainComplex<R> 
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> { 
    type R = R;
    type Generator = usize;

    fn range(&self) -> RangeInclusive<isize> {
        let n = self.generators.len() as isize;
        0..=(n - 1)
    }

    fn d_degree(&self) -> isize {
        self.d_degree
    }

    fn rank(&self, k: isize) -> usize {
        if self.range().contains(&k) { 
            let k = k as usize;
            self.generators[k].len()
        } else {
            0
        }
   }

    fn generators(&self, k: isize) -> Vec<&Self::Generator> {
        if self.range().contains(&k) { 
            let k = k as usize;
            self.generators[k].iter().collect()
        } else {
            vec![]
        }
    }

    fn differentiate(&self, k: isize, x: &Self::Generator) -> Vec<(Self::Generator, R)> {
        assert!(self.range().contains(&k));
        let v = self.vectorize_x(k, x);
        let d = self.d_matrix(k);
        let w = &d * &v;

        let gens = self.generators(k + self.d_degree());

        let res = w.iter().map(|(i, a)| { 
            (gens[i].clone(), a.clone())
        }).collect();

        res
    }
}

impl<R> ChainComplexSparseD for SimpleChainComplex<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> {}

#[cfg(test)]
pub mod tests { 
    use super::{ChainComplex, ChainComplexSparseD};
    use super::super::simple_complex::*;
    use crate::math::matrix::sparse::*;
    use sprs::CsMat;

    #[test]
    fn empty_complex() {
        let c = sample_empty();
        assert_eq!(c.range(), 0..=-1);
        assert_eq!(c.d_degree(), -1);
        assert_eq!(c.rank(0), 0);
    }

    #[test]
    fn complex_d3() {
        let c = sample_d3();

        assert_eq!(c.range(), 0..=3);
        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c.rank(0), 4);
        assert_eq!(c.rank(1), 6);
        assert_eq!(c.rank(2), 4);
        assert_eq!(c.rank(3), 1);

        c.check_d_all();
    }

    #[test]
    fn complex_s2() {
        let c = sample_s2();

        assert_eq!(c.range(), 0..=2);
        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c.rank(0), 4);
        assert_eq!(c.rank(1), 6);
        assert_eq!(c.rank(2), 4);

        c.check_d_all();
    }

    #[test]
    fn complex_t2() {
        let c = sample_t2();

        assert_eq!(c.range(), 0..=2);
        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c.rank(0), 9);
        assert_eq!(c.rank(1), 27);
        assert_eq!(c.rank(2), 18);

        c.check_d_all();
    }

    #[test]
    fn complex_rp2() {
        let c = sample_rp2();

        assert_eq!(c.range(), 0..=2);
        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c.rank(0), 6);
        assert_eq!(c.rank(1), 15);
        assert_eq!(c.rank(2), 10);

        c.check_d_all();
    }

    // below : test data // 

    pub fn sample_empty() -> SimpleChainComplex<i32> {
        SimpleChainComplex::new(vec![], -1)
    }

    pub fn sample_d3() -> SimpleChainComplex<i32> {
        SimpleChainComplex::new(
            vec![
                CsMat::csc_from_vec((4, 6), vec![-1, -1, 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 1]),
                CsMat::csc_from_vec((6, 4), vec![1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1] ),
                CsMat::csc_from_vec((4, 1), vec![-1, 1, -1, 1])
            ],
            -1
        )
    }

    pub fn sample_s2() -> SimpleChainComplex<i32> {
        SimpleChainComplex::new(
            vec![
                CsMat::csc_from_vec((4, 6), vec![-1, -1, 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 1]),
                CsMat::csc_from_vec((6, 4), vec![1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1] )
            ],
            -1
        )
    }

    pub fn sample_t2() -> SimpleChainComplex<i32> {
        SimpleChainComplex::new(
            vec![
                CsMat::csc_from_vec((9, 27), vec![-1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, -1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1]),
                CsMat::csc_from_vec((27, 18), vec![1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1])
            ],
            -1
        )
    }

    pub fn sample_rp2() -> SimpleChainComplex<i32> { 
        SimpleChainComplex::new(
            vec![
                CsMat::csc_from_vec((6, 15), vec![-1, -1, 0, 0, 0, 0, 0, -1, -1, 0, -1, 0, 0, 0, 0, 1, 0, -1, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 1, 1, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1] ),
                CsMat::csc_from_vec((15, 10), vec![1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1] ),
            ],
            -1
        )
    }
}