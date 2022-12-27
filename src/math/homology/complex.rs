use std::collections::HashMap;
use sprs::{CsMat, CsVec, TriMat};
use num_traits::One;
use crate::math::traits::{Ring, RingOps};
use crate::math::matrix::sparse::*;
use crate::utils::collections::hashmap;
use super::base::{Graded, ChainGenerator};

pub trait ChainComplex: Graded
where 
    Self::R: Ring + CsMatElem, for<'x> &'x Self::R: RingOps<Self::R> 
{
    type R;
    type Generator: ChainGenerator;

    // -- must override -- //
    fn generators(&self, k: Self::Index) -> Vec<&Self::Generator>;
    fn d_degree(&self) -> Self::Index;
    fn differentiate(&self, k: Self::Index, x:&Self::Generator) -> Vec<(Self::Generator, Self::R)>;
    fn d_matrix(&self, k: Self::Index) -> CsMat<Self::R>;

    // -- convenient methods -- //
    fn rank(&self, k: Self::Index) -> usize { 
        if self.in_range(k) { 
            self.generators(k).len()
        } else {
            0
        }
    }

    fn vectorize_x(&self, k: Self::Index, x:&Self::Generator) -> CsVec<Self::R> {
        self.vectorize(k, hashmap![x => Self::R::one()])
    }

    fn vectorize(&self, k: Self::Index, z:HashMap<&Self::Generator, Self::R>) -> CsVec<Self::R> {
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

    fn impl_d_matrix_from_differentiate(&self, k: Self::Index) -> CsMat<Self::R> {
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
                if !t_ind.contains_key(&y) { continue }
                let i = t_ind[&y];
                trip.add_triplet(i, j, a);
            }
        }

        trip.to_csc()
    }

    fn impl_differentiate_from_d_matrix(&self, k: Self::Index, x: &Self::Generator) -> Vec<(Self::Generator, Self::R)> {
        assert!(self.in_range(k));
        let v = self.vectorize_x(k, x);
        let d = self.d_matrix(k);
        let w = &d * &v;

        let gens = self.generators(k + self.d_degree());

        let res = w.iter().map(|(i, a)| { 
            (gens[i].clone(), a.clone())
        }).collect();

        res
    }

    fn check_d_at(&self, k: Self::Index) { 
        let d1 = self.d_matrix(k);
        let d2 = self.d_matrix(k + self.d_degree());
        let res = &d2 * &d1;
        assert!( res.is_zero() );
    }

    fn check_d_all(&self) {
        for k in self.range() { 
            let k2 = k + self.d_degree();
            if !self.in_range(k2) { continue }
            self.check_d_at(k);
        }
    }
}

#[cfg(test)]
pub mod tests { 
    use std::iter::Rev;
    use std::ops::RangeInclusive;
    use either::Either;
    use sprs::CsMat;
    use super::{ChainComplex, Graded};
    use crate::math::traits::{Ring, RingOps};
    use crate::math::matrix::sparse::*;
    use crate::math::matrix::CsMatElem;

    #[test]
    fn zero_complex() {
        let c = TestChainComplex::<i32>::zero();
        assert_eq!(c.d_degree(), -1);
        assert_eq!(c.rank(0), 0);
    }

    #[test]
    fn d3() {
        let c = TestChainComplex::<i32>::d3();

        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c.rank(0), 4);
        assert_eq!(c.rank(1), 6);
        assert_eq!(c.rank(2), 4);
        assert_eq!(c.rank(3), 1);

        c.check_d_all();
    }

    #[test]
    fn s2() {
        let c = TestChainComplex::<i32>::s2();

        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c.rank(0), 4);
        assert_eq!(c.rank(1), 6);
        assert_eq!(c.rank(2), 4);

        c.check_d_all();
    }

    #[test]
    fn t2() {
        let c = TestChainComplex::<i32>::t2();

        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c.rank(0), 9);
        assert_eq!(c.rank(1), 27);
        assert_eq!(c.rank(2), 18);

        c.check_d_all();
    }

    #[test]
    fn rp2() {
        let c = TestChainComplex::<i32>::rp2();

        assert_eq!(c.d_degree(), -1);
        
        assert_eq!(c.rank(0), 6);
        assert_eq!(c.rank(1), 15);
        assert_eq!(c.rank(2), 10);

        c.check_d_all();
    }

    // below : test data // 

    pub type X = usize;
    pub struct TestChainComplex<R> 
    where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> {
        range: RangeInclusive<isize>,
        generators: Vec<Vec<X>>,
        d_degree: isize,
        d_matrices: Vec<CsMat<R>>,
    }

    impl<R> TestChainComplex<R>
    where R: Ring + CsMatElem + From<i32>, for<'x> &'x R: RingOps<R> {
        pub fn new(d_degree: isize, d_matrices: Vec<CsMat<i32>>) -> TestChainComplex<R> {
            assert!(d_degree == 1 || d_degree == -1);
    
            let mut count = 0;
            let mut gens = |n: usize| { 
                let g = (count .. count + n).into_iter().collect();
                count += n;
                g
            };

            // insert matrix so that `rank(C_k) = d_k.cols()`.

            let mut d_matrices = d_matrices;
            
            if d_matrices.is_empty() { 
                d_matrices.push(CsMat::zero((0, 0)))
            } else if d_degree > 0 { 
                let d_n_1 = d_matrices.last().unwrap();
                let k = d_n_1.rows();
                let d_n = CsMat::zero((0, k));
                d_matrices.push(d_n);
            } else if d_degree < 0 { 
                let d_1 = d_matrices.first().unwrap();
                let k = d_1.rows();
                let d_0 = CsMat::zero((0, k));
                d_matrices.insert(0, d_0);
            }

            let generators: Vec<_> = d_matrices.iter().map(|d| {
                gens( d.cols() )
            }).collect();
    
            let d_matrices = d_matrices.into_iter().map(|d| 
                d.map(|&a| R::from(a))
            ).collect();
    
            let n = generators.len() as isize;
            let range = 0..=(n - 1);
            
            Self { range, generators, d_matrices, d_degree }
        }
    }

    impl<R> Graded for TestChainComplex<R>
    where 
        R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> 
    {
        type Index = isize;
        type IndexRange = Either<RangeInclusive<isize>, Rev<RangeInclusive<isize>>>;
        
        fn in_range(&self, k: isize) -> bool {
            self.range.contains(&k)
        }

        fn range(&self) -> Self::IndexRange {
            let range = self.range.clone();
            if self.d_degree > 0 {
                Either::Left(range)
            } else {
                Either::Right(range.rev())
            }
        }
    }

    impl<R> ChainComplex for TestChainComplex<R> 
    where 
        R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> 
    { 
        type R = R;
        type Generator = X;

        fn d_degree(&self) -> isize {
            self.d_degree
        }

        fn rank(&self, k: isize) -> usize {
            if self.in_range(k) { 
                let k = k as usize;
                self.generators[k].len()
            } else {
                0
            }
        }

        fn generators(&self, k: isize) -> Vec<&X> {
            if self.in_range(k) { 
                let k = k as usize;
                self.generators[k].iter().collect()
            } else {
                vec![]
            }
        }

        fn d_matrix(&self, k: isize) -> CsMat<R> {
            if self.in_range(k) { 
                let k = k as usize;
                self.d_matrices[k].clone()
            } else {
                let m = self.rank(k + self.d_degree);
                let n = self.rank(k);
                CsMat::zero((m, n))
            }
        }

        fn differentiate(&self, k: isize, x: &X) -> Vec<(X, R)> {
            self.impl_differentiate_from_d_matrix(k, x)
        }
    }

    impl<R> TestChainComplex<R>
    where R: Ring + CsMatElem + From<i32>, for<'x> &'x R: RingOps<R> {
        pub fn zero() -> Self {
            Self::new(-1, vec![])
        }
    
        pub fn d3() -> Self {
            Self::new(
                -1,
                vec![
                    CsMat::csc_from_vec((4, 6), vec![-1, -1, 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 1]),
                    CsMat::csc_from_vec((6, 4), vec![1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1] ),
                    CsMat::csc_from_vec((4, 1), vec![-1, 1, -1, 1])
                ]
            )
        }
    
        pub fn s2() -> Self {
            Self::new(
                -1,
                vec![
                    CsMat::csc_from_vec((4, 6), vec![-1, -1, 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 1]),
                    CsMat::csc_from_vec((6, 4), vec![1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1] )
                ]
            )
        }
    
        pub fn t2() -> Self {
            Self::new(
                -1,
                vec![
                    CsMat::csc_from_vec((9, 27), vec![-1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, -1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1]),
                    CsMat::csc_from_vec((27, 18), vec![1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1])
                ]
            )
        }
    
        pub fn rp2() -> Self { 
            Self::new(
                -1,
                vec![
                    CsMat::csc_from_vec((6, 15), vec![-1, -1, 0, 0, 0, 0, 0, -1, -1, 0, -1, 0, 0, 0, 0, 1, 0, -1, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 1, 1, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1] ),
                    CsMat::csc_from_vec((15, 10), vec![1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1] ),
                ]
            )
        }
    }
}