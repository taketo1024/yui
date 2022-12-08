use std::ops::RangeInclusive;
use std::collections::HashMap;
use std::hash::Hash;
use sprs::{CsMat, CsVec, TriMat};
use num_traits::One;
use crate::math::traits::{Ring, RingOps};
use crate::math::matrix::sparse::*;
use crate::utils::hashmap;

pub trait ChainGenerator: Clone + PartialEq + Eq + Hash {}

impl<T> ChainGenerator for T 
where T: Clone + PartialEq + Eq + Hash {}

pub trait ChainComplex 
where Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> {
    type R;
    type Generator: ChainGenerator;

    // -- must override -- //
    fn range(&self) -> RangeInclusive<isize>;
    fn generators(&self, k: isize) -> Vec<&Self::Generator>;
    fn d_degree(&self) -> isize;
    fn differentiate(&self, k: isize, x:&Self::Generator) -> Vec<(Self::Generator, Self::R)>;

    // -- convenient methods -- //
    fn rank(&self, k: isize) -> usize { 
        if self.range().contains(&k) { 
            self.generators(k).len()
        } else {
            0
        }
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
}

pub trait ChainComplexSparseD: ChainComplex
where Self::R: Ring + CsMatElem, for<'x> &'x Self::R: RingOps<Self::R> {
    fn d_matrix(&self, k: isize) -> CsMat<Self::R> {
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

    fn check_d_at(&self, k: isize) { 
        let d1 = self.d_matrix(k);
        let d2 = self.d_matrix(k + self.d_degree());
        let res = &d2 * &d1;
        assert!( res.is_zero() );
    }

    fn check_d_all(&self) {
        for k in self.range() { 
            let k2 = k + self.d_degree();
            if !self.range().contains(&k2) { continue }

            self.check_d_at(k);
        }
    }
}