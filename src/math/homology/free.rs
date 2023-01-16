use std::fmt::Display;
use std::collections::HashMap;
use sprs::{CsVec, CsMat, TriMat};
use crate::math::free::{FreeGenerator, LinComb};
use crate::math::traits::{Ring, RingOps};
use super::base::RModStr;

pub struct FreeRModStr<R, X>
where 
    R: Ring, for<'x> &'x R: RingOps<R>, 
    X: FreeGenerator
{
    generators: Vec<X>,
    tors: Vec<R>, // always empty
}

impl<R, X> FreeRModStr<R, X>
where 
    R: Ring, for<'x> &'x R: RingOps<R>, 
    X: FreeGenerator
{
    pub fn new(generators: Vec<X>) -> Self {
        Self { generators, tors: vec![] }
    }

    pub fn generators(&self) -> &Vec<X> {
        &self.generators
    }

    pub fn vectorize_x(&self, x: X) -> CsVec<R> {
        self.vectorize(&LinComb::wrap(x))
    }

    pub fn vectorize(&self, z: &LinComb<X, R>) -> CsVec<R> {
        let n = self.generators.len();

        let mut v_ind: Vec<usize> = vec![];
        let mut v_val: Vec<R> = vec![];
    
        for (x, a) in z.iter() { 
            // MEMO non-effective
            let Some(i) = self.generators.iter().position(|z| x == z) else { 
                continue
            };
            v_ind.push(i);
            v_val.push(a.clone());
        }
    
        CsVec::new(n, v_ind, v_val)
    }

    pub fn make_matrix<F>(&self, target: &Self, f: F) -> CsMat<R>
    where 
        F: Fn(&X) -> Vec<(X, R)> 
    {
        let (m, n) = (target.generators.len(), self.generators.len());

        let t_ind = target.generators.iter()
                .enumerate()
                .map(|(i, y)| (y, i))
                .collect::<HashMap<_, _>>();

        let mut trip = TriMat::new((m, n));

        for (j, x) in self.generators.iter().enumerate() {
            let ys = f(x);
            for (y, a) in ys {
                if !t_ind.contains_key(&y) { continue }
                let i = t_ind[&y];
                trip.add_triplet(i, j, a);
            }
        }

        trip.to_csc()
    }
}

impl<R, X> RModStr for FreeRModStr<R, X>
where 
    R: Ring, for<'x> &'x R: RingOps<R>, 
    X: FreeGenerator
{
    type R = R;

    fn zero() -> Self {
        Self::new(vec![])
    }

    fn rank(&self) -> usize {
        self.generators.len()
    }

    fn tors(&self) -> &Vec<Self::R> {
        &self.tors
    }
}

impl<R, X> Display for FreeRModStr<R, X>
where 
    R: Ring, for<'x> &'x R: RingOps<R>, 
    X: FreeGenerator
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt_default(f)
    }
}