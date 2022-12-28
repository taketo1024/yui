use std::fmt::Display;
use std::hash::Hash;
use std::collections::HashMap;
use sprs::{CsVec, CsMat, TriMat};
use crate::utils::collections::hashmap;
use crate::math::traits::{Ring, RingOps};
use super::base::RModStr;

pub trait FreeGenerator: Clone + PartialEq + Eq + Hash {}
impl<T> FreeGenerator for T 
where T: Clone + PartialEq + Eq + Hash {}

pub struct FreeRModStr<R, X>
where 
    R: Ring, for<'x> &'x R: RingOps<R>, 
    X: FreeGenerator
{
    generators: Vec<X>,
    tors: Vec<R> // always empty
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

    pub fn vectorize(&self, x: &X) -> CsVec<R> {
        vectorize(self.generators(), x)
    }

    pub fn clone_filter<F>(&self, pred: F) -> Self
    where F: Fn(&X) -> bool { 
        let gens = self.generators.iter().filter(|x| pred(x)).cloned().collect();
        Self::new(gens)
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

pub fn vectorize<R, X>(generators: &Vec<X>, x:&X) -> CsVec<R>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    X: FreeGenerator
{
    vectorize_comb(generators, hashmap![x => R::one()])
}

pub fn vectorize_comb<R, X>(generators: &Vec<X>, z:HashMap<&X, R>) -> CsVec<R>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    X: FreeGenerator
{
    let n = generators.len();

    let mut v_ind: Vec<usize> = vec![];
    let mut v_val: Vec<R> = vec![];

    for (x, a) in z { 
        let Some(i) = generators.iter().position(|z| x == z) else { 
            continue 
        };
        v_ind.push(i);
        v_val.push(a);
    }

    CsVec::new(n, v_ind, v_val)
}

pub fn make_matrix<R, X, F>(source: &Vec<X>, target: &Vec<X>, f: F) -> CsMat<R>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    X: FreeGenerator,
    F: Fn(&X) -> Vec<(X, R)> 
{
    let (m, n) = (target.len(), source.len());

    let t_ind = target.iter()
            .enumerate()
            .map(|(i, y)| (y, i))
            .collect::<HashMap<_, _>>();

    let mut trip = TriMat::new((m, n));

    for (j, x) in source.iter().enumerate() {
        let ys = f(x);
        for (y, a) in ys {
            if !t_ind.contains_key(&y) { continue }
            let i = t_ind[&y];
            trip.add_triplet(i, j, a);
        }
    }

    trip.to_csc()
}