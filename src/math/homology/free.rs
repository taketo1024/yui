use std::fmt::Display;
use std::collections::HashMap;
use crate::math::matrix::sp_vec::SpVec;
use crate::math::matrix::sparse::SpMat;
use crate::math::types::lin_comb::{FreeGen, LinComb};
use yui_core::{Ring, RingOps};
use super::base::RModStr;

pub struct FreeRModStr<X, R>
where 
    X: FreeGen,
    R: Ring, for<'x> &'x R: RingOps<R>, 
{
    generators: Vec<X>,
    tors: Vec<R>, // always empty
}

impl<X, R> FreeRModStr<X, R>
where 
    X: FreeGen,
    R: Ring, for<'x> &'x R: RingOps<R>, 
{
    pub fn new(generators: Vec<X>) -> Self {
        Self { generators, tors: vec![] }
    }

    pub fn generators(&self) -> &Vec<X> {
        &self.generators
    }

    pub fn vectorize_x(&self, x: X) -> SpVec<R> {
        self.vectorize(&LinComb::wrap(x))
    }

    pub fn vectorize(&self, z: &LinComb<X, R>) -> SpVec<R> {
        let n = self.generators.len();

        SpVec::generate(n, |set| { 
            for (x, a) in z.iter() { 
                // MEMO non-effective
                let Some(i) = self.generators.iter().position(|z| x == z) else { 
                    continue
                };
                set(i, a.clone())
            }
        })
    }

    pub fn make_matrix<F>(&self, target: &Self, f: F) -> SpMat<R>
    where 
        F: Fn(&X) -> Vec<(X, R)> 
    {
        let (m, n) = (target.generators.len(), self.generators.len());

        let t_ind = target.generators.iter()
                .enumerate()
                .map(|(i, y)| (y, i))
                .collect::<HashMap<_, _>>();

        SpMat::generate((m, n), |set|
            for (j, x) in self.generators.iter().enumerate() {
                let ys = f(x);
                for (y, a) in ys {
                    if !t_ind.contains_key(&y) { continue }
                    let i = t_ind[&y];
                    set(i, j, a);
                }
            }
        )
    }
}

impl<X, R> RModStr for FreeRModStr<X, R>
where 
    X: FreeGen,
    R: Ring, for<'x> &'x R: RingOps<R>, 
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

impl<X, R> Display for FreeRModStr<X, R>
where 
    X: FreeGen,
    R: Ring, for<'x> &'x R: RingOps<R>, 
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt_default(f)
    }
}