use std::collections::HashMap;
use std::ops::Index;

use delegate::delegate;
use itertools::Itertools;

use yui_core::{Ring, RingOps, EucRing, EucRingOps, Deg, isize2, isize3, IndexList};
use yui_lin_comb::{Gen, LinComb};
use yui_matrix::sparse::{SpMat, SpVec};

use crate::{GridBase, GridIter, DisplayForGrid, ChainComplexDisplay, ChainComplexSummand};

use super::grid::GridTrait;
use super::complex::ChainComplexTrait;
use super::homology::HomologyBase;

pub type XChainComplex <X, R> = XChainComplexBase<isize,  X, R>;
pub type XChainComplex2<X, R> = XChainComplexBase<isize2, X, R>;
pub type XChainComplex3<X, R> = XChainComplexBase<isize3, X, R>;

#[derive(Default)]
pub struct XChainComplexSummand<X, R>
where 
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    gens: IndexList<X>,
    inner: ChainComplexSummand<R>
}

impl<X, R> XChainComplexSummand<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(gens: IndexList<X>, d_matrix: SpMat<R>) -> Self { 
        let inner = ChainComplexSummand::new(d_matrix);
        Self { gens, inner }
    }

    pub fn rank(&self) -> usize { 
        self.gens.len()
    }

    pub fn gens(&self) -> impl Iterator<Item = &X> { 
        self.gens.iter()
    }

    pub fn gen(&self, i: usize) -> &X {
        &self.gens[i]
    }

    pub fn index_of(&self, x: &X) -> usize {
        self.gens.index_of(x).unwrap()
    }

    pub fn vectorize(&self, z: &LinComb<X, R>) -> SpVec<R> {
        let n = self.rank();
        SpVec::generate(n, |set| { 
            for (x, a) in z.iter() { 
                let i = self.index_of(x);
                set(i, a.clone());
            }
        })
    }

    pub fn as_chain(&self, v: &SpVec<R>) -> LinComb<X, R> {
        assert_eq!(v.dim(), self.rank());

        let elems = v.iter().map(|(i, a)| 
            (self.gen(i).clone(), a.clone())
        );
        LinComb::from_iter(elems)
    }

    delegate! { 
        to self.inner { 
            pub fn d_matrix(&self) -> &SpMat<R>;
            pub fn d(&self, v: &SpVec<R>) -> SpVec<R>;
            pub fn is_cycle(&self, v: &SpVec<R>) -> bool;
            pub fn module_str(&self) -> String;
        }
    }
}

impl<X, R> DisplayForGrid for XChainComplexSummand<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn display_for_grid(&self) -> String {
        self.module_str()
    }
}

pub struct XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    d_deg: I,
    summands: GridBase<I, XChainComplexSummand<X, R>>
}

impl<I, X, R> XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    pub fn new<It, F1, F2>(support: It, d_deg: I, gens_map: F1, d_map: F2) -> Self
    where 
        It: Iterator<Item = I>, 
        F1: Fn(I) -> Vec<X>,
        F2: Fn(I, &X) -> Vec<(X, R)>
    {
        let support = support.collect_vec();

        let mut gens = support.iter().map(|&i| {
            (i, IndexList::new(gens_map(i).into_iter()))
        }).collect::<HashMap<_, _>>();

        for &i in support.iter() {
            let r = gens[&i].len();
            if r > 0 && !gens.contains_key(&(i - d_deg)) { 
                gens.insert(i - d_deg, IndexList::default());
            }
        }

        let mut mats = gens.keys().map( |&i| {
            let from = &gens[&i];
            let to = gens.get(&(i + d_deg)); // optional
            let d = Self::make_matrix(i, from, to, &d_map);
            (i, d)
        }).collect::<HashMap<_, _>>();

        let data = gens.into_iter().map(|(i, gens)| {
            let d_matrix = mats.remove(&i).unwrap();
            let summand = XChainComplexSummand::new(gens, d_matrix);
            (i, summand)
        }).collect();

        let summands = GridBase::new_raw(support, data);

        Self { d_deg, summands }
    }

    fn make_matrix<F>(i: I, from: &IndexList<X>, to_opt: Option<&IndexList<X>>, d: &F) -> SpMat<R>
    where F: Fn(I, &X) -> Vec<(X, R)> {
        let Some(to) = to_opt else { 
            return SpMat::zero((0, from.len()))
        };
        let (m, n) = (to.len(), from.len());

        SpMat::generate((m, n), |set|
            for (j, x) in from.iter().enumerate() {
                let ys = d(i, x);
                for (y, a) in ys {
                    let i = to.index_of(&y).unwrap();
                    set(i, j, a);
                }
            }
        )
    }

    pub fn d(&self, i: I, z: &LinComb<X, R>) -> LinComb<X, R> { 
        let d = self[i].d_matrix();
        let v = self[i].vectorize(z);
        let w = d * v;

        let i1 = i + self.d_deg;
        self[i1].as_chain(&w)
    }
}

impl<I, X, R> Index<I> for XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    type Output = XChainComplexSummand<X, R>;

    delegate! { 
        to self.summands { 
            fn index(&self, index: I) -> &Self::Output;
        }
    }
}

impl<I, X, R> GridTrait<I> for XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    type Itr = GridIter<I>;
    type E = XChainComplexSummand<X, R>;
    
    delegate! { 
        to self.summands { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: I) -> bool;
            fn get(&self, i: I) -> &Self::E;
        }
    }
}

impl<I, X, R> ChainComplexTrait<I> for XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    type R = R;
    fn d_deg(&self) -> I { 
        self.d_deg
    }

    fn rank(&self, i: I) -> usize {
        self[i].rank()
    }

    fn d_matrix(&self, i: I) -> &SpMat<Self::R> {
        self[i].d_matrix()
    }
}

impl<I, X, R> ChainComplexDisplay<I> for XChainComplexBase<I, X, R>
where I: Deg, X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {}

impl<I, X, R> XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: EucRing, for<'x> &'x R: EucRingOps<R>,
{
    pub fn homology(&self, with_trans: bool) -> HomologyBase<I, R> { 
        HomologyBase::compute_from(self, with_trans)
    }
}

#[cfg(test)]
mod tests { 
    use yui_lin_comb::Free;

    use super::*;

    type X = Free<i32>;
    fn e(i: isize) -> X { 
        X::from(i as i32)
    }

    #[test]
    fn test() { 
        let c = XChainComplex::<X, i32>::new(0..2, -1, 
            |i| vec![e(i)], 
            |i, x| { 
                if x.0 == 1 { 
                    vec![(e(i - 1), 1)] 
                } else {
                    vec![]
                }
            }
        );

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 1);
        assert_eq!(c[1].d_matrix(), &SpMat::from_vec((1,1), vec![1]));

        let h = c.homology(false);
        
        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[1].rank(), 0);
    }
}