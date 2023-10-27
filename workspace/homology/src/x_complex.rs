use std::ops::Index;

use bimap::BiHashMap;
use delegate::delegate;

use itertools::Itertools;
use num_traits::Zero;
use yui_core::{Ring, RingOps, EucRing, EucRingOps, Deg, isize2, isize3};
use yui_lin_comb::{Gen, LinComb};
use yui_matrix::sparse::{SpMat, SpVec};

use crate::{GridBase, GridIter, ChainComplexSummand, ChainComplexSummandTrait};

use super::grid::GridTrait;
use super::complex::{ChainComplexTrait, ChainComplexBase};
use super::homology::HomologyBase;

pub type XChainComplex <X, R> = XChainComplexBase<isize,  X, R>;
pub type XChainComplex2<X, R> = XChainComplexBase<isize2, X, R>;
pub type XChainComplex3<X, R> = XChainComplexBase<isize3, X, R>;

pub struct XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    inner: ChainComplexBase<I, R>,
    gens: GridBase<I, BiHashMap<usize, X>>
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
        let gens = GridBase::new(
            support.clone().into_iter(), 
            |i| gens_map(i).into_iter().enumerate().collect()
        );

        let inner = ChainComplexBase::new(
            support.into_iter(), 
            d_deg, 
            |i| Self::make_matrix(i, gens.get(i), gens.get(i + d_deg), &d_map)
        );

        Self { inner, gens }
    }

    fn make_matrix<F>(i: I, from: &BiHashMap<usize, X>, to: &BiHashMap<usize, X>, d: &F) -> SpMat<R>
    where F: Fn(I, &X) -> Vec<(X, R)> {
        let (m, n) = (to.len(), from.len());

        SpMat::generate((m, n), |set|
            for (&j, x) in from.iter() {
                let ys = d(i, x);
                for (y, a) in ys {
                    let &i = to.get_by_right(&y).unwrap();
                    set(i, j, a);
                }
            }
        )
    }

    pub fn gens(&self, i: I) -> impl Iterator<Item = &X> { 
        let gens = self.gens.get(i);
        let r = gens.len();

        (0..r).map(|j| 
            gens.get_by_left(&j).unwrap()
        )
    }

    pub fn gen_at(&self, i: I, j: usize) -> &X {
        self.gens.get(i).get_by_left(&j).unwrap()
    }

    pub fn index_of(&self, i: I, x: &X) -> usize {
        self.gens.get(i).get_by_right(x).unwrap().clone()
    }

    pub fn vectorize(&self, i: I, z: &LinComb<X, R>) -> SpVec<R> {
        let r = self[i].rank();
        SpVec::generate(r, |set| { 
            for (x, a) in z.iter() { 
                let j = self.index_of(i, x);
                set(j, a.clone());
            }
        })
    }

    pub fn as_chain(&self, i: I, v: &SpVec<R>) -> LinComb<X, R> {
        if self.is_supported(i) { 
            let elems = v.iter().map(|(j, a)| 
                (self.gen_at(i, j).clone(), a.clone())
            );
            LinComb::from_iter(elems)
        } else {
            LinComb::zero()
        }
    }

    pub fn d(&self, i: I, z: &LinComb<X, R>) -> LinComb<X, R> { 
        let d = self[i].d_matrix();
        let v = self.vectorize(i, z);
        let w = d * v;
        self.as_chain(i + self.d_deg(), &w)
    }

    pub fn inner(&self) -> &ChainComplexBase<I, R> { 
        &self.inner
    }
}

impl<I, X, R> Index<I> for XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    type Output = ChainComplexSummand<R>;

    delegate! { 
        to self.inner { 
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
    type E = ChainComplexSummand<R>;
    
    delegate! { 
        to self.inner { 
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

    delegate! { 
        to self.inner { 
            fn d_deg(&self) -> I;
        }
    }
}

impl<I, X, R> XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: EucRing, for<'x> &'x R: EucRingOps<R>,
{
    delegate! { 
        to self.inner { 
            pub fn homology(&self, with_trans: bool) -> HomologyBase<I, R>;
        }
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
    }
}