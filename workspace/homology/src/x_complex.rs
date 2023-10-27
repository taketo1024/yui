use std::collections::HashMap;
use bimap::BiHashMap;
use delegate::delegate;

use itertools::Itertools;
use num_traits::Zero;
use yui_core::{Ring, RingOps, EucRing, EucRingOps, Deg, isize2, isize3};
use yui_lin_comb::{Gen, LinComb};
use yui_matrix::sparse::{SpMat, SpVec};

use crate::{HomologySummand, DisplayAt};

use super::graded::Graded;
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
    gens: HashMap<I, BiHashMap<usize, X>>,
    empty_gens: BiHashMap<usize, X>
}

impl<I, X, R> XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    pub fn new<It, F1, F2>(support: It, d_deg: I, gens: F1, d_map: F2) -> Self
    where 
        It: Iterator<Item = I>, 
        F1: Fn(I) -> Vec<X>,
        F2: Fn(I, &X) -> Vec<(X, R)>
    {
        let support = support.collect_vec();
        let gens: HashMap<_, _> = support.iter().map(|&i| {
            let g: BiHashMap<_, _> = gens(i).into_iter().enumerate().collect();
            (i, g)
        }).collect();
        let empty_gens = BiHashMap::new();

        let mut mats: HashMap<_, _> = support.iter().map(|&i| {
            let from = &gens[&i];
            let to = gens.get(&(i + d_deg));
            let d = Self::make_matrix(i, from, to, &d_map);
            (i, d)
        }).collect();

        let inner = ChainComplexBase::new(support.into_iter(), d_deg, move |i| { 
            std::mem::take(mats.get_mut(&i).unwrap())
        });

        Self { inner, gens, empty_gens }
    }

    fn make_matrix<F>(i: I, from: &BiHashMap<usize, X>, to: Option<&BiHashMap<usize, X>>, d: &F) -> SpMat<R>
    where F: Fn(I, &X) -> Vec<(X, R)> {
        let Some(to) = to else { 
            return SpMat::zero((0, from.len()))
        };
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
        let gens = self.gens.get(&i).unwrap_or(&self.empty_gens);
        let r = gens.len();

        (0..r).map(|j| 
            gens.get_by_left(&j).unwrap()
        )
    }

    pub fn gen_at(&self, i: I, j: usize) -> &X {
        self.gens[&i].get_by_left(&j).unwrap()
    }

    pub fn index_of(&self, i: I, x: &X) -> usize {
        self.gens[&i].get_by_right(x).unwrap().clone()
    }

    pub fn vectorize(&self, i: I, z: &LinComb<X, R>) -> SpVec<R> {
        let r = self.rank(i);
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
        let v = self.vectorize(i, z);
        let w = self.d_matrix(i) * v;
        self.as_chain(i + self.d_deg(), &w)
    }

    pub fn inner(&self) -> &ChainComplexBase<I, R> { 
        &self.inner
    }
}

impl<I, X, R> Graded<I> for XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    type Itr = std::vec::IntoIter<I>;
    
    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Itr;
        }
    }

    fn is_supported(&self, i: I) -> bool {
        self.gens.contains_key(&i)
    }
}

impl<I, X, R> DisplayAt<I> for XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    delegate! { 
        to self.inner { 
            fn display_at(&self, i: I) -> Option<String>;
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

    fn rank(&self, i: I) -> usize {
        if let Some(gens) = self.gens.get(&i) { 
            gens.len()
        } else { 
            0
        }
    }

    delegate! { 
        to self.inner { 
            fn d_deg(&self) -> I;
            fn d_matrix(&self, i: I) -> &SpMat<Self::R>;
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
            pub fn homology_at(&self, i: I, with_trans: bool) -> HomologySummand<R>;
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

        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 1);
        assert_eq!(c.d_matrix(1), &SpMat::from_vec((1,1), vec![1]));
    }
}