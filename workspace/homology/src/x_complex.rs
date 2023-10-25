use std::collections::HashMap;
use delegate::delegate;

use itertools::Itertools;
use yui_core::{Ring, RingOps, EucRing, EucRingOps, Deg, isize2, isize3};
use yui_lin_comb::{Gen, LinComb};
use yui_matrix::sparse::{SpMat, SpVec};

use crate::HomologySummand;

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
    gens: HashMap<I, Vec<X>>,
    empty_gens: Vec<X>
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
        let gens = support.iter().map(|&i| 
            (i, gens(i))
        ).collect::<HashMap<_, _>>();
        let mut mats = support.iter().map(|&i|
            (i, Self::make_matrix(i, &gens[&i], gens.get(&(i + d_deg)).unwrap_or(&vec![]), &d_map))
        ).collect::<HashMap<_, _>>();
        let inner = ChainComplexBase::new(support.into_iter(), d_deg, move |i| { 
            std::mem::take(mats.get_mut(&i).unwrap())
        });
        Self { 
            inner,
            gens,
            empty_gens: vec![]
        }
    }

    fn make_matrix<F>(i: I, from: &Vec<X>, to: &Vec<X>, d: &F) -> SpMat<R>
    where F: Fn(I, &X) -> Vec<(X, R)> {
        let (m, n) = (to.len(), from.len());
        let dict = to.iter()
                .enumerate()
                .map(|(i, y)| (y, i))
                .collect::<HashMap<_, _>>();

        SpMat::generate((m, n), |set|
            for (j, x) in from.iter().enumerate() {
                let ys = d(i, x);
                for (y, a) in ys {
                    let Some(&i) = dict.get(&y) else { continue };
                    set(i, j, a);
                }
            }
        )
    }

    pub fn gens(&self, i: I) -> &Vec<X> { 
        if self.is_supported(i) {
            &self.gens[&i]
        } else { 
            &self.empty_gens
        }
    }

    pub fn vectorize_x(&self, i: I, x: &X) -> SpVec<R> {
        let r = self.rank(i);
        let gens = self.gens(i);

        if let Some(j) = gens.iter().position(|y| x == y) {
            SpVec::unit(r, j)
        } else { 
            SpVec::zero(r)
        }
    }

    pub fn vectorize(&self, i: I, z: &LinComb<X, R>) -> SpVec<R> {
        let r = self.rank(i);
        let gens = self.gens(i);

        SpVec::generate(r, |set| { 
            for (x, a) in z.iter() { 
                // MEMO non-effective
                let Some(j) = gens.iter().position(|y| x == y) else {
                    continue
                };
                set(j, a.clone())
            }
        })
    }

    pub fn as_chain(&self, i: I, v: &SpVec<R>) -> LinComb<X, R> {
        let gens = self.gens(i);
        let elems = v.iter().map(|(j, a)| 
            (gens[j].clone(), a.clone())
        );
        LinComb::from_iter(elems)
    }

    pub fn d_of(&self, i: I, x: &X) -> LinComb<X, R> { 
        self.d(i, &LinComb::from_gen(x.clone()))
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
            fn display(&self, i: I) -> String;
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
        self.gens(i).len()
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