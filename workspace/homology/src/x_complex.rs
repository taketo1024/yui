use std::ops::Index;

use delegate::delegate;

use yui_core::{Ring, RingOps, EucRing, EucRingOps, Deg, isize2, isize3};
use yui_lin_comb::{Gen, LinComb};
use yui_matrix::sparse::SpMat;

use crate::{GridBase, GridIter, ChainComplexDisplay, ChainComplexBase, XModStr};

use super::grid::GridTrait;
use super::complex::ChainComplexTrait;
use super::homology::HomologyBase;

pub type XChainComplex <X, R> = XChainComplexBase<isize,  X, R>;
pub type XChainComplex2<X, R> = XChainComplexBase<isize2, X, R>;
pub type XChainComplex3<X, R> = XChainComplexBase<isize3, X, R>;

pub type XChainComplexSummand<X, R> = XModStr<X, R>;

pub struct XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    summands: GridBase<I, XChainComplexSummand<X, R>>,
    inner: ChainComplexBase<I, R>
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
        let summands = GridBase::new(support, |i| 
            XChainComplexSummand::from_iter(gens_map(i))
        );

        let inner = ChainComplexBase::new(summands.support(), d_deg, |i| {
            let from = &summands[i];
            let to = &summands[i + d_deg];
            from.make_matrix(to, |x| d_map(i, x))
        });

        Self { summands, inner }
    }

    pub fn d(&self, i: I, z: &LinComb<X, R>) -> LinComb<X, R> { 
        let d = self.d_matrix(i);
        let v = self[i].vectorize(z);
        let w = d * v;

        let i1 = i + self.d_deg();
        self[i1].as_chain(&w)
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

    delegate! { 
        to self.inner { 
            fn d_deg(&self) -> I;
            fn d_matrix(&self, i: I) -> &SpMat<Self::R>;
        }
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

impl<I, X, R> Index<I> for XChainComplexBase<I, X, R>
where I: Deg, X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = XChainComplexSummand<X, R>;
    fn index(&self, i: I) -> &Self::Output {
        self.get(i)
    }
}

impl<X, R> Index<(isize, isize)> for XChainComplex2<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = XChainComplexSummand<X, R>;
    fn index(&self, i: (isize, isize)) -> &Self::Output {
        self.get(i.into())
    }
}

impl<X, R> Index<(isize, isize, isize)> for XChainComplex3<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = XChainComplexSummand<X, R>;
    fn index(&self, i: (isize, isize, isize)) -> &Self::Output {
        self.get(i.into())
    }
}

#[cfg(test)]
mod tests { 
    use yui_lin_comb::Free;

    use crate::RModStr;

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
        assert_eq!(c.d_matrix(1), &SpMat::from_vec((1,1), vec![1]));

        let h = c.homology(false);
        
        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[1].rank(), 0);
    }
}