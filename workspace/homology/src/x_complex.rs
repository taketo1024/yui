use std::marker::PhantomData;
use std::ops::Index;

use delegate::delegate;

use yui_core::{Ring, RingOps, EucRing, EucRingOps, Deg, isize2, isize3, IndexList};
use yui_lin_comb::{Gen, LinComb};
use yui_matrix::sparse::{SpMat, SpVec};

use crate::{GridBase, GridIter, DisplayForGrid, ChainComplexDisplay, ChainComplexBase, RModStr};

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
    _r: PhantomData<R>
}

impl<X, R> XChainComplexSummand<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(gens: Vec<X>) -> Self { 
        let gens = IndexList::new(gens.into_iter());
        Self { gens, _r: PhantomData }
    }

    pub fn gens(&self) -> &IndexList<X> { 
        &self.gens
    }

    pub fn gen(&self, i: usize) -> &X {
        &self.gens[i]
    }

    pub fn vectorize(&self, z: &LinComb<X, R>) -> SpVec<R> {
        let n = self.rank();
        SpVec::generate(n, |set| { 
            for (x, a) in z.iter() { 
                let i = self.gens.index_of(x).unwrap();
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

    pub fn make_matrix<Y, F>(&self, to: &XChainComplexSummand<Y, R>, f: F) -> SpMat<R>
    where Y: Gen, F: Fn(&X) -> Vec<(Y, R)> {
        make_matrix(&self.gens, &to.gens, f)
    }
}

impl<X, R> RModStr for XChainComplexSummand<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn rank(&self) -> usize {
        self.gens.len()
    }

    fn tors(&self) -> Vec<&Self::R> {
        vec![]
    }
}

impl<X, R> DisplayForGrid for XChainComplexSummand<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn display_for_grid(&self) -> String {
        self.math_symbol()
    }
}

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
            XChainComplexSummand::new(gens_map(i))
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

fn make_matrix<X, Y, R, F>(from: &IndexList<X>, to: &IndexList<Y>, f: F) -> SpMat<R>
where 
    X: Gen, Y: Gen, 
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(&X) -> Vec<(Y, R)> 
{
    let (m, n) = (to.len(), from.len());
    SpMat::generate((m, n), |set|
        for (j, x) in from.iter().enumerate() {
            let ys = f(x);
            for (y, a) in ys {
                let i = to.index_of(&y).unwrap();
                set(i, j, a);
            }
        }
    )
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