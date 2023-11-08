use std::ops::Index;

use delegate::delegate;

use yui_core::{Ring, RingOps, Deg, isize2, isize3};
use yui_lin_comb::{Gen, LinComb};
use yui_matrix::sparse::SpMat;

use crate::{GridBase, GridIter, ChainComplexDisplay, XModStr, RModStr};

use super::grid::GridTrait;
use super::complex::ChainComplexTrait;

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
    d_deg: I,
    d_map: Box<dyn Fn(I, &X) -> Vec<(X, R)>>,
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
        F2: Fn(I, &X) -> Vec<(X, R)> + 'static
    {
        let summands = GridBase::new(support, |i| 
            XChainComplexSummand::free(gens_map(i))
        );

        let d_map = Box::new(d_map);
        Self { summands, d_deg, d_map }
    }

    pub fn d(&self, i: I, z: &LinComb<X, R>) -> LinComb<X, R> { 
        z.apply(|x| (self.d_map)(i, x))
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

    fn d_matrix(&self, i: I) -> SpMat<Self::R> { 
        let i1 = i + self.d_deg;
        let (m, n) = (self[i1].rank(), self[i].rank());

        SpMat::generate((m, n), |set| { 
            for j in 0..n { 
                let z = self[i].gen_chain(j);
                let dz = self.d(i, &z);
                let w = self[i1].vectorize(&dz);
                
                for (i, a) in w.iter() {
                    set(i, j, a.clone())
                }
            }
        })
    }
}

impl<I, X, R> ChainComplexDisplay<I> for XChainComplexBase<I, X, R>
where I: Deg, X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {}

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
pub(crate) mod tests { 
    use itertools::Itertools;
    use yui_lin_comb::Free;
    use yui_matrix::sparse::SpVec;

    use crate::{RModStr, ChainComplex, Grid};

    use super::*;

    type X = Free<i64>;
    fn e(i: usize) -> X { 
        X::from(i as i64)
    }

    impl From<ChainComplex<i64>> for XChainComplex<X, i64> {
        fn from(c: ChainComplex<i64>) -> Self {
            let gens = Grid::new(c.support(), |i| { 
                let n = c[i].rank();
                let gens = (0..n).map(|j| e(j)).collect_vec();
                gens
            });

            Self::new(
                c.support(), c.d_deg(), 
                move |i| gens[i].clone(), 
                move |i, x| {
                    let n = c[i].rank();
                    let v = SpVec::unit(n, x.0 as usize);
                    let dv = c.d(i, &v);
                    dv.iter().map(|(i, a)| (e(i), a.clone())).collect()
                }
            )
        }
    }

    #[test]
    fn d3() { 
        let c = XChainComplex::from(ChainComplex::d3());

        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 4);
        assert_eq!(c[3].rank(), 1);

        c.check_d_all();
    }

    #[test]
    fn s2() { 
        let c = XChainComplex::from(ChainComplex::s2());

        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 4);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();
    }

    #[test]
    fn t2() { 
        let c = XChainComplex::from(ChainComplex::t2());

        assert_eq!(c[0].rank(), 9);
        assert_eq!(c[1].rank(), 27);
        assert_eq!(c[2].rank(), 18);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();
    }

    #[test]
    fn rp2() { 
        let c = XChainComplex::from(ChainComplex::rp2());
        
        assert_eq!(c[0].rank(), 6);
        assert_eq!(c[1].rank(), 15);
        assert_eq!(c[2].rank(), 10);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();
    }
}