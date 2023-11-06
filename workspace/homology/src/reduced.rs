use std::{ops::Index, collections::HashMap};

use delegate::delegate;

use itertools::Itertools;
use yui_core::{Deg, Ring, RingOps, EucRing, EucRingOps, isize2, isize3};
use yui_matrix::sparse::{Trans, SpMat, SpVec};

use crate::{ChainComplexBase, GridTrait, ChainComplexTrait, HomologyBase, GridBase, GridIter, ChainComplexSummand, ChainComplexDisplay, RModStr};

pub type ReducedComplex<R>  = ReducedComplexBase<isize,  R>;
pub type ReducedComplex2<R> = ReducedComplexBase<isize2, R>;
pub type ReducedComplex3<R> = ReducedComplexBase<isize3, R>;

pub struct ReducedComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    inner: ChainComplexBase<I, R>,
    trans: GridBase<I, Option<Trans<R>>>
}

impl<I, R> ReducedComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new<Itr, F>(support: Itr, d_deg: I, mut d_matrix_map: F) -> Self
    where 
        Itr: Iterator<Item = I>, 
        F: FnMut(I) -> (SpMat<R>, Option<Trans<R>>)
    {
        let support = support.collect_vec();

        let mut mats = HashMap::new();
        let mut trans = HashMap::new();

        for i in support.iter() { 
            let (d, t) = d_matrix_map(i.clone());
            
            mats.insert(i.clone(), d);
            if let Some(t) = t { 
                trans.insert(i.clone(), t);
            }
        }

        let inner = ChainComplexBase::new(
            support.clone().into_iter(), 
            d_deg, 
            move |i| mats.remove(&i).unwrap()
        );
        let trans = GridBase::new(
            inner.support(), 
            move |i| trans.remove(&i)
        );
        Self { inner, trans }
    }

    pub fn bypass<C>(complex: &C) -> Self
    where 
        C: ChainComplexTrait<I, R = R>,
        C::E: RModStr<R = R>
    {
        Self::new(
            complex.support(), 
            complex.d_deg(), 
            |i| (complex.d_matrix(i).clone(), None)
        )
    }

    pub fn trans(&self, i: I) -> Option<&Trans<R>> {
        self.trans.get(i).as_ref()
    }

    delegate! { 
        to self.inner { 
            pub fn d(&self, i: I, v: &SpVec<R>) -> SpVec<R>;
        }
    }
}

impl<I, R> GridTrait<I> for ReducedComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type E = ChainComplexSummand<R>;
    type Itr = GridIter<I>;

    delegate! { 
        to self.inner {
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: I) -> bool;
            fn get(&self, i: I) -> &Self::E;
        }
    }
}

impl<I, R> ChainComplexTrait<I> for ReducedComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    delegate! { 
        to self.inner { 
            fn d_deg(&self) -> I;
            fn d_matrix(&self, i: I) -> &SpMat<Self::R>;
        }
    }
}

impl<I, R> ChainComplexDisplay<I> for ReducedComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {}

impl<I, R> ReducedComplexBase<I, R>
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    delegate! { 
        to self.inner { 
            pub fn homology(&self, with_trans: bool) -> HomologyBase<I, R>;
        }
    }
}

impl<I, R> Index<I> for ReducedComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = ChainComplexSummand<R>;
    fn index(&self, i: I) -> &Self::Output {
        self.get(i)
    }
}

impl<R> Index<(isize, isize)> for ReducedComplex2<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = ChainComplexSummand<R>;
    fn index(&self, i: (isize, isize)) -> &Self::Output {
        self.get(i.into())
    }
}

impl<R> Index<(isize, isize, isize)> for ReducedComplex3<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = ChainComplexSummand<R>;
    fn index(&self, i: (isize, isize, isize)) -> &Self::Output {
        self.get(i.into())
    }
}

#[cfg(test)]
mod tests { 
    use crate::{ChainComplex, RModStr};

    use super::*;

    #[test]
    fn d3() { 
        let c = ChainComplex::<i64>::d3().reduced(false);

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 0);
        assert_eq!(c[2].rank(), 0);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();
    }

    #[test]
    fn s2() { 
        let c = ChainComplex::<i64>::s2().reduced(false);

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 0);
        assert_eq!(c[2].rank(), 1);

        c.check_d_all();
    }

    #[test]
    fn t2() { 
        let c = ChainComplex::<i64>::t2().reduced(false);

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 2);
        assert_eq!(c[2].rank(), 1);

        c.check_d_all();

        assert!(c.d_matrix(1).is_zero());
        assert!(c.d_matrix(2).is_zero());
    }

    #[test]
    fn rp2() { 
        let c = ChainComplex::<i64>::rp2().reduced(false);
        
        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 1);
        assert_eq!(c[2].rank(), 1);

        c.check_d_all();

        assert!( c.d_matrix(0).is_zero());
        assert!( c.d_matrix(1).is_zero());
        assert!(!c.d_matrix(2).is_zero());

        let a = c.d_matrix(2).to_dense()[[0, 0]];
        assert!(a == 2 || a == -2);
    }
    
    #[test]
    fn bypass() { 
        let c = ChainComplex::<i32>::t2();
        let r = ReducedComplex::bypass(&c);

        r.check_d_all();

        assert_eq!(c[0].rank(), r[0].rank());
        assert_eq!(c[1].rank(), r[1].rank());
        assert_eq!(c[2].rank(), r[2].rank());

        assert_eq!(c.d_matrix(0), r.d_matrix(0));
        assert_eq!(c.d_matrix(1), r.d_matrix(1));
        assert_eq!(c.d_matrix(2), r.d_matrix(2));
    }

    #[test]
    fn bypass_reduced() { 
        let c = ChainComplex::<i32>::t2();
        let r = ReducedComplex::bypass(&c);
        let r = r.reduced(true);

        r.check_d_all();

        assert_eq!(r[0].rank(), 1);
        assert_eq!(r[1].rank(), 2);
        assert_eq!(r[2].rank(), 1);
    }
}
