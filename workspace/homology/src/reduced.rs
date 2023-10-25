use std::collections::HashMap;
use delegate::delegate;

use yui_core::{Deg, Ring, RingOps, EucRing, EucRingOps, isize2, isize3};
use yui_matrix::sparse::{Trans, SpMat, SpVec};

use crate::{ChainComplexBase, Graded, ChainComplexTrait, HomologySummand, HomologyBase};

pub type ReducedComplex<R>  = ReducedComplexBase<isize,  R>;
pub type ReducedComplex2<R> = ReducedComplexBase<isize2, R>;
pub type ReducedComplex3<R> = ReducedComplexBase<isize3, R>;

pub struct ReducedComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    inner: ChainComplexBase<I, R>,
    trans: Option<HashMap<I, Trans<R>>>
}

impl<I, R> ReducedComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new<It, F1, F2>(support: It, d_deg: I, d_matrix: F1, trans: Option<F2>) -> Self
    where 
        It: Iterator<Item = I>, 
        F1: FnMut(I) -> SpMat<R>,
        F2: FnMut(I) -> Trans<R>
    {
        let inner = ChainComplexBase::new(support, d_deg, d_matrix);
        let trans = trans.map(|mut t| { 
            HashMap::from_iter( 
                inner.support().map(|i| (i, t(i))) 
            )
        });
        
        Self { inner, trans }
    }

    pub fn trans(&self, i: I) -> Option<&Trans<R>> {
        self.trans.as_ref().map(|ts| &ts[&i])
    }

    pub fn trans_forward(&self, i: I, v: &SpVec<R>) -> SpVec<R> { 
        assert!(self.trans.is_some());
        self.trans(i).unwrap().forward(v)
    }

    pub fn trans_backward(&self, i: I, v: &SpVec<R>) -> SpVec<R> { 
        assert!(self.trans.is_some());
        self.trans(i).unwrap().backward(v)
    }
}

impl<I, R> Graded<I> for ReducedComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type Itr = std::vec::IntoIter<I>;

    delegate! { 
        to self.inner {
            fn support(&self) -> Self::Itr;
            fn display(&self, i: I) -> String;
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

impl<I, R> ReducedComplexBase<I, R>
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    delegate! { 
        to self.inner { 
            pub fn homology_at(&self, i: I, with_trans: bool) -> HomologySummand<R>;
            pub fn homology(&self, with_trans: bool) -> HomologyBase<I, R>;
        }
    }
}