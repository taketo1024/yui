use std::ops::Index;

use delegate::delegate;

use yui_core::{Deg, Ring, RingOps, EucRing, EucRingOps, isize2, isize3};
use yui_matrix::sparse::{Trans, SpMat, SpVec};

use crate::{ChainComplexBase, GridTrait, ChainComplexTrait, HomologyBase, GridBase, GridIter, ChainComplexSummand, ChainComplexSummandTrait};

pub type ReducedComplex<R>  = ReducedComplexBase<isize,  R>;
pub type ReducedComplex2<R> = ReducedComplexBase<isize2, R>;
pub type ReducedComplex3<R> = ReducedComplexBase<isize3, R>;

pub struct ReducedComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    inner: ChainComplexBase<I, R>,
    trans: Option<GridBase<I, Trans<R>>>
}

impl<I, R> ReducedComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new<It, F1, F2>(support: It, d_deg: I, d_matrix: F1, trans_opt: Option<F2>) -> Self
    where 
        It: Iterator<Item = I>, 
        F1: FnMut(I) -> SpMat<R>,
        F2: FnMut(I) -> Trans<R>
    {
        let inner = ChainComplexBase::new(support, d_deg, d_matrix);
        let trans = trans_opt.map(|t| { 
            GridBase::new(inner.support(), t)
        });
        
        Self { inner, trans }
    }

    pub fn reduced_by<C, F>(complex: &C, trans_map: F) -> Self
    where 
        C: ChainComplexTrait<I, R = R>,
        C::E: ChainComplexSummandTrait<R = R>,
        F: FnMut(I) -> Trans<R>
    {
        //              d
        //     C[i] ---------> C[i+1]
        //      ^                |
        // t_in |                | t_out
        //      |                V
        //     C'[i] - - - - > C'[i+1]

        let trans = GridBase::new(complex.support(), trans_map);        
        let inner = ChainComplexBase::new(complex.support(), complex.d_deg(), |i| { 
            let d = complex.get(i).d_matrix();
            let i1 = i + complex.d_deg();
            let t_in  = trans.get(i).backward_mat();
            let t_out = trans.get(i1).forward_mat();

            t_out * d * t_in
        });

        Self { inner, trans: Some(trans) }
    }

    pub fn bypass<C>(complex: &C, with_trans: bool) -> Self
    where 
        C: ChainComplexTrait<I, R = R>,
        C::E: ChainComplexSummandTrait<R = R>
    {
        let trans = if with_trans { 
            Some(|i| Trans::id(complex.get(i).rank()))
        } else { 
            None
        };
        Self::new(
            complex.support(), 
            complex.d_deg(), 
            |i| complex.get(i).d_matrix().clone(), 
            trans
        )
    }

    pub fn trans_at(&self, i: I) -> Option<&Trans<R>> {
        self.trans.as_ref().map(|ts| ts.get(i))
    }

    pub fn trans_forward(&self, i: I, v: &SpVec<R>) -> SpVec<R> { 
        assert!(self.trans.is_some());
        self.trans_at(i).unwrap().forward(v)
    }

    pub fn trans_backward(&self, i: I, v: &SpVec<R>) -> SpVec<R> { 
        assert!(self.trans.is_some());
        self.trans_at(i).unwrap().backward(v)
    }
}

impl<I, R> Index<I> for ReducedComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = ChainComplexSummand<R>;
    fn index(&self, i: I) -> &Self::Output {
        self.get(i)
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
        }
    }
}

impl<I, R> ReducedComplexBase<I, R>
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    delegate! { 
        to self.inner { 
            pub fn homology(&self, with_trans: bool) -> HomologyBase<I, R>;
        }
    }
}

#[cfg(test)]
mod tests { 
    use crate::ChainComplex;
    use crate::utils::ChainReducer;

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

        assert!(c[1].d_matrix().is_zero());
        assert!(c[2].d_matrix().is_zero());
    }

    #[test]
    fn rp2() { 
        let c = ChainComplex::<i64>::rp2().reduced(false);
        
        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 1);
        assert_eq!(c[2].rank(), 1);

        c.check_d_all();

        assert!( c[0].d_matrix().is_zero());
        assert!( c[1].d_matrix().is_zero());
        assert!(!c[2].d_matrix().is_zero());

        let a = c[2].d_matrix().to_dense()[[0, 0]];
        assert!(a == 2 || a == -2);
    }
    
    #[test]
    fn reduced_by() { 
        let c = ChainComplex::<i32>::rp2();

        let mut f = ChainReducer::new(&c, true);
        f.process_all();

        let c = c.reduced_by(|i| f.trans(i).unwrap().clone());

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 1);
        assert_eq!(c[2].rank(), 1);

        c.check_d_all();

        assert!( c[0].d_matrix().is_zero());
        assert!( c[1].d_matrix().is_zero());
        assert!(!c[2].d_matrix().is_zero());

        let a = c[2].d_matrix().to_dense()[[0, 0]];
        assert!(a == 2 || a == -2);
    }

    #[test]
    fn bypass() { 
        let c = ChainComplex::<i32>::t2();
        let r = ReducedComplex::bypass(&c, false);

        r.check_d_all();

        assert_eq!(c[0].rank(), r[0].rank());
        assert_eq!(c[1].rank(), r[1].rank());
        assert_eq!(c[2].rank(), r[2].rank());

        assert_eq!(c[0].d_matrix(), r[0].d_matrix());
        assert_eq!(c[1].d_matrix(), r[1].d_matrix());
        assert_eq!(c[2].d_matrix(), r[2].d_matrix());
    }
}
