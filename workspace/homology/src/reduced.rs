use std::ops::Index;

use delegate::delegate;

use yui_core::{Deg, Ring, RingOps, EucRing, EucRingOps, isize2, isize3};
use yui_matrix::sparse::{Trans, SpMat, SpVec};

use crate::utils::ChainReducer;
use crate::{ChainComplexBase, GridTrait, ChainComplexTrait, HomologyBase, GridBase, GridIter, ChainComplexSummand, ChainComplexDisplay, RModStr};

pub type ReducedComplex<R>  = ReducedComplexBase<isize,  R>;
pub type ReducedComplex2<R> = ReducedComplexBase<isize2, R>;
pub type ReducedComplex3<R> = ReducedComplexBase<isize3, R>;

pub struct ReducedComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    inner: ChainComplexBase<I, R>,
    trans: GridBase<I, Trans<R>>
}

impl<I, R> ReducedComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new<It, F1, F2>(support: It, d_deg: I, d_matrix_map: F1, trans_map: F2) -> Self
    where 
        It: Iterator<Item = I>, 
        F1: FnMut(I) -> SpMat<R>,
        F2: FnMut(I) -> Trans<R>
    {
        let inner = ChainComplexBase::new(support, d_deg, d_matrix_map);
        let trans = GridBase::new(inner.support(), trans_map);
        Self { inner, trans }
    }

    pub fn reduced_by<C, F>(complex: &C, trans_map: F) -> Self
    where 
        C: ChainComplexTrait<I, R = R>,
        C::E: RModStr<R = R>,
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
            let d = complex.d_matrix(i);
            let i1 = i + complex.d_deg();
            let t_in  = trans.get(i).backward_mat();
            let t_out = trans.get(i1).forward_mat();

            t_out * d * t_in
        });

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
            |i| complex.d_matrix(i).clone(), 
            |i| Trans::id(complex.get(i).rank())
        )
    }

    pub fn trans(&self, i: I) -> &Trans<R> {
        &self.trans[i]
    }

    pub fn combine(&self, other: &ReducedComplexBase<I, R>) -> Self { 
        assert!(self.d_deg() == other.d_deg());
        debug_assert!(self.support().all(|i| 
            self.trans(i).tgt_dim() == other.trans(i).src_dim())
        );

        Self::new(
            self.support(), 
            self.d_deg(), 
            |i| other.d_matrix(i).clone(), 
            |i| self.trans(i).compose(other.trans(i))
        )
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

    fn reduced(&self, with_trans: bool) -> ReducedComplexBase<I, Self::R> {
        let red = ChainReducer::reduce(self, with_trans);
        if with_trans { 
            self.combine(&red)
        } else { 
            red
        }
    }

    fn reduced_by<F>(&self, trans: F) -> ReducedComplexBase<I, Self::R>
    where F: FnMut(I) -> Trans<Self::R> {
        let red = Self::reduced_by(self, trans);
        self.combine(&red)
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
    fn reduced_by() { 
        let c = ChainComplex::<i32>::rp2();

        let mut f = ChainReducer::new(&c, true);
        f.process_all();

        let c = c.reduced_by(|i| f.trans(i).unwrap().clone());

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
