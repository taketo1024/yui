use yui::{EucRing, EucRingOps};

use crate::utils::HomologyCalc;
use crate::{GridDeg, Grid, GridTrait, isize2, isize3};
use super::{ChainComplexTrait, SimpleRModStr, ChainComplexBase};

pub type HomologySummand<R> = SimpleRModStr<R>;
pub type HomologyBase<I, R> = Grid<I, HomologySummand<R>>;

pub type Homology<R>  = HomologyBase<isize,  R>;
pub type Homology2<R> = HomologyBase<isize2, R>;
pub type Homology3<R> = HomologyBase<isize3, R>;

pub trait ComputeHomology<I> {
    type Output;
    fn compute_homology(&self, i: I, with_trans: bool) -> Self::Output;
}

impl<I, R, C> ComputeHomology<I> for C 
where 
    I: GridDeg, 
    R: EucRing, for<'x> &'x R: EucRingOps<R>,
    C: ChainComplexTrait<I, R = R>
{
    type Output = HomologySummand<R>;

    fn compute_homology(&self, i: I, with_trans: bool) -> HomologySummand<R> {
        let i0 = i - self.d_deg();
        let d0 = self.d_matrix(i0);
        let d1 = self.d_matrix(i );
        HomologyCalc::calculate(d0, d1, with_trans)
    }
}

impl<I, R> ChainComplexBase<I, R> 
where I: GridDeg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology_at(&self, i: I, with_trans: bool) -> HomologySummand<R> {
        let mut c = self[i].clone();
        let h = self.compute_homology(i, with_trans);
        c.merge(h, true);
        c
    }

    pub fn homology(&self, with_trans: bool) -> HomologyBase<I, R> {
        HomologyBase::generate(
            self.support(), 
            |i| self.homology_at(i, with_trans)
        )
    }
}