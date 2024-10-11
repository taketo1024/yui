use log::info;
use yui::{EucRing, EucRingOps};

use crate::utils::HomologyCalc;
use crate::{isize2, isize3, rmod_str_symbol, Grid, GridDeg, GridTrait, RModStr};
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
        info!("compute H[{i}]: {} -> {} -> {} ..", 
            rmod_str_symbol(self.rank(i - self.d_deg()), &[], "0"), 
            rmod_str_symbol(self.rank(i),                &[], "0"), 
            rmod_str_symbol(self.rank(i + self.d_deg()), &[], "0"), 
        );
        
        let d0 = self.d_matrix(i - self.d_deg());
        let d1 = self.d_matrix(i);
        let h = HomologyCalc::calculate(d0, d1, with_trans);

        info!("  H[{i}] = {}.", 
            rmod_str_symbol(h.rank(), h.tors(), "0")
        );

        h
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