use yui::{EucRing, EucRingOps};
use yui::lc::Gen;

use crate::{GridDeg, isize2, isize3, Grid, GridTrait, ComputeHomology};
use super::{XModStr, XChainComplexBase};

pub type XHomologySummand<X, R> = XModStr<X, R>;
pub type XHomologyBase<I, X, R> = Grid<I, XHomologySummand<X, R>>;

pub type XHomology <X, R> = XHomologyBase<isize,  X, R>;
pub type XHomology2<X, R> = XHomologyBase<isize2, X, R>;
pub type XHomology3<X, R> = XHomologyBase<isize3, X, R>;

impl<I, X, R> XChainComplexBase<I, X, R>
where 
    I: GridDeg,
    X: Gen,
    R: EucRing, for<'x> &'x R: EucRingOps<R>
{
    pub fn homology_at(&self, i: I, with_trans: bool) -> XHomologySummand<X, R> {
        let mut ci = self[i].clone();
        let hi = self.compute_homology(i, with_trans);
        ci.merge(hi, true);
        ci
    }

    pub fn homology(&self, with_trans: bool) -> XHomologyBase<I, X, R> {
        XHomologyBase::generate(
            self.support(), 
            |i| self.homology_at(i, with_trans)
        )
    }
}