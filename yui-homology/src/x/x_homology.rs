use yui::{EucRing, EucRingOps};
use yui::lc::Gen;

use crate::{isize2, isize3, ComputeHomology, Grid, GridDeg, GridTrait, RModStr};
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
    pub fn homology_at(&self, i: I) -> XHomologySummand<X, R> {
        let c = &self[i];
        let h = self.compute_homology_at(i, true);

        XHomologySummand::new(
            c.gens().clone(),
            h.rank(),
            h.tors().iter().cloned().collect(),
            // Some(c.trans().unwrap().merged(&h.trans().unwrap()))
            todo!()
        )
    }

    pub fn homology(&self) -> XHomologyBase<I, X, R> {
        XHomologyBase::generate(
            self.support(), 
            |i| self.homology_at(i)
        )
    }
}