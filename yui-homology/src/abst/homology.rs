use yui::{EucRing, EucRingOps, Ring, RingOps};

use crate::generic::GenericSummand;
use crate::utils::HomologyCalc;
use crate::{GenericHomologyBase, Grid, GridDeg};
use super::ChainComplexTrait;

pub trait ComputeHomology<I, R>
where I: GridDeg, R: Ring, for<'x> &'x R: RingOps<R> {
    fn compute_homology_at(&self, i: I, with_trans: bool) -> GenericSummand<I, R>;
    fn compute_homology(&self, with_trans: bool) -> GenericHomologyBase<I, R>;
}

impl<I, R, C> ComputeHomology<I, R> for C 
where 
    I: GridDeg, 
    R: EucRing, for<'x> &'x R: EucRingOps<R>,
    C: ChainComplexTrait<I, R = R>
{
    fn compute_homology_at(&self, i: I, with_trans: bool) -> GenericSummand<I, R> {
        let i0 = i - self.d_deg();
        let d0 = self.d_matrix(i0);
        let d1 = self.d_matrix(i);
        let (rank, tors, trans) = HomologyCalc::calculate(d0, d1, with_trans);
        GenericSummand::generate(i, rank, tors, trans)
    }

    fn compute_homology(&self, with_trans: bool) -> GenericHomologyBase<I, R> {
        Grid::generate(self.support(), |i| self.compute_homology_at(i, with_trans))
    }
}