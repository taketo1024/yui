use yui::{EucRing, EucRingOps};
use yui::lc::Gen;

use crate::{ComputeHomology, Grid, GridDeg, GridTrait, SummandTrait};
use super::{Summand, ChainComplexBase};

impl<I, X, R> ChainComplexBase<I, X, R>
where 
    I: GridDeg,
    X: Gen,
    R: EucRing, for<'x> &'x R: EucRingOps<R>
{
    pub fn homology_at(&self, i: I) -> Summand<X, R> {
        let c = &self[i];
        let h = self.compute_homology_at(i, true);

        Summand::new(
            c.raw_gens().clone(),
            h.rank(),
            h.tors().iter().cloned().collect(),
            c.trans().merged(h.trans())
        )
    }

    pub fn homology(&self) -> Grid<I, Summand<X, R>> {
        Grid::generate(
            self.support(), 
            |i| self.homology_at(i)
        )
    }
}