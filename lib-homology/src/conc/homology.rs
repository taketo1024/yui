use yui_core::{EucRing, EucRingOps};
use yui_core::lc::Gen;

use crate::{isize2, isize3, ComputeHomology, Grid, GridDeg, GridTrait, SummandTrait};
use super::{Summand, ChainComplexBase};

pub type HomologyBase<I, X, R> = Grid<I, Summand<X, R>>;
pub type Homology <X, R> = HomologyBase<isize,  X, R>;
pub type Homology2<X, R> = HomologyBase<isize2, X, R>;
pub type Homology3<X, R> = HomologyBase<isize3, X, R>;

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

    pub fn homology(&self) -> HomologyBase<I, X, R> {
        Grid::generate(
            self.support(), 
            |i| self.homology_at(i)
        )
    }
}