use log::info;
use yui::{EucRing, EucRingOps, Ring, RingOps};

use crate::generic::GenericSummand;
use crate::utils::HomologyCalc;
use crate::{GridDeg, Grid, isize2, isize3, rmod_str_symbol};
use super::{ChainComplexTrait, SimpleRModStr};

pub type HomologySummand<R> = SimpleRModStr<R>;
pub type HomologyBase<I, R> = Grid<I, HomologySummand<R>>;

pub type Homology<R>  = HomologyBase<isize,  R>;
pub type Homology2<R> = HomologyBase<isize2, R>;
pub type Homology3<R> = HomologyBase<isize3, R>;

pub trait ComputeHomology<I, R>
where I: GridDeg, R: Ring, for<'x> &'x R: RingOps<R> {
    fn compute_homology_at(&self, i: I, with_trans: bool) -> GenericSummand<I, R>;
    fn compute_homology(&self, with_trans: bool) -> Grid<I, GenericSummand<I, R>>;
}

impl<I, R, C> ComputeHomology<I, R> for C 
where 
    I: GridDeg, 
    R: EucRing, for<'x> &'x R: EucRingOps<R>,
    C: ChainComplexTrait<I, R = R>
{
    fn compute_homology_at(&self, i: I, with_trans: bool) -> GenericSummand<I, R> {
        info!("compute H[{i}]: {} -> {} -> {} ..", 
            rmod_str_symbol(self.rank(i - self.d_deg()), &[], "0"), 
            rmod_str_symbol(self.rank(i),                &[], "0"), 
            rmod_str_symbol(self.rank(i + self.d_deg()), &[], "0"), 
        );
        
        let i0 = i - self.d_deg();
        let d0 = self.d_matrix(i0);
        let d1 = self.d_matrix(i);
        let (rank, tors, trans) = HomologyCalc::calculate(d0, d1, with_trans);

        info!("  H[{i}] = {}.", 
            rmod_str_symbol(rank, &tors, "0")
        );

        GenericSummand::generate(i, rank, tors, trans)
    }

    fn compute_homology(&self, with_trans: bool) -> Grid<I, GenericSummand<I, R>> {
        Grid::generate(self.support(), |i| self.compute_homology_at(i, with_trans))
    }
}