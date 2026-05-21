use yui_core::{EucRing, EucRingOps};
use yui_core::lc::LcKey;

use crate::generic::GenericSummand;
use crate::algo::HomologyCalc;
use crate::{isize2, isize3, GenericHomologyBase, GrMod, AddInd};
use super::{Summand, ChainComplexBase};

pub type HomologyBase<I, X, R> = GrMod<I, X, R>;
pub type Homology <X, R> = HomologyBase<isize,  X, R>;
pub type Homology2<X, R> = HomologyBase<isize2, X, R>;
pub type Homology3<X, R> = HomologyBase<isize3, X, R>;

impl<I, X, R> ChainComplexBase<I, X, R>
where
    I: AddInd,
    X: LcKey,
    R: EucRing, for<'x> &'x R: EucRingOps<R>
{
    pub fn homology_at(&self, i: I) -> Summand<X, R> {
        let c = &self[i];
        let d0 = self.d_matrix(i - self.d_deg());
        let d1 = self.d_matrix(i);
        let (rank, tors, trans) = HomologyCalc::calculate(d0, d1, true);
        let trans = trans.unwrap();

        Summand::new(
            c.raw_generators().clone(),
            rank,
            tors,
            c.trans().merged(&trans)
        )
    }

    pub fn homology(&self) -> HomologyBase<I, X, R> {
        GrMod::generate_filtered(
            self.support().copied(),
            |i| {
                let hi = self.homology_at(i);
                (!hi.is_zero()).then_some(hi)
            }
        )
    }

    /// Compute the homology at index `i` without tracking the basis change.
    /// Faster than `homology_at` when only the rank/torsion are needed.
    pub fn generic_homology_at(&self, i: I) -> GenericSummand<I, R> {
        let d0 = self.d_matrix(i - self.d_deg());
        let d1 = self.d_matrix(i);
        let (rank, tors, _) = HomologyCalc::calculate(d0, d1, false);
        GenericSummand::generate(i, rank, tors, None)
    }

    pub fn generic_homology(&self) -> GenericHomologyBase<I, R> {
        GrMod::generate_filtered(
            self.support().copied(),
            |i| {
                let hi = self.generic_homology_at(i);
                (!hi.is_zero()).then_some(hi)
            }
        )
    }
}