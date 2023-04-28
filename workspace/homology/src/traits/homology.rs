use yui_core::{Ring, RingOps};
use crate::{RModStr, RModGrid, ChainComplex};

pub trait Homology: RModGrid
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>
{}

pub trait HomologyComputable: ChainComplex
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>,
    Self::HomologySummand: RModStr<R = Self::R>,
    Self::Homology: Homology<R = Self::R, Idx = Self::Idx, IdxIter = Self::IdxIter, Output = Self::HomologySummand>
{
    type Homology;
    type HomologySummand;

    fn homology(&self) -> Self::Homology;
    fn homology_at(&self, i: Self::Idx) -> Self::HomologySummand;
}