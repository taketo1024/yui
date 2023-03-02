use yui_matrix::sparse::*;
use yui_core::{RingOps, Ring};
use crate::{RModStr, RModGrid};

pub trait ChainComplex: RModGrid
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>
{
    fn d_degree(&self) -> Self::Idx;
    fn d_matrix(&self, k: Self::Idx) -> SpMat<Self::R>;
}