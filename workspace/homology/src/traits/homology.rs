use yui_core::{Ring, RingOps, EucRing, EucRingOps};
use crate::{RModStr, RModGrid, ChainComplex, GenericRModStr};
use crate::utils::homology_calc::HomologyCalc;

pub trait Homology: RModGrid
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>
{
    fn is_zero(&self) -> bool {
        self.indices().all(|i| self[i].is_zero())
    }

    fn is_free(&self) -> bool {
        self.indices().all(|i| self[i].is_free())
    }

    fn fmt_default(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in self.indices() { 
            write!(f, "H[{}]: {}\n", i, self[i])?
        }
        Ok(())
    }
}

pub trait HomologyComputable<S>: ChainComplex
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>,
    S: RModStr<R = Self::R>
{
    fn homology_at(&self, i: Self::Idx) -> S;
}

impl<R, C> HomologyComputable<GenericRModStr<R>> for C
where 
    R: EucRing, for<'x> &'x R: EucRingOps<R>,
    C: ChainComplex<R = R>,
    C::Output: RModStr<R = R>
{ 
    fn homology_at(&self, k: C::Idx) -> GenericRModStr<Self::R> {
        let d1 = self.d_matrix(k - self.d_degree());
        let d2 = self.d_matrix(k);
        HomologyCalc::calculate(d1, d2)
    }
}

