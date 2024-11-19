use itertools::Itertools;
use yui::{Ring, RingOps};
use yui_matrix::sparse::SpMat;

use crate::generic::GenericChainComplexBase;
use crate::{GridTrait, GridDeg};
use crate::rmod_str_symbol;

pub trait ChainComplexTrait<I>: Sized + GridTrait<I>
where 
    I: GridDeg, 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> 
{ 
    type R;
    type Element;

    // required methods
    fn rank(&self, i: I) -> usize;
    fn d_deg(&self) -> I;
    fn d(&self, i: I, z: &Self::Element) -> Self::Element;
    fn d_matrix(&self, i: I) -> SpMat<Self::R>;

    // convenient methods
    fn check_d_at(&self, i0: I) { 
        let i1 = i0 + self.d_deg();
        if !(self.is_supported(i0) && self.is_supported(i1)) {
            return 
        }

        let d0 = self.d_matrix(i0);
        let d1 = self.d_matrix(i1);
        let res = d1 * d0;

        assert!( res.is_zero(), "d² is non-zero at {i0}." );
    }

    fn check_d_all(&self) {
        for i in self.support() { 
            self.check_d_at(i);
        }
    }

    fn display_d_at(&self, i: I) -> String {
        let c = |i| rmod_str_symbol(self.rank(i), &[], "0");
        let c0 = c(i);
        let c1 = c(i + self.d_deg());
        let d = self.d_matrix(i).into_dense();
        format!("d[{i}]: {c0} -> {c1}\n{d}")
    }

    fn display_d(&self) -> String { 
        self.support().filter_map(|i| 
            if self.rank(i) > 0 && self.rank(i + self.d_deg()) > 0 && !self.d_matrix(i).is_zero() {
                Some(self.display_d_at(i))
            } else { 
                None
            }
        ).join("")
    }

    fn print_d(&self) {
        println!("{}", self.display_d());
    }

    fn as_generic(&self) -> GenericChainComplexBase<I, Self::R> {
        GenericChainComplexBase::generate(self.support(), self.d_deg(), |i| self.d_matrix(i))
    }
}