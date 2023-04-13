use num_traits::Zero;
use yui_matrix::sparse::*;
use yui_core::{RingOps, Ring};
use crate::{RModStr, ChainComplex};

pub trait ChainComplexValidation: ChainComplex
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>
{
    fn is_cycle(&self, k: Self::Idx, z: &SpVec<Self::R>) -> bool { 
        let d = self.d_matrix(k);
        (&d * z).iter().all(|(_, a)| a.is_zero())
    }
    
    fn check_d_at(&self, k: Self::Idx) { 
        let d1 = self.d_matrix(k);
        let d2 = self.d_matrix(k + self.d_degree());
        let res = &d2 * &d1;
        assert!( res.is_zero(), "d² is non-zero at {k}." );
    }

    fn check_d_all(&self) {
        for k in self.indices() { 
            let k2 = k + self.d_degree();
            if !self.contains_idx(k2) { continue }
            self.check_d_at(k);
        }
    }
}

impl<R, C> ChainComplexValidation for C
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    C: ChainComplex<R = R>,
    C::Output: RModStr<R = R>
{}