use num_rational::Ratio;
use num_traits::{Zero, One};

use crate::math::traits::{Ring, RingOps, Symbol, AlgBase, MonOps, Mon, AddMon, AddMonOps, AddGrpOps, AddGrp, RingMethods};

use super::int_ext::{Integer, IntOps};

impl<R> Symbol for Ratio<R>
where R: Integer, for<'x> &'x R: IntOps<R> {
    fn symbol() -> String {
        String::from(format!("Frac({})", R::symbol()))
    }
}

impl<R> AlgBase for Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<R> MonOps<Self> for Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<'a, R> MonOps<Ratio<R>> for &'a Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<R> Mon for Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<R> AddMonOps<Self> for Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<'a, R> AddMonOps<Ratio<R>> for &'a Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<R> AddMon for Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<R> AddGrpOps<Self> for Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<'a, R> AddGrpOps<Ratio<R>> for &'a Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<R> AddGrp for Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<R> RingOps<Self> for Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<'a, R> RingOps<Ratio<R>> for &'a Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<R> RingMethods for Ratio<R>
where R: Integer, for<'x> &'x R: IntOps<R> {
    fn inv(&self) -> Option<Self> {
        if !self.is_zero() { 
            Some( num_traits::Inv::inv(self) )
        } else {
            None
        }
    }

    fn is_unit(&self) -> bool {
        !self.is_zero()
    }

    fn normalizing_unit(&self) -> Self {
        if let Some(u) = self.inv() { 
            u
        } else { 
            Self::one()
        }
    }
}

impl<R> Ring for Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

