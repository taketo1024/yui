use num_rational::Ratio;
use num_traits::{Zero, One};
use crate::math::traits::{Ring, RingOps, AlgBase, MonOps, Mon, AddMon, AddMonOps, AddGrpOps, AddGrp};
use super::int_ext::{Integer, IntOps};

macro_rules! impl_alg_ops {
    ($trait:ident) => {
        impl<R> $trait<Self> for Ratio<R> 
        where R: Integer, for<'x> &'x R: IntOps<R> {}

        impl<'a, R> $trait<Ratio<R>> for &'a Ratio<R> 
        where R: Integer, for<'x> &'x R: IntOps<R> {}
    };
}

impl_alg_ops!(AddMonOps);
impl_alg_ops!(AddGrpOps);
impl_alg_ops!(MonOps);
impl_alg_ops!(RingOps);

impl<R> AlgBase for Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {
    fn symbol() -> String {
        format!("Frac({})", R::symbol())
    }
}

impl<R> Mon for Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<R> AddMon for Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<R> AddGrp for Ratio<R> 
where R: Integer, for<'x> &'x R: IntOps<R> {}

impl<R> Ring for Ratio<R> 
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