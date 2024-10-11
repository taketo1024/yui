use yui::{Ring, RingOps};

use crate::{GridDeg, XModStr};

use super::EnumGen;

pub type GenericSummand<I, R> = XModStr<EnumGen<I>, R>;

impl<I, R> GenericSummand<I, R>
where I: GridDeg, R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn generate_free(i: I, rank: usize) -> Self { 
        Self::free((0..rank).map(|j| EnumGen(i, j)))
    }
}
