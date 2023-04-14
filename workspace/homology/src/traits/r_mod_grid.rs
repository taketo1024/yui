use yui_core::{Ring, RingOps};
use crate::{Grid, RModStr};

pub trait RModGrid: Grid
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>,
{
    type R;

    fn is_free(&self) -> bool { 
        self.iter().all(|(_, v)| v.is_free())
    }

    fn is_zero(&self) -> bool { 
        self.iter().all(|(_, v)| v.is_zero())
    }
}

impl<R, T> RModGrid for T
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    T: Grid, 
    T::Output: RModStr<R = R>
{
    type R = R;
}