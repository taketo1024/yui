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

    fn is_free_at(&self, i: Self::Idx) -> bool { 
        self.get(i).map(|v| v.is_free()).unwrap_or(true)
    }

    fn is_zero(&self) -> bool { 
        self.iter().all(|(_, v)| v.is_zero())
    }

    fn is_zero_at(&self, i: Self::Idx) -> bool { 
        self.get(i).map(|v| v.is_zero()).unwrap_or(true)
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