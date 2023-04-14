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

    fn fmt_default(&self, f: &mut std::fmt::Formatter<'_>, s: &str) -> std::fmt::Result {
        for (i, v) in self.iter() { 
            write!(f, "{s}[{}]: {}\n", i, v)?
        }
        Ok(())
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