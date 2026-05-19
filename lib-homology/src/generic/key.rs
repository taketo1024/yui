use derive_more::derive::{Display, Debug};
use yui_core::lc::LcKey;
use yui_core::MathType;

use crate::GridDeg;

#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Display, Debug, Default)]
#[display("e({},{})", _0, _1)]
#[  debug("e({},{})", _0, _1)]
pub struct GenericKey<I>(pub I, pub usize) 
where I: GridDeg;

impl<I> MathType for GenericKey<I>
where I: GridDeg {
    fn math_symbol() -> String {
        "E".into()
    }
}

impl<I> LcKey for GenericKey<I>
where I: GridDeg {}