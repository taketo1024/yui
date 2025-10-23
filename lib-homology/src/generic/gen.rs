use derive_more::derive::{Display, Debug};
use yui::lc::Gen;
use yui::Elem;

use crate::GridDeg;

#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Display, Debug, Default)]
#[display("e({},{})", _0, _1)]
#[  debug("e({},{})", _0, _1)]
pub struct EnumGen<I>(pub I, pub usize) 
where I: GridDeg;

impl<I> Elem for EnumGen<I>
where I: GridDeg {
    fn math_symbol() -> String {
        "E".into()
    }
}

impl<I> Gen for EnumGen<I>
where I: GridDeg {}