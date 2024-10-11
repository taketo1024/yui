use derive_more::derive::Display;
use yui::lc::Gen;
use yui::Elem;

#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Display, Debug, Default)]
#[display("e({},{})", _0, _1)]
pub struct EnumGen(pub isize, pub usize);

impl Elem for EnumGen {
    fn math_symbol() -> String {
        "E".into()
    }
}

impl Gen for EnumGen {}