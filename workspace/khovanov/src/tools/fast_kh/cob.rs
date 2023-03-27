use std::fmt::Display;
use derive_more::Display;
use itertools::Itertools;
use yui_core::Elem;
use yui_lin_comb::{FreeGen, OrdForDisplay};

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Display)]
pub enum Dot { 
    X, Y
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct CobComp { 
    src: Vec<usize>, // indices of comps in the source tangle. 
    tgt: Vec<usize>,
    dots: (usize, usize) // nums of X and Y dots resp. 
}

impl Display for CobComp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let dots = vec!["X"; self.dots.0].join("") + &vec!["Y"; self.dots.1].join("");
        write!(f, "{{{:?} -> {:?}}}{}", self.src, self.tgt, dots)
    }
}

#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct Cob { 
    comps: Vec<CobComp>
}

impl Cob {
}

impl Display for Cob {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let cobs = self.comps.iter().map(|c| c.to_string()).join(", ");
        write!(f, "cob[{}]", cobs)
    }
}

impl OrdForDisplay for Cob {
    fn cmp_for_display(&self, _other: &Self) -> std::cmp::Ordering {
        std::cmp::Ordering::Equal
    }
}

impl Elem for Cob {
    fn set_symbol() -> String {
        "Cob".to_string()
    }
}

impl FreeGen for Cob {}