use clap::ValueEnum;
use derive_more::Display;

#[derive(Clone, Copy, PartialEq, Eq, ValueEnum, Display, Debug, Default)]
#[clap(rename_all="verbatim")]
pub enum CType { 
    #[default] Z, 
    Q, F2, F3
}

impl CType { 
    pub fn is_field(&self) -> bool { 
        use CType::*;
        match self {
            Q | F2 | F3 => true,
            _ => false
        }
    }
}

#[derive(PartialEq, Eq)]
pub(crate) enum PolyVars { 
    H, T, HT, None
}

pub(crate) fn poly_vars(c_value: &String) -> PolyVars { 
    use std::collections::HashSet;
    
    let s: HashSet<_> = c_value.split(',').collect();
    match (s.contains("H"), s.contains("T")) { 
        (true,  true)  => PolyVars::HT,
        (true,  false) => PolyVars::H,
        (false, true)  => PolyVars::T,
        (false, false) => PolyVars::None
    }
}

#[derive(Clone, Copy, PartialEq, Eq, ValueEnum, Display, Debug, Default)]
#[clap(rename_all="lower")]
pub enum Format { 
    #[default] Unicode, 
    TeX
}