use clap::ValueEnum;
use derive_more::Display;

pub trait AppArgs { 
    fn c_type(&self) -> CType; 
    fn c_value(&self) -> &String; 
    fn log(&self) -> u8; 

    fn poly_vars(&self) -> PolyVars { 
        parse_poly_vars(self.c_value())
    }

    fn is_poly(&self) -> bool { 
        self.poly_vars() != PolyVars::None
    }

    fn is_field(&self) -> bool { 
        self.c_type().is_field() && !self.is_poly()
    }

    fn is_euc_ring(&self) -> bool { 
        self.c_type() == CType::Z && !self.is_poly() || 
        self.c_type().is_field() && self.poly_vars().nvars() == 1
    }

    fn log_level(&self) -> log::LevelFilter { 
        use log::LevelFilter::*;
        match self.log() {
            1 => Info,
            2 => Debug,
            3 => Trace,
            _ => Off,
        }
    }
}

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

impl PolyVars {
    pub fn nvars(&self) -> usize { 
        match self {
            PolyVars::H | PolyVars::T  => 1,
            PolyVars::HT => 2,
            _ => 0,
        }
    }
}

pub(crate) fn parse_poly_vars(c_value: &String) -> PolyVars { 
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