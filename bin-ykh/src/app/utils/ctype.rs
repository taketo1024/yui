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

#[derive(Clone, Copy, PartialEq, Eq, ValueEnum, Display, Debug, Default)]
#[clap(rename_all="lower")]
pub enum Format { 
    #[default] Unicode, 
    TeX
}