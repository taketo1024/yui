use clap::ValueEnum;
use serde::Deserialize;

#[allow(non_camel_case_types)]
#[derive(Clone, ValueEnum, Debug, derive_more::Display, Deserialize)]
pub enum CType { 
    Z, Q, F2, F3, Gauss, Eisen, None
}