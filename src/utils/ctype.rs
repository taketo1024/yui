use clap::ValueEnum;
use serde::Deserialize;

#[derive(Clone, ValueEnum, Debug, derive_more::Display, Deserialize)]
#[clap(rename_all="verbatim")]
pub enum CType { 
    Z, Q, F2, F3, Gauss, Eisen
}