#![doc = include_str!("../README.md")]

pub mod dense;
pub mod sparse;

pub use dense::MatTrait;

mod perm;
pub use perm::*;