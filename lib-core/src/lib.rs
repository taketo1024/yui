#![doc = include_str!("../README.md")]

// laid out as `<group>/<type>.rs`, so a module may share its parent's name.
#![allow(clippy::module_inception)]

pub mod abst;
pub mod ext;

mod conc;
pub use conc::*;

pub mod algo;
pub mod util;