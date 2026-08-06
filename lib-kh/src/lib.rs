#![doc = include_str!("../README.md")]

// laid out as `<group>/<type>.rs`, so a module may share its parent's name.
#![allow(clippy::module_inception)]

mod ext;

pub mod kh;
pub mod khi;
pub mod ss;
pub mod tng;
pub mod util;