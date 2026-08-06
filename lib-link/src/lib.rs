#![doc = include_str!("../README.md")]

// laid out as `<group>/<type>.rs`, so a module may share its parent's name.
#![allow(clippy::module_inception)]

mod link;
mod inv_link;
mod braid;
pub mod misc;

#[cfg(any(test, feature = "test-utils"))]
mod test_data;

pub use link::*;
pub use inv_link::*;
pub use braid::*;