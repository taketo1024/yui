#![doc = include_str!("../README.md")]

// laid out as `<group>/<type>.rs`, so a module may share its parent's name.
#![allow(clippy::module_inception)]

mod app;
pub use app::{App, CliArgs};
