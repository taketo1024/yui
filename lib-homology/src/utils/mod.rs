//! Shared utilities: the sparse [`Grid`] storage, plain-text/LaTeX
//! formatting helpers, and the [`ToSeqString`] / [`ToTableString`] display
//! traits.

mod matrix;
pub(crate) mod format;
mod grid;
mod to_string;

#[cfg(feature = "tex")]
pub mod tex;

pub use matrix::*;
pub use format::*;
pub use grid::*;
pub use to_string::*;
