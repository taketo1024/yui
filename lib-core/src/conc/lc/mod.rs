//! [`Lc<X, R>`]: a formal linear combination of keys `X` with coefficients in a
//! ring `R` — the element type of every chain module in this workspace.

mod lc_key;
mod lc_data;
mod lc;

pub use lc_key::*;
pub use lc::*;
pub use lc_data::{LcDataIter, LcDataIntoIter};