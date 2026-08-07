//! Concrete instantiations of the abstract hierarchy: numbers, polynomials,
//! linear combinations, and small packed types.

pub mod num;
pub mod lc;
pub mod poly;

mod misc;
pub use misc::*;