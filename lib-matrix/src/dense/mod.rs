//! Dense matrices and their decompositions.
//!
//! [`Mat<R>`] is a thin wrapper around `nalgebra::DMatrix`. The submodules
//! provide PLUQ decomposition over a ring (with a linear solver over a field),
//! Smith normal form, and the LLL algorithm with Hermite-normal-form variant.

mod mat;
pub use mat::{Mat, MatTrait};

pub mod pluq;
pub mod snf;
pub mod lll;