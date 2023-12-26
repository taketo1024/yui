pub use crate::MatTrait;

mod sp_mat;
mod sp_vec;
pub use sp_mat::SpMat;
pub use sp_vec::{SpVec, SpVecView};


mod trans;
pub use trans::*;

pub mod pivot;
pub mod schur;
pub mod triang;

