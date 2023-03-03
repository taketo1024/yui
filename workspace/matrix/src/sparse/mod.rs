mod sp_mat;
mod sp_vec;

pub use sp_mat::{SpMat, SpMatView};
pub use sp_vec::{SpVec, SpVecView};
pub use crate::dense::MatType;

pub mod pivot;
pub mod schur;
pub mod triang;