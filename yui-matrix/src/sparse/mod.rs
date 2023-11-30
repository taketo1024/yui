pub use crate::dense::MatType;

// mod sp_mat;
// mod sp_vec;
// pub use sp_mat::{SpMat, SpMatView};
// pub use sp_vec::{SpVec, SpVecView};

mod _sp_mat;
mod _sp_vec;
pub use _sp_mat::{SpMat, SpMatView};
pub use _sp_vec::{SpVec, SpVecView};


mod trans;
pub use trans::*;

pub mod pivot;
pub mod schur;
pub mod triang;

