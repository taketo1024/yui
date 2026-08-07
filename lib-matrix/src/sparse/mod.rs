//! Sparse matrices (CSC) and their decompositions.
//!
//! [`SpMat<R>`] and [`SpVec<R>`] wrap `nalgebra_sparse::CscMatrix`. The
//! submodules provide pivot-based reductions used by the homology pipeline:
//! a heuristic pivot finder ([`pivot`]), PLUQ decomposition with an
//! incremental linear solver ([`pluq`]), Schur-complement reduction
//! ([`schur`]), and triangular solves ([`triang`]).
//!
//! [`Trans<R>`] tracks the forward/backward basis change as a sparse
//! matrix is reduced; it's used to lift solutions back through a chain of
//! reductions.

pub use crate::MatTrait;

mod sp_mat;
mod sp_vec;
pub use sp_mat::SpMat;
pub use sp_vec::SpVec;

mod trans;
pub use trans::*;

pub mod pivot;
pub mod pluq;
pub mod schur;
pub mod snf;
pub mod triang;