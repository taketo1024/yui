pub mod dense;
pub mod sparse;
pub mod snf;
pub mod triang;
pub mod pivot;
pub mod schur;

pub use dense::DnsMat;
pub use sparse::CsMatElem;
pub use snf::{snf, snf_in_place};
