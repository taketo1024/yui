pub mod dense;
pub mod sparse;
pub mod snf;

pub use dense::DnsMat;
pub use sparse::CsMatElem;
pub use snf::{snf, snf_in_place};
