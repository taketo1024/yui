//! Khovanov homology of an involutive link: the mapping cone of `1 + τ`, its
//! homology, and the `ssi` invariant read from the canonical classes.

mod khi_gen;
mod complex;
mod homology;
mod tau;

pub use khi_gen::{KhIGen, KhIGenExt};
pub use complex::{KhIComplex, KhIChain};
pub use homology::KhIHomology;
