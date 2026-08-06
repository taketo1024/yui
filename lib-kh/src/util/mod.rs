//! Internal support types for the Khovanov computations.

pub(crate) mod bigraded;
pub(crate) use bigraded::Bigraded;

pub(crate) mod cached_hash;
pub(crate) use cached_hash::CachedHash;

pub(crate) mod hash_cons;



pub mod fast_poly;
pub use fast_poly::FastPoly;