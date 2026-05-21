mod conc;

pub use conc::*;

#[cfg(test)]
mod test_data;

pub mod algo;
pub mod utils;
pub use utils::{ToSeqString, ToTableString, Grid, Grid1, Grid2, Grid3, rmod_str};

#[cfg(feature = "tex")]
pub use utils::tex;
