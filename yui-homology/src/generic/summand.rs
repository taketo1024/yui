use yui::{Ring, RingOps};
use yui_matrix::sparse::Trans;

use crate::{GridDeg, XModStr};

use super::EnumGen;

pub type GenericSummand<I, R> = XModStr<EnumGen<I>, R>;

impl<I, R> GenericSummand<I, R>
where I: GridDeg, R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn generate(i: I, rank: usize, tors: Vec<R>, trans: Option<Trans<R>>) -> Self { 
        let n = if let Some(t) = &trans { 
            t.src_dim()
        } else { 
            rank + tors.len()
        };
        let gens = (0..n).map(|j| EnumGen(i, j)).collect();
        Self::new(gens, rank, tors, trans)
    }

    pub fn generate_free(i: I, rank: usize) -> Self { 
        Self::generate(i, rank, vec![], Some(Trans::id(rank)))
    }
}
