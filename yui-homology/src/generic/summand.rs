use yui::{Ring, RingOps};
use yui_matrix::sparse::Trans;

use crate::{GridDeg, Summand};

use super::EnumGen;

pub type GenericSummand<I, R> = Summand<EnumGen<I>, R>;

impl<I, R> GenericSummand<I, R>
where I: GridDeg, R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn generate(i: I, rank: usize, tors: Vec<R>, trans: Option<Trans<R>>) -> Self { 
        let (n, trans) = if let Some(t) = trans { 
            (t.src_dim(), t)
        } else { 
            let n = rank + tors.len();
            (n, Trans::id(n))
        };
        
        let gens = (0..n).map(|j| EnumGen(i, j)).collect();
        Self::new(gens, rank, tors, trans)
    }

    pub fn generate_free(i: I, rank: usize) -> Self { 
        Self::generate(i, rank, vec![], None)
    }
}
