use std::collections::HashMap;
use derive_more::derive::{Display, Debug};
use num_traits::Zero;
use yui_core::lc::{Lc, LcKey};
use yui_core::{MathType, Ring, RingOps};
use yui_matrix::MatTrait;
use yui_matrix::sparse::{SpMat, Trans};

use crate::{isize2, isize3, AddInd, ChainComplexBase, GrMod, Summand};

// --- GenericKey -------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Display, Debug, Default)]
#[display("e({},{})", _0, _1)]
#[  debug("e({},{})", _0, _1)]
pub struct GenericKey<I>(pub I, pub usize)
where I: AddInd;

impl<I> MathType for GenericKey<I>
where I: AddInd {
    fn math_symbol() -> String {
        "E".into()
    }
}

impl<I> LcKey for GenericKey<I>
where I: AddInd {}

// --- GenericSummand ---------------------------------------------------------

pub type GenericSummand<I, R> = Summand<GenericKey<I>, R>;

impl<I, R> GenericSummand<I, R>
where I: AddInd, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn generate(i: I, rank: usize, tors: Vec<R>, trans: Option<Trans<R>>) -> Self {
        let (n, trans) = if let Some(t) = trans {
            (t.src_dim(), t)
        } else {
            let n = rank + tors.len();
            (n, Trans::id(n))
        };

        let gens = (0..n).map(|j| GenericKey(i, j)).collect();
        Self::new(gens, rank, tors, trans)
    }

    pub fn generate_free(i: I, rank: usize) -> Self {
        Self::generate(i, rank, vec![], None)
    }
}

// --- GenericChainComplex / GenericHomology aliases --------------------------

pub type GenericChainComplexBase<I, R> = ChainComplexBase<I, GenericKey<I>, R>;
pub type GenericChainComplex<R>  = GenericChainComplexBase<isize,  R>;
pub type GenericChainComplex2<R> = GenericChainComplexBase<isize2, R>;
pub type GenericChainComplex3<R> = GenericChainComplexBase<isize3, R>;

pub type GenericHomologyBase<I, R> = GrMod<I, GenericKey<I>, R>;
pub type GenericHomology<R>  = GenericHomologyBase<isize,  R>;
pub type GenericHomology2<R> = GenericHomologyBase<isize2, R>;
pub type GenericHomology3<R> = GenericHomologyBase<isize3, R>;

impl<I, R> GenericChainComplexBase<I, R>
where I: AddInd, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn generate<It, F>(support: It, d_deg: I, mut d_matrix_map: F) -> Self
    where
        It: IntoIterator<Item = I>,
        F: FnMut(I) -> SpMat<R>,
    {
        let d_matrices: HashMap<I, SpMat<R>> = support.into_iter()
            .map(|i| (i, d_matrix_map(i)))
            .collect();

        let summands = GrMod::generate(
            d_matrices.keys().copied(),
            |i| {
                let r = d_matrices[&i].n_cols();
                GenericSummand::generate_free(i, r)
            }
        );

        Self::new(
            summands.clone(), d_deg,
            move |i, z| {
                let Some(d) = d_matrices.get(&i) else {
                    return Lc::zero();
                };
                let v = summands[i].vectorize(z);
                let dv = d * v;
                summands[i + d_deg].devectorize(&dv)
            }
        )
    }
}
