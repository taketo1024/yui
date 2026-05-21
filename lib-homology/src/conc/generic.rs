use std::collections::HashMap;
use derive_more::derive::{Display, Debug};
use num_traits::Zero;
use yui_core::lc::{Lc, LcKey};
use yui_core::{MathType, Ring, RingOps};
use yui_matrix::MatTrait;
use yui_matrix::sparse::{SpMat, Trans};

use crate::{isize2, isize3, AddInd, ChainComplex, GrMod, Summand};

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

// --- GenericChainComplex / GenericGrMod aliases -----------------------------

pub type GenericChainComplex<I, R> = ChainComplex<I, GenericKey<I>, R>;
pub type GenericChainComplex1<R> = GenericChainComplex<isize,  R>;
pub type GenericChainComplex2<R> = GenericChainComplex<isize2, R>;
pub type GenericChainComplex3<R> = GenericChainComplex<isize3, R>;

pub type GenericGrMod<I, R> = GrMod<I, GenericKey<I>, R>;
pub type GenericGrMod1<R> = GenericGrMod<isize,  R>;
pub type GenericGrMod2<R> = GenericGrMod<isize2, R>;
pub type GenericGrMod3<R> = GenericGrMod<isize3, R>;

impl<I, R> GenericChainComplex<I, R>
where I: AddInd, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_d_matrices(d_deg: I, matrices: impl IntoIterator<Item = (I, SpMat<R>)>) -> Self {
        let d_matrices: HashMap<I, SpMat<R>> = matrices.into_iter().collect();

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

    pub fn dual(&self) -> Self {
        let d_deg = self.d_deg();
        let matrices = self.support().map(|&i| {
            let prev = i - d_deg;
            let m = if self.is_supported(prev) {
                self.d_matrix(prev).transpose()
            } else {
                SpMat::zero((0, self[i].rank()))
            };
            (i, m)
        }).collect::<Vec<_>>();
        Self::from_d_matrices(-d_deg, matrices)
    }
}
