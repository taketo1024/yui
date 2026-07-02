use std::collections::HashMap;

use itertools::Itertools;
use yui_core::lc::{LcKey, Lc};
use yui_core::{Ring, RingOps};
use yui_homology::{isize2, GrMod1, GrMod2, Summand};

/// q-graded decomposition of a `GrMod1` into a `GrMod2`.
///
/// Implementors describe the underlying `GrMod1` (`base`) and how to extract
/// a chain's q-degree (`decomp_key`); the default `bigraded()` computes a
/// fresh `GrMod2`. Each type typically also provides an inherent
/// `cached_bigraded()` that caches the result in a `OnceLock` field.
pub(crate) trait Bigraded<X, R>
where X: LcKey, R: Ring, for<'x> &'x R: RingOps<R> {
    fn base(&self) -> &GrMod1<X, R>;
    fn decomp_key(&self, z: &Lc<X, R>) -> isize;

    fn bigraded(&self) -> GrMod2<X, R> {
        decomp_by(self.base(), |z| self.decomp_key(z))
    }

    // structurally equal as bigraded modules: same (free rank, torsion) at every bidegree.
    fn is_identical<B>(&self, other: &B) -> bool
    where B: Bigraded<X, R> {
        let strip = |base: &GrMod1<X, R>, g: &dyn Fn(&Lc<X, R>) -> isize| -> HashMap<isize2, (usize, Vec<R>)> {
            decomp_info(base, g).into_iter().map(|(k, (r, t, _))| (k, (r, t))).collect()
        };
        strip(self.base(), &|z| self.decomp_key(z)) == strip(other.base(), &|z| other.decomp_key(z))
    }
}

pub(crate) fn decomp_info<X, R, F>(grid: &GrMod1<X, R>, decomp_key: F) -> HashMap<isize2, (usize, Vec<R>, Vec<usize>)>
where X: LcKey, R: Ring, for<'x> &'x R: RingOps<R>, F: Fn(&Lc<X, R>) -> isize {
    let mut table = HashMap::new();
    let init_entry = (0, vec![], vec![]);

    for (&i, h) in grid.iter() {
        let r = h.rank();
        let t = h.tors().len();

        for k in 0..r + t {
            let z = h.generator(k);
            let q = decomp_key(&z);
            let e = table.entry(isize2(i, q)).or_insert_with(|| init_entry.clone());
            if k < r {
                e.0 += 1;
            } else {
                e.1.push(h.tors()[k - r].clone());
            }
            e.2.push(k);
        }
    }

    table
}

pub(crate) fn decomp_by<X, R, F>(grid: &GrMod1<X, R>, decomp_key: F) -> GrMod2<X, R>
where X: LcKey, R: Ring, for<'x> &'x R: RingOps<R>, F: Fn(&Lc<X, R>) -> isize {
    let info = decomp_info(grid, decomp_key);
    let support = info.keys().cloned().collect_vec();

    GrMod2::generate(support, move |idx| {
        let i = idx.0;
        let Some(e) = info.get(&idx) else {
            return Summand::zero()
        };

        let (rank, tors, indices) = e;
        let gens = grid[i].raw_generators().clone();
        let trans = grid[i].trans().sub(indices);
        Summand::new(gens, *rank, tors.clone(), trans)
    })
}
