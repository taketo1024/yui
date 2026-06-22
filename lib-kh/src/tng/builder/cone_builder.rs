//! Builds the involutive cone `Cone(1 + τ)` at the *cobordism* level — before delooping — for a
//! strongly invertible link. [`ConeBuilder`] owns a [`SymTngBuilder`], drives the symmetric build,
//! and returns the reduced cone as a [`TngComplex`] (the `KhGen → KhIGen` conversion is done by the
//! caller in `khi`). Char-2 only.
//!
//! Two paths, by `chunks`:
//! - **whole-complex**: build the symmetric complex, double it (`make_cone`), then eliminate the
//!   invertible `1`-edges and deloop.
//! - **incremental** (`cone_merge`): fuse the final chunk merge with cone construction + reduction,
//!   so the full un-delooped product is never materialized.
//!
//! Reference:
//! - T. Sano, "Involutive Khovanov homology and equivariant knots",
//!   Algebr. Geom. Topol. 25 (2025), 5059–5111.

use itertools::Itertools;
use num_traits::Zero;
use yui_core::bitseq::Bit;
use yui_core::{Ring, RingOps};
use yui_link::{Edge, Node, InvLink};

use crate::tng::{Cob, CobComp, LcCob, Tng, TngComplex, TngComplexKey, TngComplexVertex};
use super::{reachable_range, ChunkBuilder, SymTngBuilder, SymBuildConfig, TngComplexBuilder, BuildConfig, BuildMode, TauKeyMap};

pub struct ConeBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    inner: SymTngBuilder<R>,
    cone: TngComplexBuilder<R>, // the reduced cone, filled in by the final `cone_merge`
}

impl<R> ConeBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_inv_link(l: &InvLink, h: &R, t: &R, reduced: bool) -> Self {
        let inner = SymTngBuilder::from_inv_link(l, h, t, reduced);
        let cone = TngComplexBuilder::init(h, t, (0, 0), None); // replaced by `cone_merge`
        Self { inner, cone }
    }

    pub fn with_config(mut self, config: SymBuildConfig) -> Self {
        self.inner = self.inner.with_config(config);
        self
    }

    // Build the reduced chunks, merge all but the last normally, and close with the incremental
    // cone merge. The non-chunked case is just a degenerate plan (see `plan`).
    pub fn run(mut self) -> Self {
        let plan = self.plan();
        let mode = self.inner.config().mode;
        let chunks = ChunkBuilder { builder: &self.inner }.build_chunks(plan);
        let last = chunks.len().saturating_sub(1);

        for (i, (chunk, (c, key_map, elems))) in chunks.into_iter().enumerate() {
            self.inner.drop_nodes(|x| chunk.contains(x));
            if i == last {
                self.cone_merge(c, key_map, mode);
            } else {
                self.inner.merge(c, key_map, elems);
            }
        }
        self
    }

    pub fn into_tng_complex(self) -> TngComplex<R> {
        self.cone.into_tng_complex()
    }

    // The cutwidth partition when `chunks` is set, else a degenerate split of the first τ-unit off
    // the rest — so even the non-chunked knot is closed by one incremental `cone_merge`.
    fn plan(&self) -> Vec<Vec<Node>> {
        if self.inner.config().chunks.is_some() || self.inner.config().cut.is_some() {
            return ChunkBuilder { builder: &self.inner }.plan();
        }
        let x = self.inner.choose_next_node().expect("a crossing to split off").clone();
        let tx = self.inner.inv_node(&x).clone();
        let unit = if tx == x { vec![x] } else { vec![x, tx] };
        let rest: Vec<Node> = self.inner.nodes().iter().filter(|n| !unit.contains(n)).cloned().collect();
        if rest.is_empty() { vec![unit] } else { vec![unit, rest] }
    }

    // Fuse the final merge with cone construction + reduction. Per symmetric degree: merge the slice,
    // cone-ify it, then collapse the invertible `1`-edges two degrees behind (so `cone_extend` has
    // moved past the degree being removed). Delooping is deferred to the end, on the small survivors.
    fn cone_merge(&mut self, other: TngComplex<R>, other_map: TauKeyMap, mode: BuildMode) {
        let left_map = std::mem::take(self.inner.key_map_mut());
        let (left, right) = self.inner.complex_mut().prepare_merge(other);
        let range = reachable_range(self.inner.complex().h_range(), &self.inner.config().h_range, self.inner.n_nodes());
        self.cone = TngComplexBuilder::from_tng_complex(cone_shell(self.inner.complex()), BuildConfig { mode, ..Default::default() });

        let (start, top) = (*range.start(), *range.end());
        for d in range {
            self.inner.merge_slice(&left, &right, d, &left_map, &other_map);
            self.cone_extend(d);
            if d - 2 >= start {
                self.cone.eliminate_in(d - 2);
            }
            if d > start {
                self.prune(d - 1);
            }
        }
        for d in (top - 1) ..= (top + 1) {
            self.cone.eliminate_in(d);
        }
        self.deloop_all(false);
        self.deloop_all(true);
    }

    // Add the symmetric complex's degree-`d` slice to the cone: the two copies `k·0`, `k·1` of each
    // vertex, the within-layer edges *into* degree `d`, and the `1+τ` edges out of `k·0`. Processed
    // ascending, every edge lands exactly once (its target's degree).
    fn cone_extend(&mut self, d: isize) {
        let with_bit = |k: &TngComplexKey, b: Bit| {
            let mut key = *k;
            key.state.push(b);
            key
        };
        let keys = self.inner.complex().keys_of_deg(d).copied().collect_vec();

        for k in keys.iter() {
            let tng = self.inner.complex().vertex(k).tng().clone();
            self.cone.complex_mut().add_vertex(with_bit(k, Bit::Bit0), TngComplexVertex::from(tng.clone()));
            self.cone.complex_mut().add_vertex(with_bit(k, Bit::Bit1), TngComplexVertex::from(tng));
        }

        // within-layer edges into degree d (each layer copies the symmetric differential).
        for k in keys.iter() {
            for j in self.inner.complex().vertex(k).in_edges().copied().collect_vec() {
                let f = self.inner.complex().edge(&j, k).clone();
                self.cone.complex_mut().add_edge(&with_bit(&j, Bit::Bit0), &with_bit(k, Bit::Bit0), f.clone());
                self.cone.complex_mut().add_edge(&with_bit(&j, Bit::Bit1), &with_bit(k, Bit::Bit1), f);
            }
        }

        // connecting differential (1 + τ): k·0 → k·1 (id) and k·0 → τk·1 (τ-cyl).
        for k in keys.iter() {
            let tng = self.inner.complex().vertex(k).tng().clone();
            let id = LcCob::from(Cob::id(&tng));
            let tau = LcCob::from(tau_cob(&tng, |e| self.inner.inv_edge(e)));
            let tk = *self.inner.key_map().inv_key(k);
            let k0 = with_bit(k, Bit::Bit0);

            if &tk == k {
                let f = id + tau; // same target: sum over char 2
                if !f.is_zero() {
                    self.cone.complex_mut().add_edge(&k0, &with_bit(k, Bit::Bit1), f);
                }
            } else {
                self.cone.complex_mut().add_edge(&k0, &with_bit(k, Bit::Bit1), id);
                self.cone.complex_mut().add_edge(&k0, &with_bit(&tk, Bit::Bit1), tau);
            }
        }
    }

    // Deloop every degree (marked circles only when `based`), eliminating to finish the reduction.
    fn deloop_all(&mut self, based: bool) {
        let auto_elim = self.cone.config().mode.auto_elim();
        for d in self.cone.complex().h_range() {
            self.cone.deloop_in_with(d, based);
        }
        if auto_elim {
            for d in self.cone.complex().h_range() {
                self.cone.eliminate_in(d);
            }
        }
    }

    // Drop a consumed symmetric degree and its τ key-map entries — never needed again.
    fn prune(&mut self, d: isize) {
        let doomed = self.inner.complex().keys_of_deg(d).copied().collect_vec();
        self.inner.complex_mut().remove_vertices(&doomed);
        let i0 = self.inner.complex().deg_shift().0;
        self.inner.key_map_mut().drop(|k| k.weight() as isize + i0 == d);
    }
}

// ---- cobordism-level cone construction (`1 + τ`, char-2) ----

// An empty cone shell: same `deg_shift`/base point, one extra h-degree for the cone bit.
fn cone_shell<R>(c: &TngComplex<R>) -> TngComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (h, t) = c.ht();
    TngComplex::new(h, t, c.deg_shift(), c.base_pt(), c.dim() + 1, Default::default())
}

// The τ morphism on a closed tangle: one relabel cylinder `c → τc` per circle.
fn tau_cob<G>(tng: &Tng, inv_edge: G) -> Cob
where G: Fn(Edge) -> Edge {
    Cob::new(tng.comps().map(|c|
        CobComp::plain(Tng::from(c.clone()), Tng::from(c.convert_edges(&inv_edge)))
    ))
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use yui_core::num::FF2;
    use yui_link::InvLink;
    use super::*;

    // Build the reduced cone, assert d² = 0, and return its nonzero homology ranks per degree.
    // (Full homology vs. the KhI reference is checked in `khi` via the `KhGen → KhIGen` extraction.)
    fn cone_homology(l: &InvLink, reduced: bool, config: SymBuildConfig) -> Vec<(isize, usize)> {
        let c = ConeBuilder::from_inv_link(l, &FF2::zero(), &FF2::zero(), reduced)
            .with_config(config).run().into_tng_complex().into_raw_complex();
        c.check_d_all();
        let h = c.homology();
        h.support().map(|&i| (i, h[i].rank())).filter(|(_, r)| *r > 0).sorted().collect()
    }

    // The cone homology must not depend on the chunking: whole == chunked, reduced and unreduced.
    fn check_chunk_independent(l: &InvLink, chunks: usize) {
        for reduced in [false, true] {
            let whole = cone_homology(l, reduced, SymBuildConfig::default());
            let chunked = cone_homology(l, reduced, SymBuildConfig { chunks: Some(chunks), ..Default::default() });
            assert_eq!(whole, chunked, "reduced={reduced}");
        }
    }

    #[test]
    fn cone_chunk_independent_3_1() {
        check_chunk_independent(&InvLink::test_data("3_1"), 2);
    }

    #[test]
    fn cone_chunk_independent_4_1() {
        check_chunk_independent(&InvLink::test_data("4_1"), 2);
    }

    #[test]
    fn cone_chunk_independent_6_3() {
        check_chunk_independent(&InvLink::test_data("6_3"), 3);
    }

    // The cone homology must not depend on the simplification mode.
    #[test]
    fn cone_mode_independent() {
        let l = InvLink::test_data("6_3");
        let reference = cone_homology(&l, false, SymBuildConfig::default());
        for mode in [BuildMode::MinFill, BuildMode::NoElim, BuildMode::None] {
            let h = cone_homology(&l, false, SymBuildConfig { mode, ..Default::default() });
            assert_eq!(h, reference, "mode {mode:?}");
        }
    }
}
