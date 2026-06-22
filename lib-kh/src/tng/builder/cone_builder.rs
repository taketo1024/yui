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
use yui_link::{Edge, InvLink};

use crate::tng::{Cob, CobComp, LcCob, Tng, TngComplex, TngComplexKey, TngComplexVertex};
use super::{reachable_range, ChunkBuilder, SymTngBuilder, SymBuildConfig, TngComplexBuilder, BuildConfig, BuildMode, TauKeyMap};

pub struct ConeBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    inner: SymTngBuilder<R>,
    cone: Option<TngComplex<R>>,
}

impl<R> ConeBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_inv_link(l: &InvLink, h: &R, t: &R, reduced: bool) -> Self {
        Self { inner: SymTngBuilder::from_inv_link(l, h, t, reduced), cone: None }
    }

    pub fn with_config(mut self, config: SymBuildConfig) -> Self {
        self.inner = self.inner.with_config(config);
        self
    }

    pub fn run(mut self) -> Self {
        let cone = if self.inner.config().chunks.is_some() {
            self.run_chunked()
        } else {
            self.run_whole()
        };
        self.cone = Some(cone);
        self
    }

    pub fn into_tng_complex(self) -> TngComplex<R> {
        self.cone.expect("`run` must be called before `into_tng_complex`")
    }

    // Non-chunked: build the symmetric complex, double it whole, then reduce.
    fn run_whole(&mut self) -> TngComplex<R> {
        if self.inner.config().preprocess {
            self.inner.preprocess();
        }
        self.inner.process_nodes();

        let mode = self.inner.config().mode;
        let cone = {
            let key_map = self.inner.key_map().clone();
            let e_map = self.inner.e_map().clone();
            make_cone(self.inner.complex(), &|k| *key_map.inv_key(k), &|e| e_map[&e])
        };
        reduce(cone, mode)
    }

    // Chunked: build the reduced chunks, merge all but the last normally, and close with the
    // incremental cone merge.
    fn run_chunked(&mut self) -> TngComplex<R> {
        let mode = self.inner.config().mode;
        let e_map = self.inner.e_map().clone();
        let inv_edge = |e: Edge| e_map[&e];
        let chunks = ChunkBuilder { builder: &self.inner }.build_chunks();
        let last = chunks.len().saturating_sub(1);

        let mut cone = None;
        for (i, (chunk, (c, key_map, elems))) in chunks.into_iter().enumerate() {
            self.inner.drop_nodes(|x| chunk.contains(x));
            if i == last {
                cone = Some(self.cone_merge(c, key_map, &inv_edge, mode));
            } else {
                self.inner.merge(c, key_map, elems);
            }
        }
        cone.expect("at least one chunk")
    }

    // Fuse the final merge with cone construction + reduction. Per symmetric degree: merge the slice,
    // cone-ify it, then collapse the invertible `1`-edges two degrees behind (so `cone_extend` has
    // moved past the degree being removed). Delooping is deferred to the end, on the small survivors.
    fn cone_merge<G>(&mut self, other: TngComplex<R>, other_map: TauKeyMap, inv_edge: &G, mode: BuildMode) -> TngComplex<R>
    where G: Fn(Edge) -> Edge {
        let left_map = std::mem::take(self.inner.key_map_mut());
        let (left, right) = self.inner.complex_mut().prepare_merge(other);
        let range = reachable_range(self.inner.complex().h_range(), &self.inner.config().h_range, self.inner.n_nodes());
        let mut cone = TngComplexBuilder::from_tng_complex(cone_shell(self.inner.complex()), BuildConfig { mode, ..Default::default() });

        let (start, top) = (*range.start(), *range.end());
        for d in range {
            self.inner.merge_slice(&left, &right, d, &left_map, &other_map);
            {
                let key_map = self.inner.key_map();
                let tau = |k: &TngComplexKey| *key_map.inv_key(k);
                cone_extend(self.inner.complex(), cone.complex_mut(), d, &tau, inv_edge);
            }
            if d - 2 >= start {
                cone.eliminate_in(d - 2);
            }
            if d > start {
                self.prune(d - 1);
            }
        }
        for d in (top - 1) ..= (top + 1) {
            cone.eliminate_in(d);
        }
        deloop_all(&mut cone, mode, false);
        deloop_all(&mut cone, mode, true);
        cone.into_tng_complex()
    }

    // Drop a consumed symmetric degree and its τ key-map entries — never needed again.
    fn prune(&mut self, d: isize) {
        let doomed = self.inner.complex().keys_of_deg(d).copied().collect_vec();
        self.inner.complex_mut().remove_vertices(&doomed);
        let i0 = self.inner.complex().deg_shift().0;
        self.inner.key_map_mut().drop(|k| k.weight() as isize + i0 == d);
    }
}

// Whole-complex reduction: collapse invertible `1`-edges first, then deloop (MinFill re-eliminates).
fn reduce<R>(cone: TngComplex<R>, mode: BuildMode) -> TngComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let mut b = TngComplexBuilder::from_tng_complex(cone, BuildConfig { mode, ..Default::default() });
    b.eliminate_all();
    let mut b = b.run();
    if mode.auto_elim() && !mode.immediate_elim() {
        b.eliminate_all();
    }
    b.into_tng_complex()
}

// Deloop every degree (marked circles only when `based`), eliminating to finish the reduction.
fn deloop_all<R>(cone: &mut TngComplexBuilder<R>, mode: BuildMode, based: bool)
where R: Ring, for<'x> &'x R: RingOps<R> {
    for d in cone.complex().h_range() {
        cone.deloop_in_with(d, based);
    }
    if mode.auto_elim() {
        for d in cone.complex().h_range() {
            cone.eliminate_in(d);
        }
    }
}

// ---- cobordism-level cone construction ----

// Cobordism-level mapping cone of `1 + τ` for a closed complex (the involutive complex). Doubles
// every vertex — a cone bit appended to its key — and connects the two layers by `1` (identity) and
// `τ` (per-circle relabel cylinder). Char-2, so no cone signs.
fn make_cone<R, F, G>(c: &TngComplex<R>, tau_key: &F, inv_edge: &G) -> TngComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R>, F: Fn(&TngComplexKey) -> TngComplexKey, G: Fn(Edge) -> Edge {
    debug_assert!(c.is_closed());
    let mut cone = cone_shell(c);
    for d in c.h_range() {
        cone_extend(c, &mut cone, d, tau_key, inv_edge);
    }
    cone
}

// An empty cone shell: same `deg_shift`/base point, one extra h-degree for the cone bit.
fn cone_shell<R>(c: &TngComplex<R>) -> TngComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (h, t) = c.ht();
    TngComplex::new(h, t, c.deg_shift(), c.base_pt(), c.dim() + 1, Default::default())
}

// Add `c`'s degree-`d` contribution to `cone`: the two copies `k·0`, `k·1` of each degree-`d`
// vertex, the within-layer edges *into* degree `d`, and the `1+τ` edges out of `k·0`. Processed
// ascending, every edge lands exactly once (its target's degree).
fn cone_extend<R, F, G>(c: &TngComplex<R>, cone: &mut TngComplex<R>, d: isize, tau_key: &F, inv_edge: &G)
where R: Ring, for<'x> &'x R: RingOps<R>, F: Fn(&TngComplexKey) -> TngComplexKey, G: Fn(Edge) -> Edge {
    let with_bit = |k: &TngComplexKey, b: Bit| {
        let mut key = *k;
        key.state.push(b);
        key
    };
    let keys = c.keys_of_deg(d).copied().collect_vec();

    for k in keys.iter() {
        let tng = c.vertex(k).tng().clone();
        cone.add_vertex(with_bit(k, Bit::Bit0), TngComplexVertex::from(tng.clone()));
        cone.add_vertex(with_bit(k, Bit::Bit1), TngComplexVertex::from(tng));
    }

    // within-layer edges into degree d (each layer copies the symmetric differential).
    for k in keys.iter() {
        for j in c.vertex(k).in_edges().copied().collect_vec() {
            let f = c.edge(&j, k).clone();
            cone.add_edge(&with_bit(&j, Bit::Bit0), &with_bit(k, Bit::Bit0), f.clone());
            cone.add_edge(&with_bit(&j, Bit::Bit1), &with_bit(k, Bit::Bit1), f);
        }
    }

    // connecting differential (1 + τ): k·0 → k·1 (id) and k·0 → τk·1 (τ-cyl).
    for k in keys.iter() {
        let v = c.vertex(k);
        let k0 = with_bit(k, Bit::Bit0);
        let id = LcCob::from(Cob::id(v.tng()));
        let tau = LcCob::from(tau_cob(v.tng(), inv_edge));
        let tk = tau_key(k);

        if &tk == k {
            let f = id + tau; // same target: sum over char 2
            if !f.is_zero() {
                cone.add_edge(&k0, &with_bit(k, Bit::Bit1), f);
            }
        } else {
            cone.add_edge(&k0, &with_bit(k, Bit::Bit1), id);
            cone.add_edge(&k0, &with_bit(&tk, Bit::Bit1), tau);
        }
    }
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

    // The reduced cone must be a valid complex (d² = 0) — full homology vs. the reference is checked
    // in `khi`, once the `KhGen → KhIGen` extraction is in place.
    fn check_valid(l: &InvLink, chunks: Option<usize>) {
        for reduced in [false, true] {
            let config = SymBuildConfig { chunks, ..Default::default() };
            let cone = ConeBuilder::from_inv_link(l, &FF2::zero(), &FF2::zero(), reduced)
                .with_config(config).run().into_tng_complex();
            cone.into_raw_complex().check_d_all();
        }
    }

    #[test]
    fn cone_3_1_whole() {
        check_valid(&InvLink::test_data("3_1"), None);
    }

    #[test]
    fn cone_3_1_chunked() {
        check_valid(&InvLink::test_data("3_1"), Some(2));
    }
}
