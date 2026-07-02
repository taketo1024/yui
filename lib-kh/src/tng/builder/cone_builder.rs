//! Builds the involutive cone `Cone(1 + τ)` at the *cobordism* level — before delooping — for a
//! strongly invertible link. [`ConeBuilder`] owns a [`SymTngBuilder`], drives the symmetric build,
//! and converts the reduced cone at the boundary: `into_raw_complex` yields the matrix-backed
//! `KhIGen`-keyed chain complex, `eval_khi_elements` the canon classes. Char-2 only.
//!
//! The build always closes with the incremental `cone_merge`: the final chunk merge is fused with
//! cone construction + reduction per symmetric degree, so the full un-delooped product is never
//! materialized (a non-chunked build uses the degenerate two-piece plan — see `plan`).
//!
//! Reference:
//! - T. Sano, "Involutive Khovanov homology and equivariant knots",
//!   Algebr. Geom. Topol. 25 (2025), 5059–5111.

use delegate::delegate;
use itertools::Itertools;
use log::{debug, info};
use num_traits::Zero;
use yui_core::bitseq::Bit;
use yui_core::{Ring, RingOps};
use yui_link::{Edge, Node, InvLink};

use rustc_hash::FxHashMap;
use yui_homology::{ChainComplex1, GrMod1, Summand};
use yui_matrix::sparse::SpMat;

use crate::kh::{KhChain, KhGen};
use crate::khi::{KhIChain, KhIGen, KhIGenExt};
use crate::tng::{Cob, CobComp, LcCob, LcCobTrait, Tng, TngComplex, TngComplexElem, TngComplexKey, TngComplexVertex};
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
        info!("build config:\n{:#?}", self.inner.config());
        info!("cutwidth profile:\n{}", self.inner.profile_sym());

        let plan = self.plan();
        let mode = self.inner.config().mode;

        info!("cone build: {} chunks, mode: {mode:?}", plan.len());

        let chunks = ChunkBuilder { builder: &self.inner }.build_chunks(plan);
        let last = chunks.len().saturating_sub(1);

        for (i, (chunk, (c, key_map, elems))) in chunks.into_iter().enumerate() {
            info!("cone chunk {}/{}: {} nodes{}", i + 1, last + 1, chunk.len(), if i == last { " (cone merge)" } else { "" });
            self.inner.drop_nodes(|x| chunk.contains(x));
            if i == last {
                self.cone_merge(c, key_map, elems, mode);
            } else {
                self.inner.merge(c, key_map, elems);
            }
        }
        self
    }

    delegate! {
        to self.cone {
            pub fn into_tng_complex(self) -> TngComplex<R>;
            pub fn eval_elements(&self) -> Vec<KhChain<R>>;
        }
    }

    // The cutwidth partition when `chunks` is set, else a degenerate split of the first τ-unit off
    // the rest — so even the non-chunked knot is closed by one incremental `cone_merge`.
    fn plan(&self) -> Vec<Vec<Node>> {
        if self.inner.config().cut.enabled() {
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
    fn cone_merge(&mut self, other: TngComplex<R>, other_map: TauKeyMap, other_elems: Vec<TngComplexElem<R>>, mode: BuildMode) {
        let left_map = std::mem::take(self.inner.key_map_mut());
        let (left, right) = self.inner.complex_mut().prepare_merge(other);
        let range = reachable_range(self.inner.complex().h_range(), &self.inner.config().h_range, self.inner.n_nodes());
        self.cone = TngComplexBuilder::from_tng_complex(cone_shell(self.inner.complex()), BuildConfig { mode, ..Default::default() });

        info!("cone merge {} <- {}, range: {range:?}", left.stat(), right.stat());

        // complete the symmetric canon cycles, then seed their `B`/`Q` (bit-0/bit-1) copies into the
        // cone so its deloop/eliminate carry them. `cone_extend` will create the referenced vertices.
        self.inner.elements_mut().merge(other_elems);
        self.seed_cone_elements();

        let (start, top) = (*range.start(), *range.end());
        for d in range {
            debug!("cone merge C[{d}]...");
            self.inner.merge_slice(&left, &right, d, &left_map, &other_map);
            self.cone_extend(d);
            if d - 2 >= start {
                self.cone.eliminate_in(d - 2);
            }
            if d > start {
                self.prune_consumed(d - 1);
            }
            debug!("  cone C[{d}] done: {}", self.cone.stat());
        }

        // eliminate_in(top) collapses the id-pairing C[top] (Bit0) → C[top+1] (Bit1), shrinking the
        // top before deloop; C[top+1]'s remainder is dropped by `prune_isolated_top` below.
        // (top+1 itself needs no eliminate_in: its vertices have no out-edges.)
        info!("cone eliminate top C[{}..={}]: {}", top - 1, top, self.cone.stat());
        for d in (top - 1) ..= top {
            self.cone.eliminate_in(d);
        }

        self.prune_isolated_top(top);

        info!("cone deloop: {}", self.cone.stat());
        self.deloop_all(false);
        self.deloop_all(true);

        info!("cone done: {}", self.cone.stat());
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

        debug!("  cone-extend C[{d}]: +{} verts", 2 * keys.len());

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

    // Lift each completed symmetric canon cycle to its `B` (bit-0) and `Q` (bit-1) cone copies.
    fn seed_cone_elements(&mut self) {
        let lifted = self.inner.elements().content().iter().flat_map(|e|
            [lift_elem(e, Bit::Bit0), lift_elem(e, Bit::Bit1)]
        ).collect_vec();

        debug!("  seed {} cone elements", lifted.len());

        self.cone.elements_mut().set(lifted);
    }

    // Deloop every degree (marked circles only when `based`), eliminating to finish the reduction.
    fn deloop_all(&mut self, based: bool) {
        debug!("cone deloop-all (based: {based})...");

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
    fn prune_consumed(&mut self, d: isize) {
        let doomed = self.inner.complex().keys_of_deg(d).copied().collect_vec();

        debug!("  prune sym C[{d}]: -{} verts", doomed.len());

        self.inner.complex_mut().remove_vertices(&doomed);
        let i0 = self.inner.complex().deg_shift().0;
        self.inner.key_map_mut().drop(|k| k.weight() as isize + i0 == d);
    }

    // Restrict a truncated cone to what KhI within the window needs: drop the out-of-window layer
    // C[top+1] (the `Bit1` shift of C[top]), then prune no-in-edge vertices at C[top] — they fed
    // only the dropped layer. No-op for a full-range build (its C[top+1] is genuine).
    fn prune_isolated_top(&mut self, top: isize) {
        if self.inner.config().h_range.is_none() {
            return;
        }
        let doomed = self.cone.complex().keys_of_deg(top + 1).copied().collect_vec();
        debug!("cone: drop out-of-window C[{}] ({} verts)", top + 1, doomed.len());
        self.cone.complex_mut().remove_vertices(&doomed);

        let total = self.cone.complex().keys_of_deg(top).count();
        let doomed = self.cone.complex().keys_of_deg(top)
            .filter(|k| self.cone.complex().vertex(k).in_edges().next().is_none())
            .copied()
            .collect_vec();
        debug!("cone: window top C[{top}] {}/{total} no-in-edge verts pruned", doomed.len());
        self.cone.complex_mut().remove_vertices(&doomed);
    }

    /// The canon classes as KhI chains (cone bit → `B`/`Q`).
    pub fn eval_khi_elements(&self) -> Vec<KhIChain<R>> {
        self.eval_elements().into_iter()
            .map(|z| z.map_keys(|x| into_khi_gen(&x)))
            .collect()
    }

    /// Convert the reduced cone directly into the KhI chain complex: generators become `KhIGen`s
    /// (cone bit → `B`/`Q`), the scalar edges are materialized as cached differential matrices, and
    /// the cobordism-level complex is dropped — downstream reduction is pure linear algebra.
    pub fn into_raw_complex(self) -> ChainComplex1<KhIGen, R> {
        let c = self.cone.into_tng_complex();
        assert!(c.is_completely_delooped());

        let (h, t) = c.ht().clone();

        // per degree: keys sorted by their KhIGen q-degree (descending), fixing both the summand
        // generator order and the matrix row/column order.
        let keys: FxHashMap<isize, Vec<TngComplexKey>> = c.h_range().map(|i| {
            let sorted = c.keys_of_deg(i).copied()
                .sorted_by_key(|k| -into_khi_gen(&k.as_gen()).rel_q_deg())
                .collect_vec();
            (i, sorted)
        }).collect();

        let summands = GrMod1::generate(c.h_range(), |i| {
            Summand::from_raw_generators(keys[&i].iter().map(|k| into_khi_gen(&k.as_gen())))
        });

        let matrices = c.h_range().map(|i| {
            let cols = &keys[&i];
            let rows: FxHashMap<&TngComplexKey, usize> = keys.get(&(i + 1))
                .map(|ks| ks.iter().enumerate().map(|(r, k)| (k, r)).collect())
                .unwrap_or_default();

            let entries = cols.iter().enumerate().flat_map(|(j, k)| {
                let (c, rows, h, t) = (&c, &rows, &h, &t);
                c.vertex(k).out_edges().map(move |l|
                    (rows[l], j, c.edge(k, l).eval(h, t))
                )
            });

            let m = SpMat::from_entries((rows.len(), cols.len()), entries);
            (i, m)
        }).collect_vec();

        drop(c);
        ChainComplex1::new_with_d_matrices(summands, 1, matrices)
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

// Cone key → KhIGen: the cone bit is the last `state` bit (`0` → `B`/`Left`, `1` → `Q`/`Right`);
// stripping it gives the underlying `KhGen` (with the true q-grading).
fn into_khi_gen(x: &KhGen) -> KhIGen {
    let mut state = *x.state();
    let bit = state.iter().last().unwrap();
    state.remove(state.len() - 1);
    let under = KhGen::new(state, *x.tensor());
    match bit {
        Bit::Bit0 => KhIGen::from_left(under),
        Bit::Bit1 => KhIGen::from_right(under),
    }
}

// Lift a symmetric canon cycle to the cone layer `bit`: push the cone bit onto each `out_cob` key's
// state (`Bit0` = `B`, `Bit1` = `Q`); `in_cob`/`state` are unchanged.
fn lift_elem<R>(e: &TngComplexElem<R>, bit: Bit) -> TngComplexElem<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let mut e = e.clone();
    let out = std::mem::take(e.out_cob_mut());
    *e.out_cob_mut() = out.into_iter().map(|(mut k, f)| {
        k.state.push(bit);
        (k, f)
    }).collect();
    e
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use yui_core::num::FF2;
    use yui_link::InvLink;
    use super::*;
    use super::super::CutOption;

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
            let chunked = cone_homology(l, reduced, SymBuildConfig { cut: CutOption::Auto(chunks), ..Default::default() });
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

    #[test]
    fn cone_chunk_windowed_9_46() {
        let l = InvLink::from_symmetric_pd_code([[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]);
        let narrow = |h: Vec<(isize, usize)>| h.into_iter().filter(|&(d, _)| d <= 0).collect_vec();
        for reduced in [false, true] {
            let full = narrow(cone_homology(&l, reduced, SymBuildConfig::default()));
            let chunked = narrow(cone_homology(&l, reduced, SymBuildConfig { cut: CutOption::Auto(2), mode: BuildMode::MinFill, h_range: Some(-64 ..= 1), ..Default::default() }));
            assert_eq!(full, chunked, "reduced={reduced}");
        }
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

    // The cone's canon classes must give the same ssi as the matrix cone.
    #[test]
    fn cone_canon_ssi_matches_matrix() {
        use yui_core::poly::Poly;
        use crate::khi::{KhIHomology, ssi_invariants};
        use crate::util::calc::div_vec;

        type P = Poly<'H', FF2>;
        let (c, t) = (P::variable(), P::zero());
        let knots = [
            ("3_1", InvLink::test_data("3_1")),
            ("4_1", InvLink::test_data("4_1")),
            ("6_3", InvLink::test_data("6_3")),
            // 9_46 has s̲ ≠ s̄ (ssi = (0, 2)) — exercises the canon-cycle ordering.
            ("9_46", InvLink::from_symmetric_pd_code([[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]])),
        ];
        for (name, l) in knots {
            let matrix = ssi_invariants(&l, &c, false);

            let config = SymBuildConfig { h_range: Some(isize::MIN + 1 ..= 1), ..Default::default() };
            let kh = KhIHomology::from_cone(&l, &c, &t, false, config);
            let zs = kh.canon_cycles();
            assert_eq!(zs.len(), 4);

            let ds = zs.iter().map(|z| {
                let h = kh.h_deg_of_chain(z);
                div_vec(&kh[h].vectorize_euc(z).subvec(0..2), &c).expect("invalid divisibility")
            }).collect_vec();

            let (w, r) = (l.writhe(), l.seifert_circles().len() as i32);
            let cone = (2 * ds[0] + w - r + 1, 2 * ds[2] + w - r + 1);
            assert_eq!(cone, matrix, "{name}");
        }
    }
}
