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

use crate::kh::{KhChain, KhGen, KhTensor};
use crate::khi::{KhIChain, KhIGen, KhIGenExt};
use crate::tng::{Cob, CobComp, Dot, End, LcCob, LcCobTrait, Tng, TngComp, TngComplex, TngComplexElem, TngComplexKey, TngComplexVertex};
use super::{reachable_range, ChunkBuilder, SymTngBuilder, SymBuildConfig, TngComplexBuilder, BuildConfig, BuildMode, TauKeyMap};

// τ-orbit class of a symmetric-complex vertex: the representative (smaller key) and the τ-fixed
// survive into the direct cone; the non-representative column is never materialized.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum OrbitClass {
    Fixed, Rep, Drop
}

pub struct ConeBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    inner: SymTngBuilder<R>,
    cone: TngComplexBuilder<R>, // the reduced cone, filled in by the final `cone_merge`
    // vertical identity pivots `k·0 → k·1` per free τ-orbit, recorded at `cone_extend` (while the
    // τ-pairing is still alive) and consumed by `vertical_reduce` two degrees behind.
    pending_vertical: FxHashMap<isize, Vec<(TngComplexKey, TngComplexKey)>>,
}

impl<R> ConeBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_inv_link(l: &InvLink, h: &R, t: &R, reduced: bool) -> Self {
        let inner = SymTngBuilder::from_inv_link(l, h, t, reduced);
        let cone = TngComplexBuilder::init(h, t, (0, 0), None); // replaced by `cone_merge`
        Self { inner, cone, pending_vertical: FxHashMap::default() }
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
            if self.direct() {
                self.rewrite_elements(d - 1); // degree d-1's out-edges are now complete
            }
            if d - 2 >= start {
                self.vertical_reduce(d - 2);
                self.cone.eliminate_in(d - 2); // stragglers: τ-fixed `1+τ` units, correction-created units
            }
            // direct emission reads the dropped column's in-edges one degree deeper, so its prune lags one more.
            if self.direct() {
                if d - 1 > start {
                    self.prune_consumed(d - 2);
                }
            } else if d > start {
                self.prune_consumed(d - 1);
            }
            debug!("  cone C[{d}] done: {}", self.cone.stat());
        }

        // eliminate_in(top) collapses the id-pairing C[top] (Bit0) → C[top+1] (Bit1), shrinking the
        // top before deloop; C[top+1]'s remainder is dropped by `prune_isolated_top` below.
        // (top+1 itself needs no eliminate_in: its vertices have no out-edges.)
        if self.direct() {
            self.rewrite_elements(top); // the top degree has no further out-edges — Iτ pushes only
        }

        info!("cone eliminate top C[{}..={}]: {}", top - 1, top, self.cone.stat());
        for d in (top - 1) ..= top {
            self.vertical_reduce(d);
            self.cone.eliminate_in(d);
        }

        // prune only when the window truncates the KhI complex: if the requested top exceeds the
        // reachable Kh top, C[top+1] is the genuine KhI top degree (Bit1 of C[top]) and must stay.
        let window_top = self.inner.config().h_range.as_ref().map(|r| *r.end());
        if window_top.is_some_and(|t| t <= top) {
            self.prune_isolated_top(top);
        }

        info!("cone done: {}", self.cone.stat());
    }

    fn direct(&self) -> bool {
        self.inner.config().cone_direct
    }

    fn cone_extend(&mut self, d: isize) {
        if self.direct() {
            self.cone_extend_direct(d);
        } else {
            self.cone_extend_full(d);
        }
    }

    // Symmetry-breaking direct emission (Sano2026, Prop 4.6): per free τ-orbit only the
    // representative's two copies are created; the dropped column's contributions enter as
    // closed-form corrections. Corrections through two dropped corners vanish (an edge into a
    // dropped `l⁰` is an in-edge of a removed pivot source; one out of a dropped `j¹` is an
    // out-edge of a removed pivot target), so the rewrite is one elimination step deep:
    //   surv → surv : verbatim on both layers,
    //   surv → drop : (ii)  j¹ → (τk)¹ += Iτ(k)∘f,
    //   drop → surv : (iii) (τj)⁰ → k⁰ += f∘Iτ(τj),
    //   cross via dropped N at d−1 : (i) X¹ → k⁰ += (X→N)∘(N→k),
    // and only τ-fixed vertices keep a vertical `1 + τ` (a representative's identity cancels via τ²).
    fn cone_extend_direct(&mut self, d: isize) {
        let with_bit = |k: &TngComplexKey, b: Bit| {
            let mut key = *k;
            key.state.push(b);
            key
        };
        let (h, t) = self.cone.complex().ht().clone();
        let keys = self.inner.complex().keys_of_deg(d).copied().collect_vec();
        let dropped = keys.iter().filter(|k| self.orbit_class(k).1 == OrbitClass::Drop).count();

        debug!("  cone-extend-direct C[{d}]: +{} verts ({dropped} dropped)", 2 * (keys.len() - dropped));

        // vertices: representatives and τ-fixed only.
        for k in keys.iter() {
            if self.orbit_class(k).1 == OrbitClass::Drop {
                continue;
            }
            let tng = self.inner.complex().vertex(k).tng().clone();
            self.cone.complex_mut().add_vertex(with_bit(k, Bit::Bit0), TngComplexVertex::from(tng.clone()));
            self.cone.complex_mut().add_vertex(with_bit(k, Bit::Bit1), TngComplexVertex::from(tng));
        }

        // within-layer edges into degree d, rewritten per the tables above.
        for k in keys.iter() {
            let (tk, ck) = self.orbit_class(k);
            for j in self.inner.complex().vertex(k).in_edges().copied().collect_vec() {
                let f = self.inner.complex().edge(&j, k).clone();
                let (tj, cj) = self.orbit_class(&j);
                match (cj, ck) {
                    (OrbitClass::Drop, OrbitClass::Drop) => {}
                    (OrbitClass::Drop, _) => {
                        let itau = LcCob::from(tau_cob(self.inner.complex().vertex(&tj).tng(), |e| self.inner.inv_edge(e)));
                        let corr = itau.stack(&f).reduce(&h, &t);
                        self.cone.complex_mut().add_to_edge(&with_bit(&tj, Bit::Bit0), &with_bit(k, Bit::Bit0), corr);
                    }
                    (_, OrbitClass::Drop) => {
                        let itau = LcCob::from(tau_cob(self.inner.complex().vertex(k).tng(), |e| self.inner.inv_edge(e)));
                        let corr = f.stack(&itau).reduce(&h, &t);
                        self.cone.complex_mut().add_to_edge(&with_bit(&j, Bit::Bit1), &with_bit(&tk, Bit::Bit1), corr);
                    }
                    _ => {
                        self.cone.complex_mut().add_to_edge(&with_bit(&j, Bit::Bit0), &with_bit(k, Bit::Bit0), f.clone());
                        self.cone.complex_mut().add_to_edge(&with_bit(&j, Bit::Bit1), &with_bit(k, Bit::Bit1), f);
                    }
                }
            }
        }

        // (i) cross corrections through each dropped N at degree d−1.
        let dropped_prev = self.inner.complex().keys_of_deg(d - 1)
            .filter(|n| self.orbit_class(n).1 == OrbitClass::Drop)
            .copied().collect_vec();

        for n in dropped_prev {
            let ins = self.inner.complex().vertex(&n).in_edges().copied()
                .filter(|x| self.orbit_class(x).1 != OrbitClass::Drop).collect_vec();
            let outs = self.inner.complex().vertex(&n).out_edges().copied()
                .filter(|k| self.orbit_class(k).1 != OrbitClass::Drop).collect_vec();

            for x in ins.iter() {
                let b = self.inner.complex().edge(x, &n).clone();
                for k in outs.iter() {
                    let corr = b.stack(self.inner.complex().edge(&n, k)).reduce(&h, &t);
                    self.cone.complex_mut().add_to_edge(&with_bit(x, Bit::Bit1), &with_bit(k, Bit::Bit0), corr);
                }
            }
        }

        // verticals: τ-fixed only.
        for k in keys.iter() {
            if self.orbit_class(k).1 != OrbitClass::Fixed {
                continue;
            }
            let tng = self.inner.complex().vertex(k).tng().clone();
            let id = LcCob::from(Cob::id(&tng));
            let tau = LcCob::from(tau_cob(&tng, |e| self.inner.inv_edge(e)));
            let f = id + tau; // same target: sum over char 2
            if !f.is_zero() {
                self.cone.complex_mut().add_edge(&with_bit(k, Bit::Bit0), &with_bit(k, Bit::Bit1), f);
            }
        }
    }

    // Direct mode: retract canon-element components off the dropped column of degree `d`
    // (Sano2026, Prop 4.6 SDR): an entry at `N·0` is the pivot's dependent coordinate and drops;
    // an entry at `N·1` redirects to `(τN)·1` via Iτ and to `l·0` via each out-edge `N → l`
    // (mirroring `eliminate_from` with the identity pivot). Corrections landing on a dropped
    // `l·0` are removed by the next degree's rewrite, matching the sequential SDR composition.
    fn rewrite_elements(&mut self, d: isize) {
        let with_bit = |k: &TngComplexKey, b: Bit| {
            let mut key = *k;
            key.state.push(b);
            key
        };

        let dropped = self.inner.complex().keys_of_deg(d)
            .filter(|n| self.orbit_class(n).1 == OrbitClass::Drop)
            .copied().collect_vec();

        if dropped.is_empty() {
            return;
        }

        let (h, t) = self.cone.complex().ht().clone();

        let pushes = dropped.iter().map(|n| {
            let tn = self.orbit_class(n).0;
            let itau = LcCob::from(tau_cob(self.inner.complex().vertex(n).tng(), |e| self.inner.inv_edge(e)));
            let mut outs = vec![(with_bit(&tn, Bit::Bit1), itau)];
            for l in self.inner.complex().vertex(n).out_edges() {
                outs.push((with_bit(l, Bit::Bit0), self.inner.complex().edge(n, l).clone()));
            }
            (with_bit(n, Bit::Bit0), with_bit(n, Bit::Bit1), outs)
        }).collect_vec();

        let elems = self.cone.elements_mut().take().into_iter().map(|mut e| {
            e.modify_out_cob(|mut cob| {
                for (n0, n1, outs) in pushes.iter() {
                    cob.remove(n0);
                    let Some(f) = cob.remove(n1) else { continue };
                    for (target, edge) in outs.iter() {
                        let corr = f.stack(edge).reduce(&h, &t);
                        let s = if let Some(g) = cob.remove(target) { g - corr } else { -corr };
                        if !s.is_zero() {
                            cob.insert(*target, s);
                        }
                    }
                }
                cob
            });
            e
        }).collect_vec();

        self.cone.elements_mut().set(elems);
    }

    // τ-orbit classification via the inner key map; the representative is the smaller key.
    fn orbit_class(&self, k: &TngComplexKey) -> (TngComplexKey, OrbitClass) {
        let tk = *self.inner.key_map().inv_key(k);
        let class = if tk == *k {
            OrbitClass::Fixed
        } else if *k < tk {
            OrbitClass::Rep
        } else {
            OrbitClass::Drop
        };
        (tk, class)
    }

    // Add the symmetric complex's degree-`d` slice to the cone: the two copies `k·0`, `k·1` of each
    // vertex, the within-layer edges *into* degree `d`, and the `1+τ` edges out of `k·0`. Processed
    // ascending, every edge lands exactly once (its target's degree).
    fn cone_extend_full(&mut self, d: isize) {
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

                // record the vertical pivot at the orbit's non-representative (Sano2026, Prop 4.6);
                // the representative's identity cancels via the τ² correction when this eliminates.
                if *k > tk {
                    self.pending_vertical.entry(d).or_default().push((k0, with_bit(k, Bit::Bit1)));
                }
            }
        }
    }

    // Symmetry-breaking reduction (Sano2026, Prop 4.6): eliminate the recorded vertical identity
    // pivots of degree `d` directly, skipping the pivot search. The recorded pairs are hints, not
    // guarantees — the straggler `eliminate_in` passes may have consumed a vertex or spoiled the
    // edge's invertibility, and such leftovers fall back to the generic pass.
    fn vertical_reduce(&mut self, d: isize) {
        let Some(pairs) = self.pending_vertical.remove(&d) else { return };

        debug!("  vertical-reduce C[{d}]: {} orbit pairs", pairs.len());

        for (v, w) in pairs {
            let c = self.cone.complex();
            let valid = c.contains_key(&v) && c.contains_key(&w)
                && c.has_edge(&v, &w)
                && c.edge(&v, &w).is_invertible();
            if valid {
                self.cone.eliminate(&v, &w);
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

    /// The canon classes as KhI chains (cone bit → `B`/`Q`), expanded over each vertex's circles
    /// by the same pairing as `into_raw_complex`.
    pub fn eval_khi_elements(&self) -> Vec<KhIChain<R>> {
        let c = self.cone.complex();
        let (h, t) = c.ht().clone();

        self.cone.elements().content().iter().map(|e| {
            let init = LcCob::from(e.in_cob().clone());
            e.out_cob().iter().flat_map(|(k, retr)| {
                let circles = circles_of(c.vertex(k).tng());
                label_assignments(&circles).into_iter().map(|b| {
                    let g = cap_circles(retr.clone(), End::Tgt, &circles, &b, &h, &t);
                    let x = (g * &init).eval(&h, &t);
                    (into_khi_gen(&expanded_key(k, &b).as_gen()), x)
                }).collect_vec()
            }).collect()
        }).collect()
    }

    /// Convert the reduced cone directly into the KhI chain complex, without delooping: each
    /// closed vertex expands into its label assignments (one `KhAlgGen` per circle; a marked
    /// circle is fixed to `X`), and each edge contributes the scalars `⟨b|f|a⟩` — `f` with source
    /// circles cupped by `a` and target circles capped by `b` (the deloop pairing), evaluated as
    /// a closed cobordism. Downstream reduction is pure linear algebra.
    pub fn into_raw_complex(self) -> ChainComplex1<KhIGen, R> {
        let c = self.cone.into_tng_complex();
        let (h, t) = c.ht().clone();

        info!("build raw complex: {}", c.stat());

        // expanded generators (vertex, labels) per degree, sorted by their KhIGen q-degree
        // (descending) — fixing both the summand generator order and the matrix row/column order.
        let keys: FxHashMap<isize, Vec<(TngComplexKey, KhTensor)>> = c.h_range().map(|i| {
            let expanded = c.keys_of_deg(i).flat_map(|k| {
                let circles = circles_of(c.vertex(k).tng());
                label_assignments(&circles).into_iter().map(|a| (*k, a)).collect_vec()
            }).sorted_by_key(|(k, a)|
                -into_khi_gen(&expanded_key(k, a).as_gen()).rel_q_deg()
            ).collect_vec();
            (i, expanded)
        }).collect();

        let summands = GrMod1::generate(c.h_range(), |i| {
            Summand::from_raw_generators(keys[&i].iter().map(|(k, a)| into_khi_gen(&expanded_key(k, a).as_gen())))
        });

        let matrices = c.h_range().map(|i| {
            let cols = &keys[&i];
            let rows: FxHashMap<TngComplexKey, usize> = keys.get(&(i + 1))
                .map(|ks| ks.iter().enumerate().map(|(r, (k, a))| (expanded_key(k, a), r)).collect())
                .unwrap_or_default();

            let mut entries = vec![];
            for (j, (k, a)) in cols.iter().enumerate() {
                let src_circles = circles_of(c.vertex(k).tng());
                for l in c.vertex(k).out_edges() {
                    let fa = cap_circles(c.edge(k, l).clone(), End::Src, &src_circles, a, &h, &t);
                    if fa.is_zero() {
                        continue;
                    }
                    let tgt_circles = circles_of(c.vertex(l).tng());
                    for b in label_assignments(&tgt_circles) {
                        let g = cap_circles(fa.clone(), End::Tgt, &tgt_circles, &b, &h, &t);
                        entries.push((rows[&expanded_key(l, &b)], j, g.eval(&h, &t)));
                    }
                }
            }

            debug!("  raw d[{i}]: {} -> {}, nnz: {}", cols.len(), rows.len(), entries.len());
            let m = SpMat::from_entries((rows.len(), cols.len()), entries);
            (i, m)
        }).collect_vec();

        info!("raw complex done: {} gens", keys.values().map(|ks| ks.len()).sum::<usize>());

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

// The circles of a closed vertex in expansion order: unmarked first (matching the deloop order),
// the marked circle last.
fn circles_of(tng: &Tng) -> Vec<TngComp> {
    debug_assert!(tng.comps().all(|c| c.is_circle()));
    let (marked, unmarked): (Vec<_>, Vec<_>) = tng.comps().cloned().partition(|c| c.is_marked());
    unmarked.into_iter().chain(marked).collect()
}

// All label assignments for the circles (a marked circle is fixed to `X`).
fn label_assignments(circles: &[TngComp]) -> Vec<KhTensor> {
    KhTensor::generate(circles.len()).filter(|a|
        circles.iter().enumerate().all(|(i, c)| !c.is_marked() || a[i].is_X())
    ).collect()
}

// The expanded key: the vertex key with the assignment appended to its label.
fn expanded_key(k: &TngComplexKey, a: &KhTensor) -> TngComplexKey {
    let mut kk = *k;
    kk.label.append(*a);
    kk
}

// Cap the given circles by the deloop pairing: a source circle labeled `X` is cupped with `Dot::X`
// (`1` plain); a target circle labeled `X` is capped plain (`1` with `Dot::Y`).
fn cap_circles<R>(f: LcCob<R>, e: End, circles: &[TngComp], labels: &KhTensor, h: &R, t: &R) -> LcCob<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    circles.iter().enumerate().fold(f, |f, (i, c)| {
        let dot = match (e, labels[i].is_X()) {
            (End::Src, true)  => Some(Dot::X),
            (End::Src, false) => None,
            (End::Tgt, true)  => None,
            (End::Tgt, false) => Some(Dot::Y),
        };
        f.cap_off(e, c, dot).reduce(h, t)
    })
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
    e.modify_out_cob(|out| out.into_iter().map(|(mut k, f)| {
        k.state.push(bit);
        (k, f)
    }).collect());
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
    // (Full homology vs. the KhI reference is checked in `khi`.)
    fn cone_homology(l: &InvLink, reduced: bool, config: SymBuildConfig) -> Vec<(isize, usize)> {
        let c = ConeBuilder::from_inv_link(l, &FF2::zero(), &FF2::zero(), reduced)
            .with_config(config).run().into_raw_complex();
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

    // The direct symmetry-broken emission (Sano2026, Prop 4.6) is a deformation retract of the
    // doubled cone: homology must agree with the double-then-eliminate path.
    fn check_direct_matches(l: &InvLink, config: SymBuildConfig) {
        for reduced in [false, true] {
            let full = cone_homology(l, reduced, config.clone());
            let direct = cone_homology(l, reduced, SymBuildConfig { cone_direct: true, ..config.clone() });
            assert_eq!(full, direct, "reduced={reduced}");
        }
    }

    #[test]
    fn cone_direct_3_1() {
        check_direct_matches(&InvLink::test_data("3_1"), SymBuildConfig::default());
    }

    #[test]
    fn cone_direct_3_1_m() {
        check_direct_matches(&InvLink::test_data("3_1").mirror(), SymBuildConfig::default());
    }

    #[test]
    fn cone_direct_4_1() {
        check_direct_matches(&InvLink::test_data("4_1"), SymBuildConfig::default());
    }

    #[test]
    fn cone_direct_6_3_chunked() {
        check_direct_matches(&InvLink::test_data("6_3"), SymBuildConfig { cut: CutOption::Auto(3), ..Default::default() });
    }

    #[test]
    fn cone_direct_9_46_windowed() {
        let l = InvLink::from_symmetric_pd_code([[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]);
        let config = SymBuildConfig { cut: CutOption::Auto(2), mode: BuildMode::MinFill, h_range: Some(-64 ..= 1), ..Default::default() };
        for reduced in [false, true] {
            let full = cone_homology(&l, reduced, SymBuildConfig { h_range: Some(-64 ..= 1), ..Default::default() });
            let direct = cone_homology(&l, reduced, SymBuildConfig { cone_direct: true, ..config.clone() });
            let narrow = |h: Vec<(isize, usize)>| h.into_iter().filter(|&(d, _)| d <= 0).collect_vec();
            assert_eq!(narrow(full), narrow(direct), "reduced={reduced}");
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
