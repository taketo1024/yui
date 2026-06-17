//! Equivariant variant of [`TngComplexBuilder`] for strongly invertible links:
//! processes axis-symmetric crossings singly and off-axis crossings in pairs
//! `(x, τx)`, then maintains a `key_map` matching each [`TngComplexKey`] with
//! its τ-image so that [`SymTngBuilder::tau_map`] can be used to assemble the
//! involutive Khovanov complex `CKhI = Cone(1 + τ)`.
//!
//! Reference:
//! - T. Sano, "Involutive Khovanov homology and equivariant knots",
//!   Algebr. Geom. Topol. 25 (2025), 5059–5111.
//!   <https://doi.org/10.2140/agt.2025.25.5059>, <https://arxiv.org/abs/2404.08568>

use std::collections::HashSet;
use std::ops::RangeInclusive;
use rustc_hash::{FxHashMap, FxHashSet};
use itertools::Itertools;
use log::{debug, info};
use yui_core::algo::KeyedUnionFind;
use yui_core::bitseq::{Bit, BitSeq};
use yui_core::{Ring, RingOps};
use yui_link::{Node, Edge, InvLink};

use crate::kh::{KhGen, KhTensor};
use crate::tng::{End, LcCobTrait, TngComp, TngComplex, TngComplexElem, TngComplexKey};
use crate::tng::builder::{TngComplexBuilder, BuildConfig, BuildMode, NodeOrder};
use super::{reachable_range, pop_min_pivot};

/// Toggles for the automatic simplification done while building (kept separate
/// from [`BuildConfig`] so the equivariant builder can gain its own flags).
#[derive(Clone, Debug)]
pub struct SymBuildConfig {
    // crossing order: LoopGreedy (default) or MinCut (bounds cutwidth for wide knots).
    pub node: NodeOrder,
    pub mode: BuildMode,
    // build half the off-axis crossings and mirror via τ (see `preprocess`).
    pub preprocess: bool,
    // chooser handicap on an off-axis pair = coeff · current size (its extra growth).
    pub pair_penalty_coeff: f64,
    // divide-and-conquer: build the link in chunks of ≤ this many crossings (each via a
    // child builder), merging each reduced chunk into the parent. None = single pass.
    pub chunk_bound: Option<usize>,
    // literal truncation: homology at the endpoints is wrong (build `(a-1)..=(b+1)` for correct `[a, b]`).
    pub h_range: Option<RangeInclusive<isize>>,
}

impl Default for SymBuildConfig {
    fn default() -> Self {
        Self { node: NodeOrder::default(), mode: BuildMode::default(), preprocess: true, pair_penalty_coeff: 1.0, chunk_bound: None, h_range: None }
    }
}

// τ-symmetric key map: each `TngComplexKey`'s τ-image. On-axis keys (`τk = k`) are a set;
// off-axis keys form an involution stored both ways for O(1) `inv_key`. The merge iterates
// only one representative per off-axis pair (the map is symmetric), then symmetrizes.
#[derive(Clone, Default)]
struct TauKeyMap {
    on_axis: FxHashSet<TngComplexKey>,
    off_axis: FxHashMap<TngComplexKey, TngComplexKey>,
}

impl FromIterator<(TngComplexKey, TngComplexKey)> for TauKeyMap {
    fn from_iter<I: IntoIterator<Item = (TngComplexKey, TngComplexKey)>>(iter: I) -> Self {
        let mut m = Self::default();
        for (k, tk) in iter {
            m.add_pair(k, tk);
        }
        m
    }
}

impl TauKeyMap {
    fn init() -> Self {
        Self::from_iter([(TngComplexKey::init(), TngComplexKey::init())])
    }

    fn len(&self) -> usize {
        self.on_axis.len() + self.off_axis.len()
    }

    fn inv_key(&self, k: &TngComplexKey) -> &TngComplexKey {
        self.on_axis.get(k).unwrap_or_else(|| &self.off_axis[k])
    }

    fn is_sym(&self, k: &TngComplexKey) -> bool {
        self.on_axis.contains(k)
    }

    fn add_pair(&mut self, k: TngComplexKey, tk: TngComplexKey) {
        if k == tk {
            self.on_axis.insert(k);
        } else if !self.off_axis.contains_key(&k) {
            self.off_axis.insert(k, tk);
            self.off_axis.insert(tk, k);
        }
    }

    fn remove(&mut self, k: &TngComplexKey) {
        if self.on_axis.remove(k) {
            return;
        }
        let tk = self.off_axis.remove(k).unwrap();
        self.off_axis.remove(&tk);
    }

    fn keys(&self) -> impl Iterator<Item = &TngComplexKey> + '_ {
        self.on_axis.iter().chain(self.off_axis.keys())
    }

    // τ preserves weight, so if `pred` drops a key it drops its mirror too — symmetric.
    fn drop(&mut self, pred: impl Fn(&TngComplexKey) -> bool) {
        self.on_axis.retain(|k| !pred(k));
        self.off_axis.retain(|k, _| !pred(k));
    }

    // Key map of a half-complex (`keys`) tensored with its τ-mirror: `k1+k2 ↦ k2+k1`
    // (τ swaps the halves), within `band`.
    fn from_half(keys: &[TngComplexKey], band: RangeInclusive<usize>) -> Self {
        keys.iter().flat_map(|&k1| {
            let band = band.clone();
            keys.iter().filter_map(move |&k2| {
                band.contains(&(k1.weight() + k2.weight())).then(|| (k1 + k2, k2 + k1))
            })
        }).collect()
    }
}

pub struct SymTngBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    inner: TngComplexBuilder<R>,
    x_map: FxHashMap<Node, Node>,
    e_map: FxHashMap<Edge, Edge>,
    key_map: TauKeyMap,
    config: SymBuildConfig,
    real_top: isize, // n₊ = the complex's max h-degree (deg_shift.0 + #crossings)
}

impl<R> SymTngBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_inv_link(l: &InvLink, h: &R, t: &R, reduced: bool) -> SymTngBuilder<R> {
        assert!(l.nodes().all(|x| x.is_crossing()));
        assert!(!reduced || l.base_pt().is_some());

        // the inner builder is driven by `self` — disable its own auto-simplify.
        let inner = TngComplexBuilder::from_link(l.inner(), h, t, reduced)
            .with_config(BuildConfig { mode: BuildMode::None, node: NodeOrder::default(), h_range: None });

        let x_map = l.nodes().map(|x|
            (x.clone(), l.inv_node(x).clone())
        ).collect();
        let e_map = l.edges().into_iter().map(|e| (e, l.inv_edge(e))).collect();
        let key_map = TauKeyMap::init();
        let real_top = inner.complex().deg_shift().0 + l.inner().n_crossings() as isize;

        SymTngBuilder { inner, x_map, e_map, key_map, config: SymBuildConfig::default(), real_top }
    }

    pub fn with_config(mut self, config: SymBuildConfig) -> Self {
        // propagate the window to the inner builder so the preprocess merges cap
        // to it; this also drops canon cycles when the window excludes h-degree 0.
        let inner_config = BuildConfig { mode: BuildMode::None, node: config.node, h_range: config.h_range.clone() };
        self.inner = self.inner.with_config(inner_config);
        self.config = config;
        self
    }

    pub fn run(mut self) -> Self {
        if self.config.chunk_bound.is_some() {
            self.process_chunks();
        } else {
            if self.config.preprocess {
                self.preprocess();
            }
            self.process_nodes();
        }
        self.finalize();
        self
    }

    // Divide-and-conquer: `ChunkBuilder` selects + builds each reduced chunk, then merge it into
    // the parent (so the parent never materializes the full dense slice).
    fn process_chunks(&mut self) {
        // canon cycles aren't yet tracked through chunk merges (see `merge`'s TODO), so
        // drop them — a chunked build yields the complex/homology but not α / ssi.
        self.inner.set_elements(vec![]);
        while let Some(chunk) = ChunkBuilder::new(self).next_chunk() {
            info!("process chunk ({}): {}", chunk.len(), chunk.iter().join(", "));

            let (c, key_map) = ChunkBuilder::new(self).build_chunk(&chunk);

            debug!("  chunk built: {}", c.stat());

            self.inner.drop_nodes(|x| chunk.contains(x));
            self.merge(c, key_map);

            info!("  chunk merged: {}", self.inner.stat());
        }
    }

    fn preprocess(&mut self) {
        SymTngPreprocessor::run(self);
    }

    fn process_nodes(&mut self) {
        info!("process {} nodes", self.inner.nodes().count());

        while let Some(x) = self.choose_next_node_sym().cloned() {
            let tx = self.inv_node(&x).clone();
            if x == tx {
                self.append_on_axis(&x);
            } else {
                self.append_off_axis(&x, &tx);
            }

            info!("  appended: {x}, current size: {}", self.inner.stat());
        }
    }

    /// Pick the next τ-unit by node order, ties broken by earliest crossing order. MinCut scores the
    /// *combined* x+τx toggle (a shared axis edge cancels — can't be summed); LoopGreedy sums loops.
    fn choose_next_node_sym(&self) -> Option<&Node> {
        let pair_penalty = (self.config.pair_penalty_coeff * self.inner.complex().n_verts() as f64) as isize;
        self.inner.nodes().enumerate()
            .min_by_key(|(i, x)| {
                let tx = self.inv_node(x);
                let score = match self.config.node {
                    NodeOrder::LoopGreedy => {
                        if tx == *x { self.inner.loop_count(x) }
                        else { self.inner.loop_count(x) + self.inner.loop_count(tx) - pair_penalty }
                    },
                    NodeOrder::MinCut => {
                        let mut edges = x.edges().to_vec();
                        if tx != *x { edges.extend_from_slice(tx.edges()); }
                        -self.inner.cutwidth_of(edges)
                    },
                };
                (-score, *i)
            })
            .map(|(_, x)| x)
    }

    fn append_on_axis(&mut self, x: &Node) { 
        info!("({}/{}) append on-axis: {x}", 
            self.inner.complex().dim() + 1, 
            self.inner.complex().dim() + self.inner.nodes().count(), 
        );

        self.inner.prepare_append(x);

        let (h, t) = self.inner.complex().ht();
        let c = TngComplex::from_node(h, t, x, self.inner.complex().base_pt());
        let key_map = if x.is_crossing() {
            [Bit::Bit0, Bit::Bit1].map(|b| { 
                let k = TngComplexKey { state: BitSeq::from(b), label: KhTensor::empty() };
                (k, k)
            }).into_iter().collect()
        } else { 
            let k = TngComplexKey::init();
            [(k, k)].into_iter().collect()
        };

        self.merge(c, key_map);
    }

    fn append_off_axis(&mut self, x: &Node, tx: &Node) { 
        assert_eq!(self.inv_node(x), tx);

        info!("({}/{}) append off-axis: {x}, {tx}", 
            self.inner.complex().dim() + 1, 
            self.inner.complex().dim() + self.inner.nodes().count(), 
        );

        self.inner.prepare_append(x);
        self.inner.prepare_append(tx);

        let c = {
            let (h, t) = self.inner.complex().ht();
            let mut c = TngComplex::from_node(h, t, x, self.inner.complex().base_pt());
            c.append_node(tx);
            c
        };
        let key_map = if x.is_crossing() { 
            [
                ([0, 0], [0, 0]),
                ([1, 0], [0, 1]),
                ([0, 1], [1, 0]),
                ([1, 1], [1, 1])
            ].map(|(b0, b1)| { 
                let k = TngComplexKey { state: BitSeq::from_iter(b0), label: KhTensor::empty() };
                let l = TngComplexKey { state: BitSeq::from_iter(b1), label: KhTensor::empty() };
                (k, l)
            }).into_iter().collect()
        } else {
            let k = TngComplexKey::init();
            [(k, k)].into_iter().collect()
        };

        self.merge(c, key_map);
    }

    fn merge(&mut self, c: TngComplex<R>, right_map: TauKeyMap) {
        debug!("merge {} + {}", self.inner.stat(), c.stat());
        debug!("  key_map: {} × {}", self.key_map.len(), right_map.len());

        // build the merged τ key-map per degree (next to merge_vertices) rather than as one
        // up-front cartesian — for large knots that product never fits in memory.
        let left_map = std::mem::take(&mut self.key_map);

        let (left, right) = self.inner.complex_mut().prepare_merge(c);
        let range = reachable_range(self.inner.complex().h_range(), &self.config.h_range, self.inner.nodes().count());
        debug!("  merge range: {:?}", range);

        match self.config.mode {
            BuildMode::None    => for i in range { self.merge_slice(&left, &right, i, &left_map, &right_map); },
            BuildMode::MinFill => self.merge_deferred(&left, &right, range, &left_map, &right_map),
            _                  => self.merge_immediate(&left, &right, range, &left_map, &right_map),
        }

        self.prune_h_range();

        // TODO merge elements

        debug!("  key_map built: {}", self.key_map.len());
        debug!("  merged: {} + {} -> {}", left.stat(), right.stat(), self.inner.stat());
    }

    // Immediate: per degree, eliminate then deloop (deloop inline-eliminates each new vertex). Deloop
    // may be selective (defer non-productive circles, then full-deloop + re-pass to a fixpoint at the end).
    fn merge_immediate(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, range: RangeInclusive<isize>, left_map: &TauKeyMap, right_map: &TauKeyMap) {
        let top = *range.end();
        let selective = self.config.mode.is_selective();

        for i in range {
            debug!("build C[{i}]...");
            self.merge_slice(left, right, i, left_map, right_map);
            self.eliminate_in(i - 1);
            self.deloop_in(i - 1, false, selective);
            debug!("  built C[{i}]: {}", self.inner.complex().rank(i));
        }

        self.prune_isolated_top(top);
        self.deloop_in(top, false, selective);

        // re-run selective to a fixpoint (catch loops turned productive by later equiv elims), then full-deloop the rest.
        if selective {
            self.deloop_selective();
            self.deloop_all(false, false);
        }
    }

    // Repeatedly deloop productive circles until none remain — each equiv elimination can turn a
    // previously-deferred loop productive.
    fn deloop_selective(&mut self) {
        let mut step = 0;
        loop {
            let before = self.inner.complex().n_verts();
            debug!("selective re-pass {step}: start ({before} verts)");
            self.deloop_all(false, true);
            let after = self.inner.complex().n_verts();
            debug!("  selective re-pass {step}: {before} -> {after} verts (diff {})",
                after as isize - before as isize);
            if after == before { break }
            step += 1;
        }
    }

    // MinFill: per degree, deloop then eliminate i-1, i by global equivariant min-fill (incremental
    // Markowitz). Always full-deloops — combining it with selective delooping only balloons the transient.
    fn merge_deferred(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, range: RangeInclusive<isize>, left_map: &TauKeyMap, right_map: &TauKeyMap) {
        let top = *range.end();

        for i in range {
            debug!("build C[{i}]...");
            self.merge_slice(left, right, i, left_map, right_map);
            self.deloop_in(i - 1, false, false);
            self.eliminate_in(i - 2);
            self.eliminate_in(i - 1);
            debug!("  built C[{i}]: {}", self.inner.complex().rank(i));
        }

        self.prune_isolated_top(top);
        self.deloop_in(top, false, false);
        self.eliminate_in(top - 1);
        // no eliminate_in(top): top has no outgoing edges, so it would be a no-op.
    }

    // Build degree `i`: the τ key-map slice, then its vertices and the edges into it.
    fn merge_slice(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, i: isize, left_map: &TauKeyMap, right_map: &TauKeyMap) {
        self.merge_key_slice(left, right, i, left_map, right_map);
        self.inner.merge_slice(left, right, i);
    }

    // Degree-i slice of the merged τ key-map: k1+k2 ↦ τk1+τk2 (τ preserves degree,
    // so the slice is self-contained — see `merge`).
    fn merge_key_slice(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, i: isize, left_map: &TauKeyMap, right_map: &TauKeyMap) {
        for (k1, k2) in TngComplex::collect_keys(left, right, i) {
            self.key_map.add_pair(k1 + k2, left_map.inv_key(k1) + right_map.inv_key(k2));
        }
    }

    // No-in-edge vertices at a TRUNCATED window-top (`top < real_top`) feed only the discarded
    // boundary — drop them and their τ-pairs (at the real top they're genuine generators).
    // The no-in-edge set is τ-closed, so `key_map` stays an involution.
    fn prune_isolated_top(&mut self, top: isize) {
        let truncated = self.config.h_range.as_ref().is_some_and(|w| top == *w.end()) && top < self.real_top;
        if !truncated { return; }

        let doomed_verts = self.inner.complex().keys_of_deg(top)
            .filter(|k| self.inner.complex().vertex(k).in_edges().next().is_none())
            .copied()
            .collect_vec();
        if !doomed_verts.is_empty() {
            debug!("prune {} isolated verts in C[{top}].", doomed_verts.len());
        }
        let doomed: FxHashSet<_> = doomed_verts.iter().copied().collect();
        self.inner.complex_mut().remove_vertices(&doomed_verts);
        self.key_map.drop(|k| doomed.contains(k));
    }

    // Combined off-axis weight that can still reach the window, given `r` pending crossings.
    fn weight_band(&self, r: usize) -> RangeInclusive<usize> {
        let s = self.inner.complex().deg_shift().0;
        let window = self.config.h_range.as_ref().map(|w| (*w.start() - s) ..= (*w.end() - s));
        let band = reachable_range(0 ..= isize::MAX, &window, r);
        (*band.start()).max(0) as usize ..= (*band.end()).max(0) as usize
    }

    /// Prune doomed vertices *and* their `key_map` entries. τ preserves weight,
    /// so `deg(k) == deg(τk)` — pruning is symmetric and `key_map` stays a valid
    /// involution. The `key_map` is pruned independently of the complex because
    /// the cartesian merge over-generates it past the vertex cap.
    fn prune_h_range(&mut self) {
        let Some(h_range) = self.config.h_range.clone() else { return };
        let r = self.inner.nodes().count() as isize;
        let i0 = self.inner.complex().deg_shift().0;

        // a vertex of degree `d` reaches `[d, d + r]`, so it stays relevant iff
        // `d ∈ [a - r, b]` — current degrees that can still land in `h_range`.
        let live = (*h_range.start() - r) ..= *h_range.end();
        let doomed = |k: &TngComplexKey|
            !live.contains(&(k.weight() as isize + i0));

        let doomed_verts = self.inner.complex().keys_of(&doomed).copied().collect_vec();
        if !doomed_verts.is_empty() {
            debug!("prune {} verts outside h_range.", doomed_verts.len());
        }
        self.inner.complex_mut().remove_vertices(&doomed_verts);

        self.key_map.drop(doomed);
    }

    fn deloop_all(&mut self, allow_based: bool, selective: bool) {
        for i in self.inner.complex().h_range() {
            self.deloop_in(i, allow_based, selective);
        }
    }

    // A circle whose no-dot cap yields an equiv-invertible edge (deloop+elim fires, no doubling).
    // Marked circles are skipped — they're delooped in the final full sweep.
    fn find_productive_loop(&self, k: &TngComplexKey) -> Option<usize> {
        let c = self.inner.complex();
        let v = c.vertex(k);
        v.tng().comps().enumerate()
            .filter(|(_, comp)| comp.is_circle() && !comp.is_marked())
            .find(|(_, comp)|
                v.out_edges().any(|l|
                    c.edge(k, l).is_invertible_after_cap(End::Src, comp) && self.is_equiv_edge(k, l))
                || v.in_edges().any(|j|
                    c.edge(j, k).is_invertible_after_cap(End::Tgt, comp) && self.is_equiv_edge(j, k)))
            .map(|(r, _)| r)
    }

    // Selective mode deloops only productive circles; otherwise any circle.
    fn choose_loop(&self, k: &TngComplexKey, allow_based: bool, selective: bool) -> Option<usize> {
        if selective {
            self.find_productive_loop(k)
        } else {
            self.inner.find_loop_in(k, allow_based, false)
        }
    }

    // Pivot cost for selection: off-axis pivots are eliminated in τ-pairs, so count ~2× the fill.
    fn pivot_weight(&self, k: &TngComplexKey) -> usize {
        let w = self.inner.complex().vertex(k).c_weight();
        if self.key_map.is_sym(k) { w } else { 2 * w }
    }

    // Markowitz cost of eliminating `k → l`; off-axis pivots eliminate in τ-pairs (~2× the fill).
    fn equiv_edge_weight(&self, k: &TngComplexKey, l: &TngComplexKey) -> usize {
        let w = self.inner.complex().edge_weight(k, l);
        if self.key_map.is_sym(k) { w } else { 2 * w }
    }

    // Least Markowitz cost to eliminate `k`, over its equiv-invertible incident edges — the
    // equivariant pivot priority (vs the cruder `pivot_weight`).
    fn equiv_elim_cost(&self, k: &TngComplexKey) -> usize {
        let v = self.inner.complex().vertex(k);
        let outs = v.out_edges().filter(|l| self.is_equiv_inv_edge(k, l)).map(|l| self.equiv_edge_weight(k, l));
        let ins = v.in_edges().filter(|j| self.is_equiv_inv_edge(j, k)).map(|j| self.equiv_edge_weight(j, k));
        outs.chain(ins).min().unwrap_or(0)
    }

    fn deloop_in(&mut self, i: isize, allow_based: bool, selective: bool) {
        let mut keys = self.inner.collect_keys(i,
            |k| self.choose_loop(k, allow_based, selective).is_some(),
            |k| self.pivot_weight(k),
        );
        if keys.is_empty() { return }

        debug!("deloop in C[{i}], targets: {} (selective: {selective}).", keys.len());

        let before = self.inner.complex().rank(i) as isize;

        while let Some(k) = pop_min_pivot(&mut keys, |k|
            self.inner.complex().contains_key(k).then(|| self.pivot_weight(k))
        ) {
            let Some(r) = self.choose_loop(&k, allow_based, selective) else { continue };

            for new_key in self.deloop_equiv(&k, r) {
                if self.choose_loop(&new_key, allow_based, selective).is_some() {
                    let w = self.pivot_weight(&new_key);
                    keys.push((new_key, w));
                }
            }
        }

        let after = self.inner.complex().rank(i) as isize;

        debug!("  delooped C[{i}]: {} (diff: {}).", after, after - before);
    }

    fn deloop_equiv(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> { 
        let mut added = if self.key_map.is_sym(k) {
            let c = self.inner.complex().vertex(k).tng().comp(r);
            if self.is_sym_comp(c) {
                // symmetric loop on symmetric key
                self.deloop_on_axis_sym(k, r)
            } else {
                // asymmetric loop on symmetric key
                self.deloop_on_axis_asym(k, r)
            }
        } else {
            // (symmetric or asymmetric) loop on asymmetric key
            self.deloop_off_axis(k, r)
        };

        // immediate elim eliminates each new vertex now; min-fill leaves them for the post-deloop
        // global pass, None leaves them entirely.
        if self.config.mode.immediate_elim() {
            added.retain(|k|
                self.inner.complex().contains_key(k) &&
                !self.try_eliminate_equiv_at(k)
            );
            // an equivariant elim removes the whole τ-pair, so a key kept above may since have
            // been eliminated as another's τ-partner — re-retain so only live keys are returned.
            added.retain(|k| self.inner.complex().contains_key(k));
        }
        added
    }

    fn deloop_on_axis_sym(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> {
        let c = self.inner.complex().vertex(k).tng().comp(r);

        debug_assert!(self.key_map.is_sym(k));
        debug_assert!(self.is_sym_comp(c));

        let updated = self.inner.deloop(k, r);

        self.key_map.remove(k);

        for &k_new in updated.iter() { 
            self.key_map.add_pair(k_new, k_new);
        }

        updated
    }

    #[allow(non_snake_case)]
    fn deloop_on_axis_asym(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> {
        let c = self.inner.complex().vertex(k).tng().comp(r);

        debug_assert!(self.key_map.is_sym(k));
        debug_assert!(!self.is_sym_comp(c));
        debug_assert!(!c.is_marked());

        //          ⚪︎1 | ⚪︎1
        //  ⚪︎1 | ⚪︎X  <-->  ⚪︎X | ⚪︎1
        //          ⚪︎X | ⚪︎X

        let tc = c.convert_edges(|e| self.inv_edge(e));

        let ks = self.inner.deloop(k, r);

        let (k_X, k_1) = (ks[0], ks[1]);
        let (k_XX, k_X1) = { 
            let tr = self.inner.complex().vertex(&k_X).tng().index_of(&tc).unwrap();
            let tks = self.inner.deloop(&k_X, tr);
            (tks[0], tks[1])
        };
        let (k_1X, k_11) = { 
            let tr = self.inner.complex().vertex(&k_1).tng().index_of(&tc).unwrap();
            let tks = self.inner.deloop(&k_1, tr);
            (tks[0], tks[1])
        };

        self.key_map.remove(k);

        self.key_map.add_pair(k_XX, k_XX);
        self.key_map.add_pair(k_X1, k_1X);
        self.key_map.add_pair(k_11, k_11);

        vec![k_XX, k_X1, k_1X, k_11]
    }

    #[allow(non_snake_case)]
    fn deloop_off_axis(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> {
        let c = self.inner.complex().vertex(k).tng().comp(r);

        debug_assert!(!self.key_map.is_sym(k));

        //  ⚪︎1 | ..  <-->  .. | ⚪︎1
        //  ⚪︎X | ..  <-->  .. | ⚪︎X

        let tk = *self.key_map.inv_key(k);
        let tc = c.convert_edges(|e| self.inv_edge(e));
        let tr = self.inner.complex().vertex(&tk).tng().index_of(&tc).unwrap();

        let mut ks = self.inner.deloop(k, r);
        let mut tks = self.inner.deloop(&tk, tr);

        self.key_map.remove(k);

        for (&k_new, &tk_new) in Iterator::zip(ks.iter(), tks.iter()) { 
            self.key_map.add_pair(k_new, tk_new);
        }

        ks.append(&mut tks);
        ks
    }

    fn eliminate_in(&mut self, i: isize) {
        let mut keys = self.inner.collect_keys(i,
            |k| self.inner.complex().vertex(k).out_edges()
                .any(|l| self.is_equiv_inv_edge(k, l)),
            |k| self.equiv_elim_cost(k),
        );
        if keys.is_empty() { return }

        debug!("eliminate in C[{i}], targets: {}", keys.len());

        let before = self.inner.complex().rank(i) as isize;

        while let Some(k) = pop_min_pivot(&mut keys, |k|
            self.inner.complex().contains_key(k).then(|| self.equiv_elim_cost(k))
        ) {
            self.try_eliminate_equiv_at(&k);
        }

        let after = self.inner.complex().rank(i) as isize;

        debug!("  eliminated C[{i}]: {} (diff: {}).", after, after - before);
    }

    fn choose_equiv_inv_edge_into(&self, k: &TngComplexKey) -> Option<&TngComplexKey> { 
        self.inner.complex().vertex(k).in_edges().filter_map(|j|
            self.is_equiv_inv_edge(j, k).then_some(j)
        )
        .min_by_key(|j| (self.inner.complex().edge_weight(j, k), **j))
    }

    fn choose_equiv_inv_edge_from(&self, k: &TngComplexKey) -> Option<&TngComplexKey> {
        self.inner.complex().vertex(k).out_edges().filter_map(|l|
            self.is_equiv_inv_edge(k, l).then_some(l)
        )
        .min_by_key(|l| (self.inner.complex().edge_weight(k, l), **l))
    }

    fn try_eliminate_equiv_at(&mut self, k: &TngComplexKey) -> bool {
        if let Some(&j) = self.choose_equiv_inv_edge_into(&k) { 
            self.eliminate_equiv(&j, &k);
            true
        } else if let Some(&l) = self.choose_equiv_inv_edge_from(&k) { 
            self.eliminate_equiv(&k, &l);
            true
        } else { 
            false
        }
    }

    fn eliminate_equiv(&mut self, i: &TngComplexKey, j: &TngComplexKey) {
        debug_assert_eq!(self.key_map.is_sym(i), self.key_map.is_sym(j));
        debug_assert!(self.inner.complex().has_edge(i, j));

        if self.key_map.is_sym(i) { 
            self.inner.eliminate(i, j);
        } else { 
            let ti = *self.key_map.inv_key(i);
            let tj = *self.key_map.inv_key(j);

            debug_assert!(self.inner.complex().has_edge(&ti, &tj));

            self.inner.eliminate(i, j);
            self.inner.eliminate(&ti, &tj);
        }

        self.key_map.remove(i);
        self.key_map.remove(j);
    }

    fn is_equiv_inv_edge(&self, i: &TngComplexKey, j: &TngComplexKey) -> bool { 
        let f = self.inner.complex().edge(i, j);
        f.is_invertible() && self.is_equiv_edge(i, j)
    }

    fn is_equiv_edge(&self, i: &TngComplexKey, j: &TngComplexKey) -> bool { 
        if self.key_map.is_sym(i) && self.key_map.is_sym(j) { 
            true
        } else if !self.key_map.is_sym(i) && !self.key_map.is_sym(j) { 
            //  i - - -> j 
            //    \   /   
            //      /     : not allowed
            //    /   \   
            // ti - - -> tj
            let ti = self.key_map.inv_key(i);
            let tj = self.key_map.inv_key(j);

            !self.inner.complex().has_edge(ti, j) &&
            !self.inner.complex().has_edge(i, tj)
        } else { 
            false
        }
    }

    fn finalize(&mut self) {
        info!("finalize: {}", self.inner.stat());

        if self.inner.complex().is_completely_delooped() { 
            return
        }

        self.deloop_all(false, false);
        self.deloop_all(true,  false); // deloop marked loops

        info!("  finalized: {}", self.inner.stat());
    }

    pub fn tau_map(&self) -> impl Fn(&KhGen) -> KhGen + Send + Sync + 'static {
        let key_map = self.key_map.clone();

        move |x: &KhGen| -> KhGen {
            let k = TngComplexKey::from(x);
            let tk = *key_map.inv_key(&k);
            tk.as_gen()
        }
    }

    pub fn into_inner(self) -> TngComplexBuilder<R> { 
        self.inner
    }

    pub fn into_tng_complex(self) -> TngComplex<R> { 
        self.inner.into_tng_complex()
    }

    fn inv_node(&self, x: &Node) -> &Node { 
        &self.x_map[x]
    }

    fn inv_edge(&self, e: Edge) -> Edge { 
        self.e_map[&e]
    }

    fn is_sym_comp(&self, c: &TngComp) -> bool {
        &c.convert_edges(|e| self.inv_edge(e)) == c
    }

    #[allow(unused)]
    fn print_keys(&self) {
        let mut done = HashSet::new();
        for k in self.key_map.keys().sorted() { 
            if done.contains(&k) { continue }

            let tk = self.key_map.inv_key(k);
            if k == tk {
                println!("{}", self.inner.complex().vertex(k));
            } else { 
                println!("{} ↔ {}", self.inner.complex().vertex(k), self.inner.complex().vertex(tk));
            }

            done.insert(k);
            done.insert(tk);
        }
        println!();
    }

}

/// Divide-and-conquer chunked build for a [`SymTngBuilder`]: select τ-closed chunks at thin
/// cutwidth boundaries, build each via a child builder, and merge the reduced chunk into the
/// parent — so the parent never materializes the full dense slice.
struct ChunkBuilder<'a, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    builder: &'a SymTngBuilder<R>,
}

impl<'a, R> ChunkBuilder<'a, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn new(builder: &'a SymTngBuilder<R>) -> Self {
        Self { builder }
    }

    fn is_on_axis(&self, x: &Node) -> bool {
        self.builder.inv_node(x) == x
    }

    // Next chunk: grow a τ-closed, boundary-connected piece (≤ `chunk_bound`), then cut it at the
    // first cutwidth valley past the peak — isolating the heavy region for a thin merge interface.
    fn next_chunk(&self) -> Option<Vec<Node>> {
        if self.builder.inner.nodes().next().is_none() { return None }
        let bound = self.builder.config.chunk_bound.unwrap_or(usize::MAX);
        let (chunk, cuts) = self.grow_chunk(bound);
        let cut = Self::cut_at_valley(&cuts).unwrap_or(chunk.len());
        Some(chunk[..cut].to_vec())
    }

    // Grow a τ-closed chunk ≤ `bound` with its proportional on-axis share (off-axis-heavy chunks
    // explode in preprocess). `open` is the connectivity frontier *and* the merge-cutwidth tracker.
    fn grow_chunk(&self, bound: usize) -> (Vec<Node>, Vec<(usize, usize)>) {
        let remaining = self.builder.inner.nodes().cloned().collect_vec();
        let n_on = remaining.iter().filter(|x| self.is_on_axis(x)).count();
        let target_on = if bound >= remaining.len() { n_on } else { (bound * n_on).div_ceil(remaining.len()) };

        let mut open: FxHashSet<Edge> = self.builder.inner.complex().boundary_ends().collect();
        let mut chunk: Vec<Node> = vec![];
        let mut cuts: Vec<(usize, usize)> = vec![];
        let mut on_taken = 0;

        while chunk.len() < bound {
            let Some(x) = self.pick_next(&remaining, &chunk, &open, on_taken < target_on) else { break };
            let tx = self.builder.inv_node(&x).clone();
            if tx == x { on_taken += 1; }
            for n in [&x, &tx] {
                if !chunk.contains(n) {
                    Self::toggle(&mut open, n.edges());
                    chunk.push(n.clone());
                }
            }
            cuts.push((chunk.len(), open.len()));
        }
        (chunk, cuts)
    }

    // The next crossing to absorb: boundary-connected and honoring the on-axis target, chosen by
    // node order (MinCut → smallest resulting cutwidth; LoopGreedy → first in crossing order).
    fn pick_next(&self, remaining: &[Node], chunk: &[Node], open: &FxHashSet<Edge>, prefer_on: bool) -> Option<Node> {
        let connected = |x: &&Node| !chunk.contains(*x) && x.edges().iter().any(|e| open.contains(e));
        let pick = |pool: Vec<&Node>| match self.builder.config.node {
            NodeOrder::MinCut => pool.into_iter().min_by_key(|x| self.unit_cutwidth(x, open)).cloned(),
            NodeOrder::LoopGreedy => pool.into_iter().next().cloned(),
        };
        pick(remaining.iter().filter(|x| connected(x) && self.is_on_axis(x) == prefer_on).collect())
            .or_else(|| pick(remaining.iter().filter(connected).collect()))
            .or_else(|| remaining.iter().find(|x| !chunk.contains(x)).cloned()) // seed when disconnected
    }

    // Merge cutwidth after absorbing the τ-unit of `x` (x and τx) into the simulated `open` set.
    fn unit_cutwidth(&self, x: &Node, open: &FxHashSet<Edge>) -> usize {
        let tx = self.builder.inv_node(x);
        let mut s = open.clone();
        Self::toggle(&mut s, x.edges());
        if tx != x { Self::toggle(&mut s, tx.edges()); }
        s.len()
    }

    // Build `chunk` into a reduced sub-complex via a child builder sharing the parent's
    // τ-maps; return the complex and its τ key_map for the parent merge.
    fn build_chunk(&self, chunk: &[Node]) -> (TngComplex<R>, TauKeyMap) {
        let child = self.child_builder(chunk).run();
        let SymTngBuilder { key_map, inner, .. } = child;
        (inner.into_tng_complex(), key_map)
    }

    // A child builder over `chunk` (a sub-tangle), inheriting the parent's τ-maps and
    // `preprocess`/simplify settings; `chunk_bound = None` so it doesn't recurse.
    fn child_builder(&self, chunk: &[Node]) -> SymTngBuilder<R> {
        let (h, t) = self.builder.inner.complex().ht();
        let base_pt = self.builder.inner.complex().base_pt();
        let mut inner = TngComplexBuilder::init(h, t, (0, 0), base_pt)
            .with_config(BuildConfig { mode: BuildMode::None, node: self.builder.config.node, h_range: None });
        inner.set_nodes(chunk.iter().cloned());

        // cap the child to the chunk's reachable band: a chunk vertex of weight
        // > b - deg_shift.0 can never reach the window (weight only grows).
        let h_range = self.builder.config.h_range.as_ref().map(|r| {
            let s = self.builder.inner.complex().deg_shift().0;
            0 ..= (*r.end() - s).max(0)
        });
        let config = SymBuildConfig { chunk_bound: None, h_range, ..self.builder.config.clone() };
        let key_map = TauKeyMap::init();
        let real_top = inner.complex().deg_shift().0 + chunk.len() as isize; // child deg_shift = 0
        SymTngBuilder { inner, x_map: self.builder.x_map.clone(), e_map: self.builder.e_map.clone(), key_map, config, real_top }
    }

    // Flip `edges` in the open-boundary set: each edge appears twice across crossings, so toggling
    // closes an already-open edge and opens a fresh one.
    fn toggle(open: &mut FxHashSet<Edge>, edges: &[Edge]) {
        for &e in edges {
            if !open.insert(e) { open.remove(&e); }
        }
    }

    // First cutwidth valley past the peak in a `(chunk_len, width)` profile → the chunk length to cut
    // at. None for a monotone-decreasing profile (a closing / last chunk) — caller takes it whole.
    fn cut_at_valley(cuts: &[(usize, usize)]) -> Option<usize> {
        let peak = cuts.iter().enumerate().max_by_key(|(_, (_, w))| *w)?.0;
        cuts[peak..].iter().tuple_windows()
            .find(|((_, w0), (_, w1))| w0 <= w1)
            .map(|((len, _), _)| *len)
    }
}

/// Builds the off-axis part of a [`SymTngBuilder`] by τ-symmetry: build one
/// representative half, mirror via τ, and merge both into the builder.
struct SymTngPreprocessor<'a, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    builder: &'a mut SymTngBuilder<R>,
}

impl<'a, R> SymTngPreprocessor<'a, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn run(builder: &'a mut SymTngBuilder<R>) {
        Self { builder }.build();
    }

    fn build(&mut self) {
        assert_eq!(self.builder.inner.complex().dim(), 0, "must start from init state.");

        let elements = self.builder.inner.take_elements();
        let (half, t_half) = self.partition_off_axis();

        info!("({}) preprocess off-axis: {} + {}", self.builder.inner.stat(), half.len(), t_half.len());

        let (c, tc, key_map, elements) = self.build_from_half(&half, elements);

        // merge the half, then its τ-image.
        self.merge_half(&half, c);
        self.merge_half(&t_half, tc);

        self.builder.key_map = key_map;
        self.builder.inner.set_elements(elements);

        info!("({}) preprocess done.", self.builder.inner.stat());
    }

    fn merge_half(&mut self, nodes: &[Node], c: TngComplex<R>) {
        info!("merge half {} into {}", c.stat(), self.builder.inner.stat());
        self.builder.inner.drop_nodes(|x| nodes.contains(x));
        self.builder.inner.merge(c);
        info!("  merged: {}", self.builder.inner.stat());
    }

    // Split the off-axis crossings (`τx != x`) into two τ-mirror halves: each adjacency
    // group goes opposite the side already holding its τ-image.
    fn partition_off_axis(&self) -> (Vec<Node>, Vec<Node>) {
        let off_axis = self.builder.inner.nodes().filter(|&x| self.builder.inv_node(x) != x).collect_vec();
        let groups = self.group_by_adjacency(&off_axis);

        let (mut half, mut t_half): (Vec<&Node>, Vec<&Node>) = (vec![], vec![]);
        for group in groups {
            let Some(&rep) = group.first() else { continue };
            if half.contains(&self.builder.inv_node(rep)) {
                t_half.extend(group);
            } else {
                half.extend(group);
            }
        }
        (half.into_iter().cloned().collect(), t_half.into_iter().cloned().collect())
    }

    // Union-find grouping: two nodes are adjacent iff they share a non-axis edge.
    fn group_by_adjacency<'b>(&self, nodes: &[&'b Node]) -> Vec<Vec<&'b Node>> {
        let shares_edge = |x: &Node, y: &Node|
            x.edges().iter()
                .filter(|&&e| self.builder.inv_edge(e) != e)
                .any(|e| y.edges().contains(e));

        let mut uf = KeyedUnionFind::from_iter(nodes.iter().copied());
        for (i, &x) in nodes.iter().enumerate() {
            for &y in &nodes[..i] {
                if shares_edge(x, y) {
                    uf.union(&x, &y);
                }
            }
        }
        uf.into_disjoint()
    }

    fn build_from_half(&self, crossings: &[Node], elements: Vec<TngComplexElem<R>>) -> (TngComplex<R>, TngComplex<R>, TauKeyMap, Vec<TngComplexElem<R>>) {
        // crossings appended after this chunk: all remaining nodes minus the chunk (half + τ-half).
        let r_rest = self.builder.inner.nodes().count() - 2 * crossings.len();

        let mut b = self.half_builder();
        b.set_nodes(crossings.iter().cloned());
        b.set_elements(elements);
        b.process_nodes();
        info!("half complex built: {}", b.stat());

        let keys = b.complex().keys().cloned().collect_vec();
        let mut elements = b.take_elements();

        let c = b.into_tng_complex();
        info!("mirror half via τ...");
        let tc = c.convert_edges(|e| self.builder.inv_edge(e));

        info!("pair key_map ({}² entries)...", keys.len());
        let key_map = TauKeyMap::from_half(&keys, self.builder.weight_band(r_rest));

        info!("complete {} elements...", elements.len());
        elements.iter_mut().for_each(|e|
            self.complete_element(e)
        );

        info!("build_from_half done: {} key-pairs.", key_map.len());
        (c, tc, key_map, elements)
    }

    // Inner builder for the half-complex; caps to `b - deg_shift.0`, since a
    // half-vertex of larger weight can never land in the window.
    fn half_builder(&self) -> TngComplexBuilder<R> {
        let (h, t) = self.builder.inner.complex().ht();
        let base_pt = self.builder.inner.complex().base_pt();
        let node = self.builder.config.node;
        let h_range = self.builder.config.h_range.as_ref().map(|w|
            0 ..= (*w.end() - self.builder.inner.complex().deg_shift().0)
        );
        TngComplexBuilder::init(h, t, (0, 0), base_pt)
            .with_config(BuildConfig { node, h_range, ..Default::default() })
    }


    // Complete a half-element into the full off-axis element.
    fn complete_element(&self, e: &mut TngComplexElem<R>) {
        let out = std::mem::take(e.out_cob_mut());
        *e.out_cob_mut() = out.into_iter().map(|(k, cob)| {
            let kk = k + k;
            let cc = cob.map(|c, r| {
                let tc = c.convert_edges(|e| self.builder.inv_edge(e));
                (c.connect(&tc), &r * &r)
            });
            (kk, cc)
        }).collect();
    }
}

#[cfg(test)]
#[allow(unused)]
mod tests {
    use crate::khi::{KhIGen, KhIHomology};

    use super::*;
    use num_traits::Zero;

    use yui_core::IteratorExt;
    use yui_core::lc::Lc;
    use yui_core::num::FF2;
    use yui_core::poly::Poly;
    use yui_core::RangeExt;
    use yui_homology::{ChainComplex1, ChainMap, ToSeqString, ToTableString};

    fn make_cone(b: SymTngBuilder<FF2>) -> ChainComplex1<KhIGen, FF2> { 
        let t = b.tau_map();
        let c = b.into_inner().into_tng_complex().into_raw_complex();
        let h_range = c.support().cloned().range().unwrap().mv(0, 1);
        let one_plus_tau = ChainMap::new(&c, &c, 0, move |_, z| {
            z.clone() + z.apply(|x| Lc::from(t(x)))
        });
        one_plus_tau.cone(h_range, false)
    }

    #[test]
    fn build_modes_agree() {
        // all build modes on the sym builder must agree on homology.
        let l = InvLink::from_symmetric_pd_code(
            [[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]
        );
        let (h, t) = (FF2::zero(), FF2::zero());
        let build = |mode| {
            let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
            b.config.mode = mode;
            b.run().into_tng_complex().into_raw_complex()
        };

        let ref_c = build(BuildMode::Greedy);
        let range = ref_c.support().cloned().range().unwrap();
        let ref_h = ref_c.homology();
        for mode in [BuildMode::Selective, BuildMode::MinFill, BuildMode::None] {
            let c = build(mode);
            c.check_d_all();
            let h = c.homology();
            for i in range.clone() {
                assert_eq!(h[i].rank(), ref_h[i].rank(), "rank at {i}, {mode:?}");
            }
        }
    }

    #[test]
    fn preprocess_matches() {
        let l = InvLink::test_data("6_3");
        let (h, t) = (FF2::zero(), FF2::zero());

        let build = |preprocess: bool, window: Option<RangeInclusive<isize>>| {
            let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
            b.config.preprocess = preprocess;
            b.config.h_range = window;
            b.run().into_tng_complex().into_raw_complex()
        };

        // full build: preprocess on/off agree on every degree.
        let (c_on, c_off) = (build(true, None), build(false, None));
        let range = c_off.support().cloned().range().unwrap();
        let (h_on, h_off) = (c_on.homology(), c_off.homology());
        for i in range.clone() {
            assert_eq!(h_on[i].rank(), h_off[i].rank(), "full rank at {i}");
        }

        // windowed build: agree on the interior of the window.
        let (a, b) = (*range.start() + 1, *range.end() - 1);
        let (c_on, c_off) = (build(true, Some(a..=b)), build(false, Some(a..=b)));
        let (h_on, h_off) = (c_on.homology(), c_off.homology());
        for i in (a + 1)..=(b - 1) {
            assert_eq!(h_on[i].rank(), h_off[i].rank(), "windowed rank at {i}");
        }
    }

    #[test]
    fn chunk_build_matches() {
        // k9_46 has enough crossings for several chunks at small bounds.
        let l = InvLink::from_symmetric_pd_code(
            [[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]
        );
        let (h, t) = (FF2::zero(), FF2::zero());

        let build = |chunk_bound: Option<usize>, window: Option<RangeInclusive<isize>>| {
            let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
            b.config.chunk_bound = chunk_bound;
            b.config.h_range = window;
            b.run().into_tng_complex().into_raw_complex()
        };

        let normal = build(None, None);
        let range = normal.support().cloned().range().unwrap();
        let hn = normal.homology();

        // full chunked builds match on every degree.
        for bound in [Some(2), Some(4), Some(100)] {
            let hc = build(bound, None).homology();
            for i in range.clone() {
                assert_eq!(hc[i].rank(), hn[i].rank(), "rank at {i} (bound {bound:?})");
            }
        }

        // windowed chunked build matches on the interior (exercises the band-cap).
        let (a, b) = (*range.start() + 1, *range.end() - 1);
        let hw = build(Some(4), Some(a..=b)).homology();
        for i in (a + 1)..=(b - 1) {
            assert_eq!(hw[i].rank(), hn[i].rank(), "windowed rank at {i}");
        }
    }

    #[test]
    fn test_kh_3_1() {
        let l = InvLink::test_data("3_1");
        let (h, t) = (FF2::zero(), FF2::zero());

        let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
        b.config.mode = BuildMode::None;
        let c = b.run().into_tng_complex().into_raw_complex();
        c.check_d_all();

        let h = c.homology();

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 2);
    }

    #[test]
    fn test_kh_3_1_h_range_full() {
        // A range covering the whole complex must reproduce the full homology.
        let l = InvLink::test_data("3_1");
        let (h, t) = (FF2::zero(), FF2::zero());

        let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
        b.config.mode = BuildMode::None;
        b.config.h_range = Some(0..=3);
        let c = b.run().into_tng_complex().into_raw_complex();
        c.check_d_all();

        let h = c.homology();

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 2);
    }

    #[test]
    fn test_kh_6_3_h_range() {
        let l = InvLink::test_data("6_3");
        let (h, t) = (FF2::zero(), FF2::zero());

        // full build for reference.
        let full = SymTngBuilder::from_inv_link(&l, &h, &t, false).run()
            .into_tng_complex().into_raw_complex();
        let range = full.support().cloned().range().unwrap();
        let (lo, hi) = (*range.start(), *range.end());
        let full_h = full.homology();

        // restrict to a strict sub-window.
        let (a, b) = (lo + 1, hi - 1);
        let mut bld = SymTngBuilder::from_inv_link(&l, &h, &t, false);
        bld.config.h_range = Some(a..=b);
        let bld = bld.run();

        // key_map stays in sync with the truncated complex.
        let verts: HashSet<_> = bld.inner.complex().keys().copied().collect();
        let kmap: HashSet<_> = bld.key_map.keys().copied().collect();
        assert_eq!(verts, kmap);

        let c = bld.into_tng_complex().into_raw_complex();
        c.check_d_all();
        let trunc_h = c.homology();

        // chain groups vanish outside [a, b].
        assert_eq!(c[a - 1].rank(), 0);
        assert_eq!(c[b + 1].rank(), 0);

        // strict interior degrees of the window are correct.
        for i in (a + 1)..=(b - 1) {
            assert_eq!(trunc_h[i].rank(), full_h[i].rank(), "rank mismatch at h-deg {i}");
        }
    }

    #[test]
    fn test_khi_3_1() {
        let l = InvLink::test_data("3_1");
        let (h, t) = (FF2::zero(), FF2::zero());

        let b = SymTngBuilder::from_inv_link(&l, &h, &t, false).run();
        let c = make_cone(b);
        c.check_d_all();

        let h = c.homology();

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 4);
        assert_eq!(h[4].rank(), 2);
    }

    #[test]
    fn no_auto_deloop() {
        // BuildMode::None defers all deloop to finalize and never eliminates, so the
        // finalized complex is delooped but unreduced.
        let l = InvLink::test_data("3_1");
        let (h, t) = (FF2::zero(), FF2::zero());

        let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
        b.config.mode = BuildMode::None;
        b.process_nodes();

        assert!(!b.inner.complex().is_completely_delooped());

        b.finalize();

        assert!(b.inner.complex().is_completely_delooped());

        let c = b.into_tng_complex().into_raw_complex();
        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 12);
        assert_eq!(c[3].rank(), 8);

        let h = c.homology();
        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 2);
    }
}