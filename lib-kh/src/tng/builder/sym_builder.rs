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
use delegate::delegate;
use rustc_hash::{FxHashMap, FxHashSet};
use itertools::Itertools;
use log::{debug, info};
use yui_core::algo::KeyedUnionFind;
use yui_core::bitseq::{Bit, BitSeq};
use yui_core::{Ring, RingOps};
use yui_link::{Node, Edge, InvLink};

use crate::kh::{KhGen, KhTensor};
use crate::tng::{LcCobTrait, TngComp, TngComplex, TngComplexElem, TngComplexKey};
use crate::tng::builder::{TngComplexBuilder, BuildConfig, BuildMode, NodeOrder};
use std::fmt;
use super::{reachable_range, pop_min_pivot, sparkline, cutwidth_after, toggle_boundary, boundary_edges};

/// Toggles for the automatic simplification done while building (kept separate
/// from [`BuildConfig`] so the equivariant builder can gain its own flags).
#[derive(Clone, Debug)]
pub struct SymBuildConfig {
    // crossing order: MinCut (default; bounds cutwidth for wide knots) or Given (PD order, debug).
    pub node_order: NodeOrder,
    pub mode: BuildMode,
    // build half the off-axis crossings and mirror via τ (see `preprocess`).
    pub preprocess: bool,
    // divide-and-conquer: partition the link into this many chunks at the deepest cutwidth
    // valleys, build each via a child builder, and merge into the parent. None = single pass.
    pub chunks: Option<usize>,
    // literal truncation: homology at the endpoints is wrong (build `(a-1)..=(b+1)` for correct `[a, b]`).
    pub h_range: Option<RangeInclusive<isize>>,
}

impl Default for SymBuildConfig {
    fn default() -> Self {
        Self { node_order: NodeOrder::default(), mode: BuildMode::default(), preprocess: true, chunks: None, h_range: None }
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
            .with_config(BuildConfig { mode: BuildMode::None, node_order: NodeOrder::default(), h_range: None });

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
        let inner_config = BuildConfig { mode: BuildMode::None, node_order: config.node_order, h_range: config.h_range.clone() };
        self.inner = self.inner.with_config(inner_config);
        self.config = config;
        self
    }

    delegate! {
        to self.inner {
            pub fn complex(&self) -> &TngComplex<R>;
            fn complex_mut(&mut self) -> &mut TngComplex<R>;
            pub fn set_elements<I>(&mut self, elements: I) where I: IntoIterator<Item = TngComplexElem<R>>;
            pub fn nodes(&self) -> &[Node];
            pub fn n_nodes(&self) -> usize;
            fn drop_nodes<F>(&mut self, pred: F) where F: Fn(&Node) -> bool;
            fn prepare_append(&mut self, x: &Node);
            fn cutwidth_of(&self, edges: impl IntoIterator<Item = Edge>) -> isize;
            fn find_loop_in(&self, k: &TngComplexKey, allow_based: bool) -> Option<usize>;
            fn deloop(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey>;
            fn eliminate(&mut self, i: &TngComplexKey, j: &TngComplexKey);
            fn collect_keys<F, W>(&self, i: isize, pred: F, weight: W) -> Vec<(TngComplexKey, usize)>
                where F: Fn(&TngComplexKey) -> bool, W: Fn(&TngComplexKey) -> usize;
            fn current_step(&self) -> String;
            pub(crate) fn stat(&self) -> String;
        }
    }

    pub fn run(mut self) -> Self {
        info!("build config:\n{:#?}", self.config);
        info!("cutwidth profile:\n{}", self.profile_sym());
        if self.config.chunks.is_some() {
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

    fn process_chunks(&mut self) {
        ChunkBuilder::run(self);
    }

    fn preprocess(&mut self) {
        SymTngPreprocessor::run(self);
    }

    fn process_nodes(&mut self) {
        info!("{} process {} nodes", self.current_step(), self.n_nodes());

        while let Some(x) = self.choose_next_node().cloned() {
            let tx = self.inv_node(&x).clone();
            if x == tx {
                self.append_on_axis(&x);
            } else {
                self.append_off_axis(&x, &tx);
            }
        }
    }

    /// Pick the next τ-unit by node order, ties broken by earliest crossing order. MinCut scores the
    /// *combined* x+τx toggle (a shared axis edge cancels — can't be summed).
    fn choose_next_node(&self) -> Option<&Node> {
        self.nodes().iter().enumerate()
            .min_by_key(|(i, x)| {
                let tx = self.inv_node(x);
                let score = match self.config.node_order {
                    NodeOrder::MinCut => {
                        let mut edges = x.edges().to_vec();
                        if tx != *x { edges.extend_from_slice(tx.edges()); }
                        -self.cutwidth_of(edges)
                    },
                    NodeOrder::Given => 0, // constant → ties broken by earliest index = given order
                };
                (-score, *i)
            })
            .map(|(_, x)| x)
    }

    fn append_on_axis(&mut self, x: &Node) { 
        info!("{} append on-axis: {x}", self.current_step());

        self.prepare_append(x);

        let (h, t) = self.complex().ht();
        let c = TngComplex::from_node(h, t, x, self.complex().base_pt());
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
        info!("{} append off-axis: {x}, {tx}", self.current_step());

        self.prepare_append(x);
        self.prepare_append(tx);

        let c = {
            let (h, t) = self.complex().ht();
            let mut c = TngComplex::from_node(h, t, x, self.complex().base_pt());
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
        // build the merged τ key-map per degree (next to merge_vertices) rather than as one
        // up-front cartesian — for large knots that product never fits in memory.
        let left_map = std::mem::take(&mut self.key_map);
        let (left, right) = self.complex_mut().prepare_merge(c);
        let range = reachable_range(self.complex().h_range(), &self.config.h_range, self.n_nodes());
        
        debug!("{} merge {} <- {}", self.current_step(), left.stat(), right.stat());
        debug!("  key_map: {} × {}", left_map.len(), right_map.len());
        debug!("  merge range: {:?}", range);

        match self.config.mode {
            BuildMode::None => for i in range { self.merge_slice(&left, &right, i, &left_map, &right_map); },
            _               => self.merge_incremental(&left, &right, range, &left_map, &right_map),
        }

        self.prune_h_range();

        debug!("{} merged: {}", self.current_step(), self.stat());
    }

    // Per degree: deloop, then (if the mode eliminates) sweep i-2,i-1 by equivariant Markowitz cost.
    // Greedy also inline-eliminates during deloop; the sweep just catches what it missed.
    fn merge_incremental(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, range: RangeInclusive<isize>, left_map: &TauKeyMap, right_map: &TauKeyMap) {
        debug_assert!(self.config.mode.auto_deloop()); // None is dispatched to merge_slice
        let top = *range.end();

        for i in range {
            debug!("{} build C[{i}]...", self.current_step());
            self.merge_slice(left, right, i, left_map, right_map);
            self.deloop_in(i - 1);
            if self.config.mode.auto_elim() {
                self.eliminate_in(i - 2);
                self.eliminate_in(i - 1);
            }
            debug!("{} built C[{i}]: {}", self.current_step(), self.complex().rank(i));
        }

        self.prune_isolated_top(top);
        self.deloop_in(top);
        if self.config.mode.auto_elim() {
            self.eliminate_in(top - 1);
        }
    }

    // Build degree `i`: the τ key-map slice, then its vertices and the edges into it.
    fn merge_slice(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, i: isize, left_map: &TauKeyMap, right_map: &TauKeyMap) {
        for (k1, k2) in TngComplex::collect_keys(left, right, i) {
            self.key_map.add_pair(k1 + k2, left_map.inv_key(k1) + right_map.inv_key(k2));
        }
        self.inner.merge_slice(left, right, i);
    }

    // No-in-edge vertices at a TRUNCATED window-top (`top < real_top`) feed only the discarded
    // boundary — drop them and their τ-pairs (at the real top they're genuine generators).
    // The no-in-edge set is τ-closed, so `key_map` stays an involution.
    fn prune_isolated_top(&mut self, top: isize) {
        let truncated = self.config.h_range.as_ref().is_some_and(|w| top == *w.end()) && top < self.real_top;
        if !truncated { return; }

        let doomed_verts = self.complex().keys_of_deg(top)
            .filter(|k| self.complex().vertex(k).in_edges().next().is_none())
            .copied()
            .collect_vec();
        if !doomed_verts.is_empty() {
            debug!("prune {} isolated verts in C[{top}].", doomed_verts.len());
        }
        let doomed: FxHashSet<_> = doomed_verts.iter().copied().collect();
        self.complex_mut().remove_vertices(&doomed_verts);
        self.key_map.drop(|k| doomed.contains(k));
    }

    /// Prune doomed vertices *and* their `key_map` entries. τ preserves weight,
    /// so `deg(k) == deg(τk)` — pruning is symmetric and `key_map` stays a valid
    /// involution. The `key_map` is pruned independently of the complex because
    /// the cartesian merge over-generates it past the vertex cap.
    fn prune_h_range(&mut self) {
        let Some(h_range) = self.config.h_range.clone() else { return };
        let r = self.n_nodes() as isize;
        let i0 = self.complex().deg_shift().0;

        // a vertex of degree `d` reaches `[d, d + r]`, so it stays relevant iff
        // `d ∈ [a - r, b]` — current degrees that can still land in `h_range`.
        let live = (*h_range.start() - r) ..= *h_range.end();
        let doomed = |k: &TngComplexKey|
            !live.contains(&(k.weight() as isize + i0));

        let doomed_verts = self.complex().keys_of(&doomed).copied().collect_vec();
        if !doomed_verts.is_empty() {
            debug!("prune {} verts outside h_range.", doomed_verts.len());
        }
        self.complex_mut().remove_vertices(&doomed_verts);

        self.key_map.drop(doomed);
    }

    // Deloop unmarked loops over all degrees.
    fn deloop_all(&mut self) {
        for i in self.complex().h_range() {
            self.deloop_in(i);
        }
    }

    // Deloop the marked (based) loops into a single summand. All unmarked loops must already
    // be gone — else delooping only the marked ones would break the complex.
    fn deloop_all_marked(&mut self) {
        debug_assert!(
            self.complex().keys().all(|k| self.find_loop_in(k, false).is_none()),
            "deloop_all_marked: unmarked loops remain"
        );
        for i in self.complex().h_range() {
            self.deloop_in_with(i, true);
        }
    }

    fn deloop_in(&mut self, i: isize) {
        self.deloop_in_with(i, false);
    }

    // Pivot cost for selection: off-axis pivots are eliminated in τ-pairs, so count ~2× the fill.
    fn pivot_weight(&self, k: &TngComplexKey) -> usize {
        let w = self.complex().vertex(k).c_weight();
        if self.key_map.is_sym(k) { w } else { 2 * w }
    }

    // Markowitz cost of eliminating `k → l`; off-axis pivots eliminate in τ-pairs (~2× the fill).
    fn equiv_edge_weight(&self, k: &TngComplexKey, l: &TngComplexKey) -> usize {
        let w = self.complex().edge_weight(k, l);
        if self.key_map.is_sym(k) { w } else { 2 * w }
    }

    // Least Markowitz cost to eliminate `k`, over its equiv-invertible incident edges — the
    // equivariant pivot priority (vs the cruder `pivot_weight`).
    fn equiv_elim_cost(&self, k: &TngComplexKey) -> usize {
        let v = self.complex().vertex(k);
        let outs = v.out_edges().filter(|l| self.is_equiv_inv_edge(k, l)).map(|l| self.equiv_edge_weight(k, l));
        let ins = v.in_edges().filter(|j| self.is_equiv_inv_edge(j, k)).map(|j| self.equiv_edge_weight(j, k));
        outs.chain(ins).min().unwrap_or(0)
    }

    fn deloop_in_with(&mut self, i: isize, allow_based: bool) {
        let mut keys = self.collect_keys(i,
            |k| self.find_loop_in(k, allow_based).is_some(),
            |k| self.pivot_weight(k),
        );
        if keys.is_empty() { return }

        debug!("{} deloop in C[{i}], targets: {}.", self.current_step(), keys.len());

        let before = self.complex().rank(i) as isize;

        while let Some(k) = pop_min_pivot(&mut keys, |k|
            self.complex().contains_key(k).then(|| self.pivot_weight(k))
        ) {
            let Some(r) = self.find_loop_in(&k, allow_based) else { continue };

            for new_key in self.deloop_equiv(&k, r) {
                if self.find_loop_in(&new_key, allow_based).is_some() {
                    let w = self.pivot_weight(&new_key);
                    keys.push((new_key, w));
                }
            }
        }

        let after = self.complex().rank(i) as isize;

        debug!("{}   delooped C[{i}]: {} (diff: {}).", self.current_step(), after, after - before);
    }

    fn deloop_equiv(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> { 
        let mut added = if self.key_map.is_sym(k) {
            let c = self.complex().vertex(k).tng().comp(r);
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
                self.complex().contains_key(k) &&
                !self.try_eliminate_equiv_at(k)
            );
            // an equivariant elim removes the whole τ-pair, so a key kept above may since have
            // been eliminated as another's τ-partner — re-retain so only live keys are returned.
            added.retain(|k| self.complex().contains_key(k));
        }
        added
    }

    fn deloop_on_axis_sym(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> {
        let c = self.complex().vertex(k).tng().comp(r);

        debug_assert!(self.key_map.is_sym(k));
        debug_assert!(self.is_sym_comp(c));

        let updated = self.deloop(k, r);

        self.key_map.remove(k);

        for &k_new in updated.iter() { 
            self.key_map.add_pair(k_new, k_new);
        }

        updated
    }

    #[allow(non_snake_case)]
    fn deloop_on_axis_asym(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> {
        let c = self.complex().vertex(k).tng().comp(r);

        debug_assert!(self.key_map.is_sym(k));
        debug_assert!(!self.is_sym_comp(c));
        debug_assert!(!c.is_marked());

        //          ⚪︎1 | ⚪︎1
        //  ⚪︎1 | ⚪︎X  <-->  ⚪︎X | ⚪︎1
        //          ⚪︎X | ⚪︎X

        let tc = c.convert_edges(|e| self.inv_edge(e));

        let ks = self.deloop(k, r);

        let (k_X, k_1) = (ks[0], ks[1]);
        let (k_XX, k_X1) = { 
            let tr = self.complex().vertex(&k_X).tng().index_of(&tc).unwrap();
            let tks = self.deloop(&k_X, tr);
            (tks[0], tks[1])
        };
        let (k_1X, k_11) = { 
            let tr = self.complex().vertex(&k_1).tng().index_of(&tc).unwrap();
            let tks = self.deloop(&k_1, tr);
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
        let c = self.complex().vertex(k).tng().comp(r);

        debug_assert!(!self.key_map.is_sym(k));

        //  ⚪︎1 | ..  <-->  .. | ⚪︎1
        //  ⚪︎X | ..  <-->  .. | ⚪︎X

        let tk = *self.key_map.inv_key(k);
        let tc = c.convert_edges(|e| self.inv_edge(e));
        let tr = self.complex().vertex(&tk).tng().index_of(&tc).unwrap();

        let mut ks = self.deloop(k, r);
        let mut tks = self.deloop(&tk, tr);

        self.key_map.remove(k);

        for (&k_new, &tk_new) in Iterator::zip(ks.iter(), tks.iter()) { 
            self.key_map.add_pair(k_new, tk_new);
        }

        ks.append(&mut tks);
        ks
    }

    fn eliminate_in(&mut self, i: isize) {
        let mut keys = self.collect_keys(i,
            |k| self.complex().vertex(k).out_edges()
                .any(|l| self.is_equiv_inv_edge(k, l)),
            |k| self.equiv_elim_cost(k),
        );
        if keys.is_empty() { return }

        debug!("{} eliminate in C[{i}], targets: {}", self.current_step(), keys.len());

        let before = self.complex().rank(i) as isize;

        while let Some(k) = pop_min_pivot(&mut keys, |k|
            self.complex().contains_key(k).then(|| self.equiv_elim_cost(k))
        ) {
            self.try_eliminate_equiv_at(&k);
        }

        let after = self.complex().rank(i) as isize;

        debug!("{}   eliminated C[{i}]: {} (diff: {}).", self.current_step(), after, after - before);
    }

    fn choose_equiv_inv_edge_into(&self, k: &TngComplexKey) -> Option<&TngComplexKey> { 
        self.complex().vertex(k).in_edges().filter_map(|j|
            self.is_equiv_inv_edge(j, k).then_some(j)
        )
        .min_by_key(|j| (self.complex().edge_weight(j, k), **j))
    }

    fn choose_equiv_inv_edge_from(&self, k: &TngComplexKey) -> Option<&TngComplexKey> {
        self.complex().vertex(k).out_edges().filter_map(|l|
            self.is_equiv_inv_edge(k, l).then_some(l)
        )
        .min_by_key(|l| (self.complex().edge_weight(k, l), **l))
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
        debug_assert!(self.complex().has_edge(i, j));

        if self.key_map.is_sym(i) { 
            self.eliminate(i, j);
        } else { 
            let ti = *self.key_map.inv_key(i);
            let tj = *self.key_map.inv_key(j);

            debug_assert!(self.complex().has_edge(&ti, &tj));

            self.eliminate(i, j);
            self.eliminate(&ti, &tj);
        }

        self.key_map.remove(i);
        self.key_map.remove(j);
    }

    fn is_equiv_inv_edge(&self, i: &TngComplexKey, j: &TngComplexKey) -> bool { 
        let f = self.complex().edge(i, j);
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

            !self.complex().has_edge(ti, j) &&
            !self.complex().has_edge(i, tj)
        } else { 
            false
        }
    }

    fn finalize(&mut self) {
        if self.complex().is_completely_delooped() { 
            info!("{} completely delooped: {}", self.current_step(), self.stat());
            return
        }

        info!("{} finalize: {}", self.current_step(), self.stat());

        self.deloop_all();
        self.deloop_all_marked(); // deloop marked loops

        info!("{} finalized: {}", self.current_step(), self.stat());
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
                println!("{}", self.complex().vertex(k));
            } else { 
                println!("{} ↔ {}", self.complex().vertex(k), self.complex().vertex(tk));
            }

            done.insert(k);
            done.insert(tk);
        }
        println!();
    }

    /// Boundary-cutwidth profile of the symmetric (τ-equivariant) MinCut order: on-axis crossings
    /// singly, off-axis in `(x, τx)` pairs scored by combined cutwidth (shared axis edges cancel).
    pub(crate) fn profile_sym(&self) -> SymBuildProfile {
        let nodes = self.nodes();
        let tau = &self.x_map;
        let idx_of: FxHashMap<Node, usize> = nodes.iter().enumerate()
            .map(|(i, x)| (x.clone(), i))
            .collect();

        // each pair is kept once, at its lower index; on-axis crossings stay singletons
        let units: Vec<Vec<usize>> = (0..nodes.len()).filter_map(|i| {
            let j = idx_of[&tau[&nodes[i]]];
            match j {
                _ if j == i => Some(vec![i]),
                _ if i < j  => Some(vec![i, j]),
                _           => None,
            }
        }).collect();

        let on_axis = units.iter().filter(|u| u.len() == 1).count();
        let off_axis = nodes.len() - on_axis;

        let unit_nodes = |u: &[usize]| -> Vec<&Node> { u.iter().map(|&i| &nodes[i]).collect() };

        let mut remaining: Vec<usize> = (0..units.len()).collect();
        let mut open: FxHashSet<Edge> = FxHashSet::default();

        let (order, widths): (Vec<Vec<usize>>, Vec<usize>) = std::iter::from_fn(|| {
            let pos = remaining.iter()
                .position_min_by_key(|&&u| (cutwidth_after(&open, &unit_nodes(&units[u])), units[u][0]))?;
            let u = remaining.swap_remove(pos);
            toggle_boundary(&mut open, &unit_nodes(&units[u]));
            Some((units[u].clone(), open.len()))
        }).unzip();

        let peak = widths.iter().copied().max().unwrap_or(0);
        SymBuildProfile { on_axis, off_axis, order, widths, peak }
    }

}

/// Boundary-cutwidth profile of the symmetric (τ-equivariant) MinCut order.
pub(crate) struct SymBuildProfile {
    pub on_axis: usize,
    pub off_axis: usize,        // node count; pairs = off_axis / 2
    #[allow(dead_code)] // replayed only by the faithfulness test
    pub order: Vec<Vec<usize>>, // each unit = 1 (on-axis) or 2 (τ-pair) node indices
    pub widths: Vec<usize>,
    pub peak: usize,
}

impl fmt::Display for SymBuildProfile {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "on-axis:  {}", self.on_axis)?;
        writeln!(f, "off-axis: {} ({} pairs)", self.off_axis, self.off_axis / 2)?;
        writeln!(f, "peak:     {}", self.peak)?;
        write!(f, "{}", sparkline(&self.widths, self.peak))
    }
}

/// Divide-and-conquer chunked build for a [`SymTngBuilder`]: plan k τ-closed chunks up front at
/// thin cutwidth interfaces, build each via a child builder, and merge the reduced chunk into the
/// parent — so the parent never materializes the full dense slice.
struct ChunkBuilder<'a, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    builder: &'a mut SymTngBuilder<R>,
}

impl<'a, R> ChunkBuilder<'a, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn run(builder: &'a mut SymTngBuilder<R>) {
        Self { builder }.build();
    }

    // Plan k chunks up front, then build each reduced chunk and merge it into the parent.
    fn build(&mut self) {
        // canon cycles aren't yet tracked through chunk merges (see `merge`'s TODO), so
        // drop them — a chunked build yields the complex/homology but not α / ssi.
        self.builder.set_elements(vec![]);
        let plan = self.plan();
        info!("{} chunk plan: {} pieces {:?}", self.builder.current_step(), plan.len(),
            plan.iter().map(|c| c.len()).collect_vec());
        for chunk in plan {
            let (c, key_map) = self.build_chunk(&chunk);
            self.builder.drop_nodes(|x| chunk.contains(x));
            self.builder.merge(c, key_map);
            info!("{} chunk merged: {}", self.builder.current_step(), self.builder.stat());
        }
    }

    // Partition the crossings into `chunks` pieces at the deepest cutwidth valleys of the MinCut
    // order. Each piece is τ-closed (τ-units stay whole) and contiguous in that order, so merging
    // them in sequence keeps a thin interface at every step.
    fn plan(&self) -> Vec<Vec<Node>> {
        let prof = self.builder.profile_sym();
        let k = self.builder.config.chunks.unwrap_or(1).max(1);
        let nodes = self.builder.inner.nodes();
        let cuts = Self::select_cuts(&prof.widths, k - 1);

        // segment `prof.order` after each cut position; expand each unit to its nodes.
        let starts = std::iter::once(0).chain(cuts.iter().map(|&v| v + 1));
        let ends = cuts.iter().map(|&v| v + 1).chain(std::iter::once(prof.order.len()));
        starts.zip(ends)
            .map(|(s, e)| prof.order[s..e].iter().flatten().map(|&i| nodes[i].clone()).collect())
            .filter(|c: &Vec<Node>| !c.is_empty())
            .collect()
    }

    // Valley positions (a descent that turns back up) in `widths`, ascending. The `scan` carries
    // a `descended` flag so flats are ignored — a mid-descent plateau isn't mistaken for the bottom.
    fn valleys(widths: &[usize]) -> Vec<usize> {
        widths.windows(2).enumerate()
            .scan(false, |descended, (i, w)| Some(
                if w[1] < w[0] {
                    *descended = true;
                    None
                } else if w[1] > w[0] && *descended {
                    *descended = false;
                    Some(i)
                } else {
                    None
                }
            ))
            .flatten()
            .collect()
    }

    // `n_cuts` cut positions in `widths`: take the deepest valleys first; if more cuts than
    // valleys are needed, keep every valley and split the widest pieces evenly (greedy, which
    // minimizes the largest piece).
    fn select_cuts(widths: &[usize], n_cuts: usize) -> Vec<usize> {
        let n = widths.len();
        let n_cuts = n_cuts.min(n.saturating_sub(1));
        if n_cuts == 0 {
            return vec![];
        }

        let vs = Self::valleys(widths); // ascending positions
        if vs.len() >= n_cuts {
            // enough valleys: keep the `n_cuts` deepest, back in position order.
            return vs.into_iter()
                .sorted_by_key(|&p| widths[p])
                .take(n_cuts)
                .sorted()
                .collect();
        }

        // too few valleys: every valley is a cut; spend the rest splitting the widest pieces
        // evenly. `alloc[i]` is the number of pieces segment `i` is divided into.
        let segs: Vec<(usize, usize)> = std::iter::once(0)
            .chain(vs.iter().map(|&v| v + 1))
            .chain(std::iter::once(n))
            .tuple_windows()
            .collect();

        let alloc = (vs.len()..n_cuts).fold(vec![1usize; segs.len()], |mut alloc, _| {
            let widest = (0..segs.len())
                .max_by(|&a, &b| ((segs[a].1 - segs[a].0) * alloc[b]).cmp(&((segs[b].1 - segs[b].0) * alloc[a])))
                .unwrap();
            alloc[widest] += 1;
            alloc
        });

        let even = segs.iter().zip(&alloc)
            .flat_map(|(&(lo, hi), &a)| (1..a).map(move |j| lo + j * (hi - lo) / a - 1));
        vs.into_iter().chain(even).sorted().collect()
    }

    // Build `chunk` into a reduced sub-complex via a child builder, returning it with its τ key-map.
    fn build_chunk(&self, chunk: &[Node]) -> (TngComplex<R>, TauKeyMap) {
        let step = self.builder.current_step();
        let ends = boundary_edges(&chunk.iter().collect::<Vec<_>>()).into_iter().sorted().collect_vec();
        info!("{step} build chunk (n: {}, nb: {} {:?}): {}", chunk.len(), ends.len(), ends, chunk.iter().join(", "));

        let child = self.child_builder(chunk).run();
        let SymTngBuilder { key_map, inner, .. } = child;
        let c = inner.into_tng_complex();
        info!("{step} chunk built: {}", c.stat());
        (c, key_map)
    }

    // A child builder over `chunk` (a sub-tangle), inheriting the parent's τ-maps and
    // `preprocess`/simplify settings; chunking always uses the MinCut order, never recursing.
    fn child_builder(&self, chunk: &[Node]) -> SymTngBuilder<R> {
        let (h, t) = self.builder.inner.complex().ht();
        let base_pt = self.builder.inner.complex().base_pt();
        let mut inner = TngComplexBuilder::init(h, t, (0, 0), base_pt)
            .with_config(BuildConfig { mode: BuildMode::None, node_order: NodeOrder::MinCut, h_range: None });
        inner.set_nodes(chunk.iter().cloned());

        // cap the child to the chunk's reachable band: a chunk vertex of weight
        // > b - deg_shift.0 can never reach the window (weight only grows).
        let h_range = self.builder.config.h_range.as_ref().map(|r| {
            let s = self.builder.inner.complex().deg_shift().0;
            0 ..= (*r.end() - s).max(0)
        });
        let config = SymBuildConfig { chunks: None, node_order: NodeOrder::MinCut, h_range, ..self.builder.config.clone() };
        let key_map = TauKeyMap::init();
        let real_top = inner.complex().deg_shift().0 + chunk.len() as isize; // child deg_shift = 0
        SymTngBuilder { inner, x_map: self.builder.x_map.clone(), e_map: self.builder.e_map.clone(), key_map, config, real_top }
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

        info!("{} preprocess off-axis: {} + {}", self.builder.current_step(), half.len(), t_half.len());

        let (c, tc, key_map, elements) = self.build_from_half(&half, elements);

        // merge the half, then its τ-image.
        self.merge_half(&half, c);
        self.merge_half(&t_half, tc);

        self.builder.key_map = key_map;
        self.builder.inner.set_elements(elements);

        info!("{} preprocess done: {}", self.builder.current_step(), self.builder.stat());
    }

    fn merge_half(&mut self, nodes: &[Node], c: TngComplex<R>) {
        self.builder.inner.drop_nodes(|x| nodes.contains(x));
        self.builder.inner.merge(c);
    }

    // Split the off-axis crossings (`τx != x`) into two τ-mirror halves: each adjacency
    // group goes opposite the side already holding its τ-image.
    fn partition_off_axis(&self) -> (Vec<Node>, Vec<Node>) {
        let off_axis = self.builder.inner.nodes().iter().filter(|&x| self.builder.inv_node(x) != x).collect_vec();
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
        let r_rest = self.builder.n_nodes() - 2 * crossings.len();

        let mut b = self.half_builder();
        b.set_nodes(crossings.iter().cloned());
        b.set_elements(elements);
        b.process_nodes();

        info!("{} half complex built: {}", self.builder.current_step(), b.stat());

        let keys = b.complex().keys().cloned().collect_vec();
        let mut elements = b.take_elements();

        let c = b.into_tng_complex();
        debug!("  mirror half via τ...");
        let tc = c.convert_edges(|e| self.builder.inv_edge(e));

        debug!("  pair key_map ({}² entries)...", keys.len());
        let key_map = TauKeyMap::from_half(&keys, self.weight_band(r_rest));

        debug!("  complete {} elements...", elements.len());
        elements.iter_mut().for_each(|e|
            self.complete_element(e)
        );

        (c, tc, key_map, elements)
    }

    // Combined off-axis weight that can still reach the window, given `r` pending crossings.
    fn weight_band(&self, r: usize) -> RangeInclusive<usize> {
        let s = self.builder.complex().deg_shift().0;
        let window = self.builder.config.h_range.as_ref().map(|w| (*w.start() - s) ..= (*w.end() - s));
        let band = reachable_range(0 ..= isize::MAX, &window, r);
        (*band.start()).max(0) as usize ..= (*band.end()).max(0) as usize
    }

    // Inner builder for the half-complex; caps to `b - deg_shift.0`, since a
    // half-vertex of larger weight can never land in the window.
    fn half_builder(&self) -> TngComplexBuilder<R> {
        let (h, t) = self.builder.inner.complex().ht();
        let base_pt = self.builder.inner.complex().base_pt();
        let node_order = self.builder.config.node_order;
        let h_range = self.builder.config.h_range.as_ref().map(|w|
            0 ..= (*w.end() - self.builder.inner.complex().deg_shift().0)
        );
        TngComplexBuilder::init(h, t, (0, 0), base_pt)
            .with_config(BuildConfig { node_order, h_range, ..Default::default() })
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

    // `profile_sym`'s dry-run widths must match the real complex's `boundary_ends`. Boundary depends
    // only on the processed *set*, not the pairing, so a plain builder is a valid oracle.
    #[test]
    fn sym_dry_run_matches_boundary() {
        for name in ["3_1", "6_3"] {
            let l = InvLink::test_data(name);
            let prof = SymTngBuilder::<i32>::from_inv_link(&l, &0, &0, false).profile_sym();
            let nodes: Vec<Node> = l.nodes().cloned().collect();

            assert_eq!(prof.on_axis + prof.off_axis, nodes.len(), "{name}: unit counts");

            let mut b = TngComplexBuilder::<i32>::init(&0, &0, (0, 0), None)
                .with_config(BuildConfig { mode: BuildMode::None, ..Default::default() });

            let mut open: FxHashSet<Edge> = FxHashSet::default();
            for (step, unit) in prof.order.iter().enumerate() {
                unit.iter().for_each(|&i| b.append_node(&nodes[i]));
                let u: Vec<&Node> = unit.iter().map(|&i| &nodes[i]).collect();
                toggle_boundary(&mut open, &u);
                let real: FxHashSet<Edge> = b.complex().boundary_ends().collect();
                assert_eq!(open, real, "{name} step {step}: open-set vs boundary_ends");
                assert_eq!(prof.widths[step], open.len(), "{name} step {step}: width");
            }
            assert_eq!(*prof.widths.last().unwrap(), 0, "{name} should close up");
        }
    }

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
        for mode in [BuildMode::MinFill, BuildMode::NoElim, BuildMode::None] {
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
        // k9_46 has enough crossings to split into several chunks.
        let l = InvLink::from_symmetric_pd_code(
            [[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]
        );
        let (h, t) = (FF2::zero(), FF2::zero());

        let build = |chunks: Option<usize>, window: Option<RangeInclusive<isize>>| {
            let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
            b.config.chunks = chunks;
            b.config.h_range = window;
            b.run().into_tng_complex().into_raw_complex()
        };

        let normal = build(None, None);
        let range = normal.support().cloned().range().unwrap();
        let hn = normal.homology();

        // chunked builds match on every degree.
        for k in [Some(2), Some(3), Some(4)] {
            let hc = build(k, None).homology();
            for i in range.clone() {
                assert_eq!(hc[i].rank(), hn[i].rank(), "rank at {i} (chunks {k:?})");
            }
        }

        // windowed chunked build matches on the interior (exercises the band-cap).
        let (a, b) = (*range.start() + 1, *range.end() - 1);
        let hw = build(Some(3), Some(a..=b)).homology();
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

    #[test]
    fn no_auto_elim() {
        // NoElim deloops every circle during the build but never eliminates → delooped, unreduced.
        let l = InvLink::test_data("3_1");
        let (h, t) = (FF2::zero(), FF2::zero());

        let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
        b.config.mode = BuildMode::NoElim;
        b.process_nodes();

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