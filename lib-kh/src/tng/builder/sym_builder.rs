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
use itertools::{iproduct, Itertools};
use log::{debug, info};
use yui_core::algo::KeyedUnionFind;
use yui_core::bitseq::{Bit, BitSeq};
use yui_core::abst::{Ring, RingOps};
use yui_link::{Node, Edge, InvLink};

use crate::kh::{KhAlgGen, KhGen, KhTensor};
use crate::tng::{ElimDir, LcCob, LcCobTrait, TngComp, TngComplex, TngComplexElem, TngComplexKey};
use crate::tng::builder::{TngComplexBuilder, TngElemBuilder, BuildConfig, Strategy, NodeOrder};
use std::fmt;
use super::{assert_supported_symmetry, reachable_range, pop_min_pivot, pivot_pool, push_pivot, sparkline, fill_cost_sparkline, cutwidth_after, toggle_boundary, BuildPlanner, CutOption};
use super::builder::PROGRESS_LOG_STEP;
use log::Level;
use yui_core::util::log::log_progress;

/// Toggles for the automatic simplification done while building (kept separate
/// from [`BuildConfig`] so the equivariant builder can gain its own flags).
#[derive(Clone, Debug)]
pub struct SymBuildConfig {
    // crossing order: MinCut (default; bounds cutwidth for wide knots) or Given (PD order, debug).
    pub node_order: NodeOrder,
    pub strategy: Strategy,
    // build half the off-axis crossings and mirror via τ (see `preprocess`).
    pub preprocess: bool,
    // literal truncation: homology at the endpoints is wrong (build `(a-1)..=(b+1)` for correct `[a, b]`).
    pub h_range: Option<RangeInclusive<isize>>,
    // drop generators outside this q-range. Applied once the diagram is closed (exact q); τ is
    // q-homogeneous, so this stays consistent through the cone. `None` = no q-truncation.
    pub q_range: Option<RangeInclusive<isize>>,
    // divide-and-conquer chunking (auto cutwidth or manual edge-cuts); None = single pass.
    pub cut: CutOption,
    // cone only: cap the per-elimination fill cost during cone_merge; survivors defer to the
    // matrix reduction (two-pass: cheap cobordism elim, then scalar F2[H] reduction). None = no cap.
    pub max_elim_cost: Option<usize>,
    // skip the final deloop (the last merge and `finalize`): remaining circles defer to
    // `into_raw_complex`'s matrix-level expansion + the `ChainReducer`. See `should_deloop`.
    pub no_full_deloop: bool,
}

impl Default for SymBuildConfig {
    fn default() -> Self {
        Self { node_order: NodeOrder::default(), strategy: Strategy::default(), preprocess: true, h_range: None, q_range: None, cut: CutOption::None, max_elim_cost: None, no_full_deloop: false }
    }
}

impl SymBuildConfig {
    // Config for the inner merge builder: `strategy: None` (the sym builder drives deloop/elim, the inner
    // never orders nodes), so only `h_range` carries over — to drop out-of-window canon cycles.
    pub(crate) fn inner_build_config(&self) -> BuildConfig {
        BuildConfig { strategy: Strategy::None, h_range: self.h_range.clone(), q_range: self.q_range.clone(), ..Default::default() }
    }
}

// τ-symmetric key map: each `TngComplexKey`'s τ-image. On-axis keys (`τk = k`) are a set;
// off-axis keys form an involution stored both ways for O(1) `inv_key`. The merge iterates
// only one representative per off-axis pair (the map is symmetric), then symmetrizes.
#[derive(Clone, Default)]
pub(crate) struct TauKeyMap {
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
    pub(crate) fn init() -> Self {
        Self::from_iter([(TngComplexKey::init(), TngComplexKey::init())])
    }

    fn len(&self) -> usize {
        self.on_axis.len() + self.off_axis.len()
    }

    pub(crate) fn inv_key(&self, k: &TngComplexKey) -> &TngComplexKey {
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
    pub(crate) fn drop(&mut self, pred: impl Fn(&TngComplexKey) -> bool) {
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
        assert_supported_symmetry(l);
        assert!(
            !reduced || l.base_pt().is_some_and(|e| l.is_on_axis(e)),
            "reduced requires a base point on the axis"
        );

        // the inner builder is driven by `self` — disable its own auto-simplify.
        let config = SymBuildConfig::default();
        let inner = TngComplexBuilder::from_link(l.inner(), h, t, reduced)
            .with_config(config.inner_build_config());

        let x_map = l.nodes().map(|x|
            (x.clone(), l.inv_node(x).clone())
        ).collect();
        let e_map = l.edges().into_iter().map(|e| (e, l.inv_edge(e))).collect();
        let key_map = TauKeyMap::init();
        let real_top = inner.complex().deg_shift().0 + l.inner().n_crossings() as isize;

        SymTngBuilder { inner, x_map, e_map, key_map, config, real_top }
    }

    pub fn with_config(mut self, config: SymBuildConfig) -> Self {
        // propagate the window to the inner builder so the preprocess merges cap
        // to it; this also drops canon cycles when the window excludes h-degree 0.
        self.inner = self.inner.with_config(config.inner_build_config());
        self.config = config;
        self
    }

    pub fn config(&self) -> &SymBuildConfig {
        &self.config
    }

    delegate! {
        to self.inner {
            pub fn complex(&self) -> &TngComplex<R>;
            pub(crate) fn complex_mut(&mut self) -> &mut TngComplex<R>;
            pub(crate) fn elements(&self) -> &TngElemBuilder<R>;
            pub fn nodes(&self) -> &[Node];
            pub fn n_nodes(&self) -> usize;
            pub(crate) fn drop_nodes<F>(&mut self, pred: F) where F: Fn(&Node) -> bool;
            fn prepare_append(&mut self, x: &Node);
            fn find_loop_in(&self, k: &TngComplexKey, allow_based: bool) -> Option<&TngComp>;
            fn deloop(&mut self, k: &TngComplexKey, c: &TngComp) -> Vec<TngComplexKey>;
            fn eliminate(&mut self, i: &TngComplexKey, j: &TngComplexKey);
            fn collect_keys<F, W>(&self, i: isize, pred: F, weight: W) -> Vec<(TngComplexKey, usize)>
                where F: Fn(&TngComplexKey) -> bool, W: Fn(&TngComplexKey) -> usize;
            fn current_step(&self) -> String;
            pub(crate) fn stat(&self) -> String;
        }
    }

    // Pivot cost for selection: off-axis pivots are eliminated in τ-pairs, so count ~2× the fill.
    fn pivot_weight(&self, k: &TngComplexKey) -> usize {
        let w = self.complex().vertex(k).c_weight();
        if self.key_map.is_sym(k) { w } else { 2 * w }
    }

    pub fn run(mut self) -> Self {
        info!("build config:\n{:#?}", self.config);
        info!("cutwidth profile:\n{}", self.profile_sym());

        if self.config.cut.enabled() {
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
        let planner = BuildPlanner::new(
            self.nodes(), self.tau_units(), &self.config.cut,
            self.config.node_order, self.complex().boundary_ends()
        );
        let plan = planner.plan();

        info!("chunk plan: {} pieces", plan.len());
        debug!("chunks: {:?}", plan.iter().map(|c| c.len()).collect_vec());

        let chunks = plan.into_iter().map(|chunk| {
            let built = Self::build_chunk(self.init_child(&chunk));
            (chunk, built)
        }).collect_vec();

        for (chunk, (c, key_map, elems)) in chunks {
            self.drop_nodes(|x| chunk.contains(x));
            self.merge(c, key_map, elems);
            info!("{} chunk merged: {}", self.current_step(), self.stat());
        }
    }

    // Run a child builder over its chunk and extract the reduced complex with its τ key-map
    // and transformed elements.
    fn build_chunk(child: Self) -> (TngComplex<R>, TauKeyMap, Vec<TngComplexElem<R>>) {
        info!("build chunk (n: {}): {}", child.nodes().len(), child.nodes().iter().join(", "));
        let child = child.run();
        info!("chunk built: {}", child.stat());
        
        let SymTngBuilder { key_map, mut inner, .. } = child;
        let elems = inner.elements_mut().take();
        (inner.into_tng_complex(), key_map, elems)
    }

    pub fn preprocess(&mut self) {
        SymTngPreprocessor::run(self);
    }

    pub fn process_nodes(&mut self) {
        info!("{} process {} nodes", self.current_step(), self.n_nodes());

        let planner = BuildPlanner::new(
            self.nodes(), self.tau_units(), &CutOption::None,
            self.config.node_order, self.complex().boundary_ends()
        );
        let Some(order) = planner.plan().pop() else {
            return;
        };
        debug!("node order: {}", order.iter().join(", "));

        for x in order {
            if !self.nodes().contains(&x) {
                continue; // already consumed as a τ-partner
            }
            let tx = self.inv_node(&x).clone();
            if x == tx {
                self.append_on_axis(&x);
            } else {
                self.append_off_axis(&x, &tx);
            }
        }
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

        self.merge(c, key_map, vec![]);
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

        self.merge(c, key_map, vec![]);
    }

    // Whether the automatic deloop runs. `false` for `Strategy::None` and, under
    // `no_full_deloop`, once all nodes are merged — remaining circles then defer to `into_raw_complex`.
    pub(crate) fn should_deloop(&self) -> bool {
        self.config.strategy.auto_deloop()
            && !(self.config.no_full_deloop && self.n_nodes() == 0)
    }

    pub(crate) fn merge(&mut self, c: TngComplex<R>, right_map: TauKeyMap, right_elements: Vec<TngComplexElem<R>>) {
        // build the merged τ key-map per degree (next to merge_vertices) rather than as one
        // up-front cartesian — for large knots that product never fits in memory.
        let left_map = std::mem::take(&mut self.key_map);
        let (left, right) = self.complex_mut().prepare_merge(c);
        let range = reachable_range(self.complex().h_range(), &self.config.h_range, self.n_nodes());

        // merge elements before delooping/eliminating, so the per-degree hooks transform them too.
        self.inner.elements_mut().merge(right_elements);
        
        debug!("{} merge {} <- {}", self.current_step(), left.stat(), right.stat());
        debug!("  key_map: {} × {}", left_map.len(), right_map.len());
        debug!("  merge range: {:?}", range);

        self.merge_incremental(&left, &right, range, &left_map, &right_map);
        self.prune_h_range();

        debug!("{} merged: {}", self.current_step(), self.stat());
    }

    // Per degree: deloop (when enabled), then (if the strategy eliminates) sweep i-2,i-1 by equivariant
    // Markowitz cost. Greedy also inline-eliminates during deloop; the sweep just catches what it
    // missed. Without delooping the sweep still applies — invertible pivots need no delooping,
    // shared circles pass through as cylinders.
    fn merge_incremental(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, range: RangeInclusive<isize>, left_map: &TauKeyMap, right_map: &TauKeyMap) {
        let top = *range.end();

        for i in range {
            debug!("{} build C[{i}]...", self.current_step());
            self.merge_slice(left, right, i, left_map, right_map);
            if self.should_deloop() {
                self.deloop_in(i - 1);
            }
            if self.config.strategy.auto_elim() {
                self.eliminate_in(i - 2);
                self.eliminate_in(i - 1);
            }
            debug!("{} built C[{i}]: {}", self.current_step(), self.complex().rank(i));
        }

        self.prune_isolated_top(top);
        if self.should_deloop() {
            self.deloop_in(top);
        }
        if self.config.strategy.auto_elim() {
            self.eliminate_in(top - 1);
        }
    }

    // Build degree `i`: the τ key-map slice, then its vertices and the edges into it.
    pub(crate) fn merge_slice(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, i: isize, left_map: &TauKeyMap, right_map: &TauKeyMap) {
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

    fn deloop_in(&mut self, i: isize) {
        self.deloop_in_with(i, false);
    }

    fn deloop_in_with(&mut self, i: isize, allow_based: bool) {
        let keys = self.collect_keys(i,
            |k| self.find_loop_in(k, allow_based).is_some(),
            |k| self.pivot_weight(k),
        );
        if keys.is_empty() { return }

        let total = keys.len();
        debug!("{} deloop in C[{i}]: {}, targets: {}", self.current_step(), self.complex().rank(i), total);

        let before = self.complex().rank(i) as isize;
        let mut done = 0;

        let mut pool = pivot_pool(keys);
        while let Some(k) = pop_min_pivot(&mut pool, |k|
            self.complex().contains_key(k).then(|| self.pivot_weight(k))
        ) {
            let Some(c) = self.find_loop_in(&k, allow_based).cloned() else { continue };
            let added = self.deloop_equiv(&k, &c);

            for nk in added {
                if self.find_loop_in(&nk, allow_based).is_some() {
                    let w = self.pivot_weight(&nk);
                    push_pivot(&mut pool, nk, w);
                }
            }
            
            done += 1;
            log_progress(Level::Debug, done, done - 1, done + pool.len(), PROGRESS_LOG_STEP, 2);
        }

        let after = self.complex().rank(i) as isize;

        debug!("{}   delooped C[{i}]: {} (delooped: {done}, diff: {})", self.current_step(), after, after - before);
        debug!("{}   neighbors: C[{}] {} / C[{}] {}",
            self.current_step(), i - 1, self.complex().rank(i - 1), i + 1, self.complex().rank(i + 1));
    }

    fn deloop_equiv(&mut self, k: &TngComplexKey, c: &TngComp) -> Vec<TngComplexKey> {
        let mut added = if self.key_map.is_sym(k) {
            if self.is_sym_comp(c) {
                // symmetric loop on symmetric key
                self.deloop_on_axis_sym(k, c)
            } else {
                // asymmetric loop on symmetric key
                self.deloop_on_axis_asym(k, c)
            }
        } else {
            // (symmetric or asymmetric) loop on asymmetric key
            self.deloop_off_axis(k, c)
        };

        // immediate elim eliminates each new vertex now; min-fill leaves them for the post-deloop
        // global pass, None leaves them entirely.
        if self.config.strategy.immediate_elim() {
            added.retain(|k|
                self.complex().contains_key(k) &&
                self.try_eliminate_equiv_at(k, ElimDir::Both) == 0
            );
            // an equivariant elim removes the whole τ-pair, so a key kept above may since have
            // been eliminated as another's τ-partner — re-retain so only live keys are returned.
            added.retain(|k| self.complex().contains_key(k));
        }
        added
    }

    fn deloop_on_axis_sym(&mut self, k: &TngComplexKey, c: &TngComp) -> Vec<TngComplexKey> {
        debug_assert!(self.key_map.is_sym(k));
        debug_assert!(self.is_sym_comp(c));

        let updated = self.deloop(k, c);

        self.key_map.remove(k);

        for &k_new in updated.iter() { 
            self.key_map.add_pair(k_new, k_new);
        }

        updated
    }

    #[allow(non_snake_case)]
    fn deloop_on_axis_asym(&mut self, k: &TngComplexKey, c: &TngComp) -> Vec<TngComplexKey> {
        debug_assert!(self.key_map.is_sym(k));
        debug_assert!(!self.is_sym_comp(c));
        debug_assert!(!c.is_marked());

        //          ⚪︎1 | ⚪︎1
        //  ⚪︎1 | ⚪︎X  <-->  ⚪︎X | ⚪︎1
        //          ⚪︎X | ⚪︎X

        let tc = c.convert_edges(|e| self.inv_edge(e));

        // deloop c, then tc on each surviving branch (the q-filter may have dropped some).
        self.deloop(k, c);
        for a in [KhAlgGen::X, KhAlgGen::I] {
            let ka = k + a;
            if self.complex().contains_key(&ka) {
                self.deloop(&ka, &tc);
            }
        }

        self.key_map.remove(k);

        // τ swaps the circles: τ(k+a+b) = k+b+a; τ preserves q ⇒ a pair survives or drops together.
        let added = iproduct!([KhAlgGen::X, KhAlgGen::I], [KhAlgGen::X, KhAlgGen::I])
            .map(|(a, b)| (k + a + b, k + b + a))
            .filter(|(ka, _)| self.complex().contains_key(ka))
            .collect_vec();

        added.into_iter().map(|(ka, kb)| {
            debug_assert!(self.complex().contains_key(&kb));
            self.key_map.add_pair(ka, kb);
            ka
        }).collect()
    }

    #[allow(non_snake_case)]
    fn deloop_off_axis(&mut self, k: &TngComplexKey, c: &TngComp) -> Vec<TngComplexKey> {
        debug_assert!(!self.key_map.is_sym(k));

        //  ⚪︎1 | ..  <-->  .. | ⚪︎1
        //  ⚪︎X | ..  <-->  .. | ⚪︎X

        let tk = *self.key_map.inv_key(k);
        let tc = c.convert_edges(|e| self.inv_edge(e));

        let mut ks = self.deloop(k, c);
        let mut tks = self.deloop(&tk, &tc);

        self.key_map.remove(k);

        for (&k_new, &tk_new) in Iterator::zip(ks.iter(), tks.iter()) { 
            self.key_map.add_pair(k_new, tk_new);
        }

        ks.append(&mut tks);
        ks
    }

    pub(crate) fn eliminate_in(&mut self, i: isize) {
        // pivot = equiv-invertible outgoing edge at its source (see the non-sym `eliminate_in`).
        let keys = self.collect_keys(i,
            |k| self.complex().vertex(k).out_edges()
                .any(|l| self.is_equiv_inv_edge(k, l)),
            |k| self.equiv_elim_cost(k, ElimDir::Outgoing),
        );
        if keys.is_empty() { return }

        // `targets` counts only the pivots the cap will actually eliminate (equiv cost ≤ cap); the
        // rest defer to the matrix. The sparkline shows the whole eliminatable distribution.
        let targets = keys.iter().filter(|(_, c)| !self.exceeds_elim_cap(*c)).count();
        
        debug!("{} eliminate in C[{i}]: {}, targets: {}", self.current_step(), self.complex().rank(i), targets);
        debug!("{}   fill: {}", self.current_step(), fill_cost_sparkline(&keys, self.config.max_elim_cost));

        let before = self.complex().rank(i) as isize;
        let mut pool = pivot_pool(keys);
        let mut done = 0;
        
        while let Some(k) = pop_min_pivot(&mut pool, |k|
            self.complex().contains_key(k).then(|| self.equiv_elim_cost(k, ElimDir::Outgoing))
        ) {
            // `pop_min_pivot` returns the cheapest pivot; once it exceeds the cap, so do all the
            // rest — stop and defer them (with the whole remaining frontier) to the matrix reduction.
            let cost = self.equiv_elim_cost(&k, ElimDir::Outgoing);
            if self.exceeds_elim_cap(cost) {
                debug!("{}   deferred {} pivots", self.current_step(), pool.len() + 1);
                break;
            }
            let n = self.try_eliminate_equiv_at(&k, ElimDir::Outgoing);
            if n > 0 {
                let prev = done;
                done += n; // off-axis events consume the pivot and its τ-mirror.
                log_progress(Level::Debug, done, prev, targets, PROGRESS_LOG_STEP, 2);
            }
        }

        let after = self.complex().rank(i) as isize;

        debug!("{}   eliminated C[{i}]: {} (diff: {})", self.current_step(), after, after - before);
        debug!("{}   neighbors: C[{}] {} / C[{}] {}",
            self.current_step(), i - 1, self.complex().rank(i - 1), i + 1, self.complex().rank(i + 1));
    }

    fn exceeds_elim_cap(&self, cost: usize) -> bool {
        self.config.max_elim_cost.is_some_and(|max| cost > max)
    }

    // Returns the number of degree-`i` pivot targets consumed: 1 on-axis, 2 off-axis (`k` and `τk`).
    // `Both` prefers the incoming side, matching the non-sym `try_eliminate_at`.
    fn try_eliminate_equiv_at(&mut self, k: &TngComplexKey, dir: ElimDir) -> usize {
        if matches!(dir, ElimDir::Incoming | ElimDir::Both) {
            if let Some(&j) = self.choose_equiv_inv_edge_into(k) {
                return self.eliminate_equiv(&j, k);
            }
        }
        if matches!(dir, ElimDir::Outgoing | ElimDir::Both) {
            if let Some(&l) = self.choose_equiv_inv_edge_from(k) {
                return self.eliminate_equiv(k, &l);
            }
        }
        0
    }

    fn eliminate_equiv(&mut self, i: &TngComplexKey, j: &TngComplexKey) -> usize {
        debug_assert_eq!(self.key_map.is_sym(i), self.key_map.is_sym(j));
        debug_assert!(self.complex().has_edge(i, j));

        let n = if self.key_map.is_sym(i) {
            self.eliminate(i, j);
            1
        } else {
            let ti = *self.key_map.inv_key(i);
            let tj = *self.key_map.inv_key(j);

            debug_assert!(self.complex().has_edge(&ti, &tj));

            self.eliminate(i, j);
            self.eliminate(&ti, &tj);
            2
        };

        self.key_map.remove(i);
        self.key_map.remove(j);
        n
    }

    // Markowitz cost of eliminating `k → l`; off-axis pivots eliminate in τ-pairs (~2× the fill).
    fn equiv_edge_weight(&self, k: &TngComplexKey, l: &TngComplexKey) -> usize {
        let w = self.complex().edge_weight(k, l);
        if self.key_map.is_sym(k) { w } else { 2 * w }
    }

    // Least Markowitz cost to eliminate `k`, over its equiv-invertible edges in the given
    // direction — the equivariant pivot priority (vs the cruder `pivot_weight`).
    fn equiv_elim_cost(&self, k: &TngComplexKey, dir: ElimDir) -> usize {
        let v = self.complex().vertex(k);
        let outs = || v.out_edges().filter(|l| self.is_equiv_inv_edge(k, l)).map(|l| self.equiv_edge_weight(k, l)).min();
        let ins = || v.in_edges().filter(|j| self.is_equiv_inv_edge(j, k)).map(|j| self.equiv_edge_weight(j, k)).min();
        let min = match dir {
            ElimDir::Outgoing => outs(),
            ElimDir::Incoming => ins(),
            ElimDir::Both => Iterator::chain(outs().into_iter(), ins()).min(),
        };
        min.unwrap_or(0)
    }

    // Cheapest invertible in-/out-edge within the cost cap (over-cap pivots are left for the matrix
    // pass — gates greedy's inline elim too). The cap is compared against the equiv fill cost.
    fn choose_equiv_inv_edge_into(&self, k: &TngComplexKey) -> Option<&TngComplexKey> {
        let cap = self.config.max_elim_cost;
        self.complex().vertex(k).in_edges().filter_map(|j|
            self.is_equiv_inv_edge(j, k).then_some(j)
        )
        .filter(|j| cap.map_or(true, |max| self.equiv_edge_weight(j, k) <= max))
        .min_by_key(|j| (self.complex().edge_weight(j, k), **j))
    }

    fn choose_equiv_inv_edge_from(&self, k: &TngComplexKey) -> Option<&TngComplexKey> {
        let cap = self.config.max_elim_cost;
        self.complex().vertex(k).out_edges().filter_map(|l|
            self.is_equiv_inv_edge(k, l).then_some(l)
        )
        .filter(|l| cap.map_or(true, |max| self.equiv_edge_weight(k, l) <= max))
        .min_by_key(|l| (self.complex().edge_weight(k, l), **l))
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
        if !self.should_deloop() {
            info!("{} skip finalize (deloop deferred): {}", self.current_step(), self.stat());
            return
        }

        if self.complex().is_completely_delooped() {
            info!("{} completely delooped: {}", self.current_step(), self.stat());
            return
        }

        info!("{} finalize: {}", self.current_step(), self.stat());

        self.deloop_all();

        // Deloop marked circles only when there are no other unmarked components left. 
        if self.complex().is_closed() {
            for i in self.complex().h_range() {
                self.deloop_in_with(i, true);
            }
        }

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

    // --- accessors for the cobordism-cone builder ---

    pub(crate) fn key_map(&self) -> &TauKeyMap {
        &self.key_map
    }

    pub(crate) fn key_map_mut(&mut self) -> &mut TauKeyMap {
        &mut self.key_map
    }

    pub(crate) fn inv_node(&self, x: &Node) -> &Node {
        &self.x_map[x]
    }

    pub(crate) fn inv_edge(&self, e: Edge) -> Edge {
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

    // The τ-units of the crossings: on-axis singly, off-axis as `(x, τx)` pairs kept once at
    // the lower index.
    fn tau_units(&self) -> Vec<Vec<usize>> {
        let nodes = self.nodes();
        let tau = &self.x_map;
        let idx_of: FxHashMap<Node, usize> = nodes.iter().enumerate()
            .map(|(i, x)| (x.clone(), i))
            .collect();

        (0..nodes.len()).filter_map(|i| {
            let j = idx_of[&tau[&nodes[i]]];
            match j {
                _ if j == i => Some(vec![i]),
                _ if i < j  => Some(vec![i, j]),
                _           => None,
            }
        }).collect()
    }

    /// Boundary-cutwidth profile of the symmetric (τ-equivariant) MinCut order: on-axis crossings
    /// singly, off-axis in `(x, τx)` pairs scored by combined cutwidth (shared axis edges cancel).
    pub(crate) fn profile_sym(&self) -> SymBuildProfile {
        let nodes = self.nodes();
        let units = self.tau_units();

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

    // A child builder over `chunk` (a sub-tangle), inheriting the parent's τ-maps and
    // `preprocess`/simplify settings; chunking always uses the MinCut order, never recursing.
    fn init_child(&self, chunk: &[Node]) -> Self {
        let (h, t) = self.inner.complex().ht();
        let base_pt = self.inner.complex().base_pt();

        // No per-chunk h_range: it would under-cover preprocess's off-axis key_map. Applied at the merge.
        // No no_full_deloop either: "final" means root-final — a chunk must deloop before the cross-chunk merge.
        // max_elim_cost: None so the root's "eliminate only free pivots at final" (max_elim_cost = 0) does NOT cascade — chunks eliminate fully.
        // No q_range either: q is exact only at the root's final merge (a chunk is an open sub-tangle).
        let config = SymBuildConfig { cut: CutOption::None, node_order: NodeOrder::MinCut, h_range: None, no_full_deloop: false, max_elim_cost: None, q_range: None, ..self.config.clone() };

        let mut inner = TngComplexBuilder::init(h, t, (0, 0), base_pt)
            .with_config(config.inner_build_config());
        inner.set_nodes(chunk.iter().cloned());
        inner.elements_mut().set(self.inner.elements().content().to_vec());

        let key_map = TauKeyMap::init();
        let real_top = inner.complex().deg_shift().0 + chunk.len() as isize; // child deg_shift = 0
        SymTngBuilder { inner, x_map: self.x_map.clone(), e_map: self.e_map.clone(), key_map, config, real_top }
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

        let elements = self.builder.inner.elements().content().to_vec();
        let (half, t_half) = self.partition_off_axis();

        info!("{} preprocess off-axis: {} + {}", self.builder.current_step(), half.len(), t_half.len());

        let (c, tc, key_map, e_half, e_t_half) = self.build_from_half(&half, elements);

        // merge the half (combines with the resident seeds = e_half), then its τ-image (completes them).
        self.merge_half(&half, c, e_half);
        self.merge_half(&t_half, tc, e_t_half);

        self.builder.key_map = key_map;

        info!("{} preprocess done: {}", self.builder.current_step(), self.builder.stat());
    }

    fn merge_half(&mut self, nodes: &[Node], c: TngComplex<R>, elements: Vec<TngComplexElem<R>>) {
        self.builder.inner.drop_nodes(|x| nodes.contains(x));
        self.builder.inner.merge(c, elements);
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

    fn build_from_half(&self, crossings: &[Node], elements: Vec<TngComplexElem<R>>) -> (TngComplex<R>, TngComplex<R>, TauKeyMap, Vec<TngComplexElem<R>>, Vec<TngComplexElem<R>>) {
        // crossings appended after this chunk: all remaining nodes minus the chunk (half + τ-half).
        let r_rest = self.builder.n_nodes() - 2 * crossings.len();

        let n = elements.len();

        // Transport each element with its τ-image; applying τ again to (τe)_h gives the τ-half element
        // (e_h's state/in_cob via τ² = id, τ-converted out_cob), which the bilinear merge then completes.
        let t_elements = elements.iter().map(|e| self.tau_element(e)).collect_vec();

        let mut b = self.half_builder();
        b.set_nodes(crossings.iter().cloned());
        b.elements_mut().set(elements.into_iter().chain(t_elements));
        b.process_nodes();

        info!("{} half complex built: {}", self.builder.current_step(), b.stat());

        let keys = b.complex().keys().cloned().collect_vec();
        let mut e_half = b.elements_mut().take();
        let tau_half = e_half.split_off(n);
        let e_t_half = tau_half.iter().map(|e| self.tau_element(e)).collect_vec();

        debug!("  build complexes...");
        let c = b.into_tng_complex();
        let tc = c.convert_edges(|e| self.builder.inv_edge(e));

        debug!("  pair key_map ({}² entries)...", keys.len());
        let key_map = TauKeyMap::from_half(&keys, self.weight_band(r_rest));

        (c, tc, key_map, e_half, e_t_half)
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

    // The τ-image of an element: mirror the state (inv_node), in_cob/out_cob/base_pt (inv_edge).
    fn tau_element(&self, e: &TngComplexElem<R>) -> TngComplexElem<R> {
        let tau_cob = |f: &LcCob<R>| f.map_ref(|c, r|
            (c.convert_edges(|e| self.builder.inv_edge(e)), r.clone())
        );
        let state = e.state().iter().map(|(x, b)|
            (self.builder.inv_node(x).clone(), *b)
        ).collect();
        let base_pt = e.base_pt().map(|b| self.builder.inv_edge(b));
        let in_cob = e.in_cob().convert_edges(&|e| self.builder.inv_edge(e));

        let mut te = TngComplexElem::new(state, in_cob, base_pt);
        te.set_out_cob(e.out_cob().iter().map(|(k, f)| (*k, tau_cob(f))));
        te
    }
}

#[cfg(test)]
#[allow(unused)]
mod tests {
    use crate::khi::{KhIGen, KhIHomology};

    use super::*;
    use num_traits::Zero;

    use yui_core::ext::IteratorExt;
    use yui_core::lc::Lc;
    use yui_core::num::FF2;
    use yui_core::poly::Poly;
    use yui_core::ext::RangeExt;
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
                .with_config(BuildConfig { strategy: Strategy::None, cut: CutOption::None, ..Default::default() });

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
    fn strategies_agree() {
        // every strategy on the sym builder must agree on homology.
        let l = InvLink::from_symmetric_pd_code(
            [[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]
        );
        let (h, t) = (FF2::zero(), FF2::zero());
        let build = |strategy| {
            let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
            b.config.strategy = strategy;
            b.run().into_tng_complex().into_raw_complex()
        };

        let ref_c = build(Strategy::Greedy);
        let range = ref_c.support().cloned().range().unwrap();
        let ref_h = ref_c.homology();
        for strategy in [Strategy::MinFill, Strategy::NoElim, Strategy::None] {
            let c = build(strategy);
            c.check_d_all();
            let h = c.homology();
            for i in range.clone() {
                assert_eq!(h[i].rank(), ref_h[i].rank(), "rank at {i}, {strategy:?}");
            }
        }
    }

    #[test]
    fn no_full_deloop_agrees() {
        // no_full_deloop must not change homology: into_raw_complex redoes the deferred deloop.
        let l = InvLink::test_data("6_3");
        let (h, t) = (FF2::zero(), FF2::zero());
        let build = |skip| {
            let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
            b.config.no_full_deloop = skip;
            b.run().into_tng_complex().into_raw_complex()
        };

        let ref_c = build(false);
        let range = ref_c.support().cloned().range().unwrap();
        let ref_h = ref_c.homology();
        let c = build(true);
        c.check_d_all();
        let h = c.homology();
        for i in range {
            assert_eq!(h[i].rank(), ref_h[i].rank(), "rank at {i}");
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
            b.config.cut = chunks.map_or(CutOption::None, CutOption::Auto);
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
    fn chunk_build_reduced_matches() {
        // reduced (based) chunked build must reproduce the non-chunked reduced homology — the
        // based component is capped once, by the top-level finalize, not per chunk.
        let l = InvLink::from_symmetric_pd_code(
            [[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]
        ).with_base_pt(1); // edge 1 is on-axis (inv(1) = 1)
        let (h, t) = (FF2::zero(), FF2::zero());

        let build = |chunks: Option<usize>| {
            let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, true);
            b.config.cut = chunks.map_or(CutOption::None, CutOption::Auto);
            b.run().into_tng_complex().into_raw_complex()
        };

        let normal = build(None);
        let range = normal.support().cloned().range().unwrap();
        let hn = normal.homology();

        for k in [Some(2), Some(3)] {
            let hc = build(k).homology();
            for i in range.clone() {
                assert_eq!(hc[i].rank(), hn[i].rank(), "reduced rank at {i} (chunks {k:?})");
            }
        }
    }

    // chunked builds must track the canon cycles too: the Lee-class divisibility (the ssi
    // ingredient) computed from each chunked KhI homology must match the non-chunked one.
    #[test]
    fn chunk_elements_match() {
        use crate::ss::div_vec;
        type P = Poly<'H', FF2>;

        let l = InvLink::from_symmetric_pd_code(
            [[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]
        );
        let c = P::variable();
        let (h, t) = (c.clone(), P::zero());

        let div = |chunks: Option<usize>| {
            let config = SymBuildConfig { cut: chunks.map_or(CutOption::None, CutOption::Auto), ..Default::default() };
            let kh = KhIHomology::new_with_config(&l, &h, &t, false, config);
            kh.canon_cycles().iter().map(|z| {
                let v = kh[kh.h_deg_of_chain(z)].vectorize_euc(z);
                div_vec(&v.subvec(0..2), &c).unwrap()
            }).collect_vec()
        };

        let ref_d = div(None);
        for k in [Some(2), Some(3)] {
            assert_eq!(div(k), ref_d, "divisibility, chunks {k:?}");
        }
    }

    #[test]
    fn test_kh_3_1() {
        let l = InvLink::test_data("3_1");
        let (h, t) = (FF2::zero(), FF2::zero());

        let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
        b.config.strategy = Strategy::None;
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
        b.config.strategy = Strategy::None;
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
        // Strategy::None never deloops (not even in finalize); into_raw_complex expands the
        // remaining circles at the matrix level, giving the same generators.
        let l = InvLink::test_data("3_1");
        let (h, t) = (FF2::zero(), FF2::zero());

        let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
        b.config.strategy = Strategy::None;
        b.process_nodes();

        assert!(!b.inner.complex().is_completely_delooped());

        b.finalize();

        assert!(!b.inner.complex().is_completely_delooped());

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
        b.config.strategy = Strategy::NoElim;
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