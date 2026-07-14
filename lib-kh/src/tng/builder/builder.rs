//! Incremental builder for [`TngComplex`]: scan crossings one at a time,
//! tensor-merge with the new crossing's small complex, then deloop newborn
//! circles and gauss-eliminate invertible edges to keep the complex small.
//!
//! References:
//! - BN05 — D. Bar-Natan, "Khovanov's homology for tangles and cobordisms",
//!   Geom. Topol. 9 (2005), 1443–1499.
//!   <https://doi.org/10.2140/gt.2005.9.1443>, <https://arxiv.org/abs/math/0410495>
//! - BN07 — D. Bar-Natan, "Fast Khovanov homology computations",
//!   J. Knot Theory Ramif. 16 (2007), 243–255.
//!   <https://doi.org/10.1142/S0218216507005294>, <https://arxiv.org/abs/math/0606318>

use std::fmt;
use std::ops::RangeInclusive;

use rustc_hash::FxHashSet;
use itertools::Itertools;
use log::{debug, info, trace, log_enabled, Level};
use yui_core::{Ring, RingOps};
use yui_link::{Node, Edge, Link};

use yui_homology::ChainComplex1;

use crate::kh::{KhChain, KhComplex, KhGen};
use crate::tng::{MAX_EDGE, TngComp, TngComplexElem, LcCobTrait, TngComplex, TngComplexKey};
use super::{reachable_range, pop_min_pivot, pivot_pool, push_pivot, sparkline, fill_cost_sparkline, cutwidth_after, toggle_boundary, boundary_edges, select_cuts, cut_components, merge_order, TngElemBuilder};

// Progress logging for the long per-op build loops (eliminate / deloop / asymmetric elimination):
// emit a line every `PROGRESS_LOG_STEP` ops, but only for rounds larger than `PROGRESS_LOG_MIN`.
pub(super) const PROGRESS_LOG_STEP: usize = 20_000;
pub(super) const PROGRESS_LOG_MIN: usize = 50_000;

/// How the next crossing to append is chosen.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum NodeOrder {
    #[default]
    MinCut, // minimize the boundary cutwidth (default; bounds dense-slice memory, wins on wide knots)
    Given,  // process crossings in the given (PD) order — no reordering
}

/// How the complex is simplified while building.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum BuildMode {
    #[default]
    Greedy,    // deloop every circle, eliminate immediately
    MinFill,   // deloop a whole degree, then eliminate by global min-fill (Markowitz)
    NoElim,    // deloop every circle but don't eliminate (delooped, unreduced complex)
    None,      // don't deloop, don't eliminate (raw merge; finalize still deloops to a valid complex)
}

impl BuildMode {
    // whether the build deloops at all (None = raw merge, deloop deferred to finalize).
    pub fn auto_deloop(&self) -> bool {
        *self != BuildMode::None
    }

    // whether the build eliminates at all (Greedy inline, MinFill swept; NoElim/None don't).
    pub fn auto_elim(&self) -> bool {
        matches!(self, BuildMode::Greedy | BuildMode::MinFill)
    }

    // whether each newly-delooped vertex is eliminated inline (vs swept after).
    pub fn immediate_elim(&self) -> bool {
        *self == BuildMode::Greedy
    }
}

/// Divide-and-conquer chunking: `None` = single pass; `Auto(k)` cuts the MinCut order at its `k-1`
/// deepest cutwidth valleys; `Manual` severs the given edge-cut(s) (a union) into pieces.
#[derive(Clone, Debug, Default)]
pub enum CutOption {
    #[default]
    None,
    Auto(usize),
    // cut after the unit positions closest to the given cumulative crossing counts —
    // direct control over chunk balance (Auto cuts only at cutwidth valleys).
    AtCrossings(Vec<usize>),
    Manual(Vec<Vec<Edge>>),
}

impl CutOption {
    pub fn enabled(&self) -> bool {
        !matches!(self, CutOption::None)
    }
}

/// Toggles for the automatic simplification done while building.
#[derive(Clone, Debug)]
pub struct BuildConfig {
    pub node_order: NodeOrder,
    pub mode: BuildMode,
    // divide-and-conquer chunking (auto cutwidth or manual edge-cuts); None = single pass.
    pub cut: CutOption,
    pub h_range: Option<RangeInclusive<isize>>,
    // drop generators outside this q-range. Only applied once the diagram is closed (see
    // `should_drop`), where q-degrees are exact. `None` = no q-truncation.
    pub q_range: Option<RangeInclusive<isize>>,
    // skip eliminations whose fill cost (`edge_weight` = Schur block size) exceeds this; the
    // survivors defer to the matrix reduction. `None` = eliminate everything (current behavior).
    pub max_elim_cost: Option<usize>,
    // skip the final deloop (the last merge and `finalize`): remaining circles are deferred to
    // `into_raw_complex`'s matrix-level expansion + the `ChainReducer`. See `should_deloop`.
    pub no_full_deloop: bool,
}

impl Default for BuildConfig {
    fn default() -> Self {
        Self { node_order: NodeOrder::default(), mode: BuildMode::default(), cut: CutOption::None, h_range: None, q_range: None, max_elim_cost: None, no_full_deloop: false }
    }
}


pub struct TngComplexBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    complex: TngComplex<R>,
    nodes: Vec<Node>,
    loops: Vec<Edge>,
    elements: TngElemBuilder<R>,
    config: BuildConfig,
}

impl<R> TngComplexBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_link(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        Self::assert_max_edge(l);
        let base_pt = if reduced { l.base_pt() } else { None };
        let deg_shift = KhComplex::deg_shift_for(l, reduced);

        let mut b = Self::init(h, t, deg_shift, base_pt);
        b.set_nodes(l.nodes().cloned());
        b.set_loops(l.loops().iter().cloned());

        if t.is_zero() && l.is_knot() {
            let canon = TngComplexElem::canon_cycles(l, base_pt);
            b.elements_mut().set(canon);
        }

        b
    }

    // note: the `EdgeSet` bitmap wraps silently on overflow in release builds — fail loudly up front.
    fn assert_max_edge(l: &Link) {
        if let Some(e) = l.edges().into_iter().max() {
            assert!(e <= MAX_EDGE, "edge label {e} exceeds the EdgeSet capacity ({MAX_EDGE}); enable the `big-link` feature");
        }
    }

    pub fn init(h: &R, t: &R, deg_shift: (isize, isize), base_pt: Option<Edge>) -> Self { 
        let complex = TngComplex::init(h, t, deg_shift, base_pt);
        Self {
            complex,
            nodes: vec![],
            loops: vec![],
            elements: TngElemBuilder::new(),
            config: BuildConfig::default(),
        }
    }

    pub fn with_config(mut self, config: BuildConfig) -> Self {
        // drop canon cycles whose h-degree falls outside the range (computed, not assumed h0).
        if let Some(range) = &config.h_range {
            let shift = self.complex.deg_shift().0;
            self.elements.retain(|e| range.contains(&(shift + e.rel_h_deg())));
        }
        self.config = config;
        self
    }

    pub(crate) fn from_tng_complex(complex: TngComplex<R>, config: BuildConfig) -> Self {
        Self { complex, nodes: vec![], loops: vec![], elements: TngElemBuilder::new(), config }
    }

    pub fn config(&self) -> &BuildConfig {
        &self.config
    }

    pub fn complex(&self) -> &TngComplex<R> {
        &self.complex
    }

    pub(crate) fn complex_mut(&mut self) -> &mut TngComplex<R> {
        &mut self.complex
    }

    pub fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn set_nodes<I>(&mut self, nodes: I)
    where I: IntoIterator<Item = Node> {
        self.nodes = nodes.into_iter().collect_vec();
    }

    pub(crate) fn drop_nodes<F>(&mut self, pred: F)
    where F: Fn(&Node) -> bool {
        self.nodes.retain(|x| !pred(x));
    }

    pub fn loops(&self) -> &[Edge] {
        &self.loops
    }

    pub fn set_loops<I>(&mut self, loops: I)
    where I: IntoIterator<Item = Edge> {
        self.loops = loops.into_iter().collect_vec();
    }

    pub(crate) fn elements(&self) -> &TngElemBuilder<R> {
        &self.elements
    }

    pub(crate) fn elements_mut(&mut self) -> &mut TngElemBuilder<R> {
        &mut self.elements
    }

    /// Keys at degree `i` matching `pred`, each paired with `weight(k)`
    pub(crate) fn collect_keys<F, W>(&self, i: isize, pred: F, weight: W) -> Vec<(TngComplexKey, usize)>
    where F: Fn(&TngComplexKey) -> bool, W: Fn(&TngComplexKey) -> usize {
        self.complex.keys_of_deg(i)
            .filter(|k| pred(k))
            .map(|k| (*k, weight(k)))
            .collect_vec()
    }

    // "(committed/total)" crossing progress, for log prefixes.
    pub(crate) fn current_step(&self) -> String {
        format!("({}/{})", self.complex.dim(), self.complex.dim() + self.n_nodes())
    }

    pub fn run(mut self) -> Self {
        info!("build config:\n{:#?}", self.config);
        info!("cutwidth profile:\n{}", self.profile());
        if self.config.cut.enabled() {
            self.process_chunks();
        } else {
            self.process_nodes();
        }
        self.process_free_loops();
        self.finalize();
        self
    }

    fn process_chunks(&mut self) {
        let chunks = ChunkBuilder { builder: self }.build_chunks();
        for (chunk, (c, elems)) in chunks {
            self.drop_nodes(|x| chunk.contains(x));
            self.merge(c, elems);
            info!("{} chunk merged: {}", self.current_step(), self.stat());
        }
    }

    // See [BN07, §7] (scan-and-cancel algorithm).
    pub fn process_nodes(&mut self) {
        info!("{} process {} nodes", self.current_step(), self.n_nodes());

        while let Some(x) = self.choose_next_node().cloned() {
            self.append_node(&x)
        }
    }

    /// Pick the next node by node order, ties broken by earliest crossing order
    /// (`self.nodes` keeps PD order — `prepare_append` removes via order-preserving `Vec::remove`).
    pub(crate) fn choose_next_node(&self) -> Option<&Node> {
        self.nodes.iter().enumerate()
            .min_by_key(|(i, x)| {
                let score = match self.config.node_order {
                    NodeOrder::MinCut => self.cutwidth(x),
                    NodeOrder::Given => 0, // constant → ties broken by earliest index = given order
                };
                (score, *i)
            })
            .map(|(_, x)| x)
    }

    /// Boundary cutwidth (open-edge count) after appending `x`. `boundary_ends` is cheap, so
    /// recomputing it per call is fine.
    pub(crate) fn cutwidth(&self, x: &Node) -> isize {
        let open: FxHashSet<Edge> = self.complex.boundary_ends().collect();
        cutwidth_after(&open, &[x]) as isize
    }

    pub fn append_node(&mut self, x: &Node) {
        info!("{} append: {x}", self.current_step());

        self.prepare_append(x);
        
        let (h, t) = self.complex.ht();
        let cx = TngComplex::from_node(h, t, x, self.complex.base_pt());
        self.merge(cx, vec![]);
    }

    pub(crate) fn prepare_append(&mut self, x: &Node) { 
        if let Some(i) = self.nodes.iter().find_position(|&e| e == x) { 
            self.nodes.remove(i.0);
        }

        self.elements.append_node(x);
    }

    // Whether the automatic deloop runs. `false` for `BuildMode::None` and, under
    // `no_full_deloop`, once all nodes are merged — remaining circles then defer to `into_raw_complex`.
    pub(crate) fn should_deloop(&self) -> bool {
        self.config.mode.auto_deloop()
            && !(self.config.no_full_deloop && self.n_nodes() == 0)
    }

    pub fn merge(&mut self, other: TngComplex<R>, other_elements: Vec<TngComplexElem<R>>) {
        let (left, right) = self.complex.prepare_merge(other);
        let range = reachable_range(self.complex.h_range(), &self.config.h_range, self.n_nodes());

        // merge elements before delooping/eliminating, so the per-degree hooks transform them too.
        self.elements.merge(other_elements);

        debug!("{} merge {} <- {}", self.current_step(), left.stat(), right.stat());
        debug!("  merge range: {:?}", range);

        if self.should_deloop() {
            self.merge_incremental(&left, &right, range);
        } else {
            self.complex.merge_with(&left, &right);
        }

        self.prune_h_range();
        debug!("{} merged: {}", self.current_step(), self.stat());
    }

    // Per degree: deloop, then (if the mode eliminates) sweep i-2,i-1 by Markowitz cost.
    // Greedy also inline-eliminates during deloop; the sweep just catches what it missed.
    fn merge_incremental(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, range: RangeInclusive<isize>) {
        debug_assert!(self.config.mode.auto_deloop()); // None is dispatched to merge_with
        let top = *range.end();

        for i in range {
            debug!("{} build C[{i}]...", self.current_step());
            self.merge_slice(left, right, i);
            self.deloop_in(i - 1);
            if self.config.mode.auto_elim() {
                self.eliminate_in(i - 2);
                self.eliminate_in(i - 1);
            }
            debug!("{} built C[{i}]: {}", self.current_step(), self.complex.rank(i));
        }

        self.prune_isolated_top(top);
        self.deloop_in(top);
        if self.config.mode.auto_elim() {
            self.eliminate_in(top - 1);
        }
    }

    // Build degree `i`: merge in its vertices and the edges into it.
    pub(super) fn merge_slice(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, i: isize) {
        if log_enabled!(Level::Debug) {
            let (nv, ne) = slice_estimate(left, right, i);
            // no sources below the already-built bottom — the raw estimate would overcount there.
            let ne = if self.complex.rank(i - 1) == 0 { 0 } else { ne };
            debug!("  merge C[{i}]: ({nv} verts, ~{ne} edges)");
        }
        let nv = self.complex.merge_vertices(left, right, i);
        debug!("  +{nv} verts");
        let ne = self.complex.merge_edges(left, right, i - 1);
        debug!("  +{ne} edges");

        // drop merged vertices already outside q_range (before we deloop/eliminate them).
        self.prune_q_range(i);
    }

    // Drop vertices in degree `i` whose q-degree can't reach `config.q_range`.
    fn prune_q_range(&mut self, i: isize) {
        if self.config.q_range.is_none() {
            return;
        }
        let doomed = self.complex.keys_of_deg(i).filter(|k| self.should_drop(k)).copied().collect_vec();
        if !doomed.is_empty() {
            debug!("  -{} verts (q_range)", doomed.len());
        }
        self.complex.remove_vertices(&doomed);
    }

    /// Drop vertices that can't end up in `config.h_range`: degree `d` ends in
    /// `[d, d + r]` (`r` = pending crossings), so doomed iff `d > b` or `d + r < a`.
    fn prune_h_range(&mut self) {
        let Some(h_range) = self.config.h_range.clone() else { return };
        let i0 = self.complex.deg_shift().0;
        let r = self.n_nodes() as isize;

        // a vertex of degree `d` reaches `[d, d + r]`, so it stays relevant iff
        // `d ∈ [a - r, b]` — current degrees that can still land in `h_range`.
        let live = (*h_range.start() - r) ..= *h_range.end();

        let doomed = self.complex.keys_of(|k|
            !live.contains(&(k.weight() as isize + i0))
        ).copied().collect_vec();

        if !doomed.is_empty() {
            debug!("prune {} verts outside h_range.", doomed.len());
        }

        self.complex.remove_vertices(&doomed);
    }

    // Drop no-in-edge vertices at a TRUNCATED window-top (`top < real_top`): they feed only the
    // discarded `top+1` homology. At the real top they're genuine generators, so skip.
    fn prune_isolated_top(&mut self, top: isize) {
        // real top = deg_shift + total crossings (dim + remaining nodes).
        let real_top = self.complex.deg_shift().0 + (self.complex.dim() + self.n_nodes()) as isize;
        let truncated = self.config.h_range.as_ref().is_some_and(|w| top == *w.end()) && top < real_top;
        if !truncated { return; }

        let doomed = self.complex.keys_of_deg(top)
            .filter(|k| self.complex.vertex(k).in_edges().next().is_none())
            .copied()
            .collect_vec();
        if !doomed.is_empty() {
            debug!("prune {} isolated verts in C[{top}].", doomed.len());
        }
        self.complex.remove_vertices(&doomed);
    }

    // The first unmarked (or based, if `allow_based`) circle in `k`'s tangle.
    pub(crate) fn find_loop_in(&self, k: &TngComplexKey, allow_based: bool) -> Option<&TngComp> {
        let v = self.complex.vertex(k);
        v.tng().comps()
            .find(|c| c.is_circle() && (allow_based || !c.is_marked()))
    }

    // Deloop unmarked loops over all degrees.
    pub fn deloop_all(&mut self) {
        for i in self.complex.h_range() {
            self.deloop_in(i);
        }
    }

    pub fn deloop_in(&mut self, i: isize) {
        self.deloop_in_with(i, false);
    }

    pub fn deloop_in_with(&mut self, i: isize, allow_based: bool) {
        let keys = self.collect_keys(i,
            |k| self.find_loop_in(k, allow_based).is_some(),
            |k| self.complex.vertex(k).c_weight(),
        );
        if keys.is_empty() { return }

        let total = keys.len();
        debug!("{} deloop in C[{i}]: {}, targets: {}", self.current_step(), self.complex.rank(i), total);

        let before = self.complex.rank(i) as isize;

        let mut pool = pivot_pool(keys);
        let (mut done, mut elim) = (0, 0);
        while let Some(k) = pop_min_pivot(&mut pool, |k|
            self.complex.contains_key(k).then(|| self.complex.vertex(k).c_weight())
        ) {
            let Some(&c) = self.find_loop_in(&k, allow_based) else { continue };

            let new_keys = self.deloop(&k, &c);
            // a normal circle yields 2 branches; a based circle yields 1, and in greedy mode each
            // branch may be eliminated on the spot. Count branches short of 2 as eliminated, so
            // `delooped - eliminated = diff` holds in both greedy and min-fill.
            elim += 2usize.saturating_sub(new_keys.len());
            for new_key in new_keys {
                if self.find_loop_in(&new_key, allow_based).is_some() {
                    let w = self.complex.vertex(&new_key).c_weight();
                    push_pivot(&mut pool, new_key, w);
                }
            }
            done += 1;
            if done % PROGRESS_LOG_STEP == 0 {
                debug!("{}   ... delooped {done} ({}% eliminated, remain: {})", self.current_step(), elim * 100 / done, pool.len());
            }
        }

        let after = self.complex.rank(i) as isize;

        debug!("{}   delooped C[{i}]: {} (delooped: {done}, eliminated: {elim}, diff: {})", self.current_step(), after, after - before);
        debug!("{}   neighbors: C[{}] {} / C[{}] {}",
            self.current_step(), i - 1, self.complex.rank(i - 1), i + 1, self.complex.rank(i + 1));
    }

    pub fn deloop(&mut self, k: &TngComplexKey, c: &TngComp) -> Vec<TngComplexKey> {
        trace!("{} deloop {c} in {}", self.stat(), self.complex.vertex(k));

        self.elements.deloop(k, c);

        let mut added = self.complex.deloop(k, c);

        // drop delooped branches outside `config.q_range` (exact once closed, see `should_drop`).
        // The sym builder's paired deloop tolerates this — it guards its nested deloops on existence.
        if self.config.q_range.is_some() {
            let (keep, doomed): (Vec<_>, Vec<_>) = added.into_iter().partition(|k| !self.should_drop(k));
            self.complex.remove_vertices(&doomed);
            added = keep;
        }

        // immediate elim eliminates each new vertex now; min-fill leaves them for the post-deloop
        // global pass, None leaves them entirely. `try_eliminate_at` skips over-cap pivots.
        if self.config.mode.immediate_elim() {
            // retain only the keys that weren't eliminated
            added.retain(|k| self.try_eliminate_at(k).is_none());
        }
        added
    }

    // A delooped branch is doomed if its q-degree can't land in `config.q_range`. q is exact only
    // once the diagram is closed (no open arcs); then each remaining circle shifts q by ±1.
    pub(crate) fn should_drop(&self, k: &TngComplexKey) -> bool {
        let Some(q_range) = self.config.q_range.as_ref() else {
            return false;
        };
        if !self.complex.is_closed() {
            return false;
        }

        let q0 = self.complex.deg_shift().1 + k.as_gen().rel_q_deg();
        let nc = self.complex.vertex(k).tng().comps().filter(|c| c.is_circle()).count() as isize;

        q0 + nc < *q_range.start() || q0 - nc > *q_range.end()
    }

    pub fn eliminate_in(&mut self, i: isize) {
        let keys = self.collect_keys(i,
            |k| self.complex.vertex(k).out_edges().any(|l|
                self.complex.edge(k, l).is_invertible()
            ),
            |k| self.complex.elim_cost(k),
        );
        if keys.is_empty() { return }

        // `targets` counts only the pivots the cap will actually eliminate (cost ≤ cap); the rest
        // defer to the matrix. The sparkline shows the *whole* eliminatable distribution for context.
        let targets = self.config.max_elim_cost
            .map_or(keys.len(), |max| keys.iter().filter(|(_, c)| *c <= max).count());
        debug!("{} eliminate in C[{i}]: {}, targets: {}", self.current_step(), self.complex.rank(i), targets);
        debug!("{}   fill: {}", self.current_step(), fill_cost_sparkline(&keys, self.config.max_elim_cost));

        let before = self.complex.rank(i) as isize;

        let mut pool = pivot_pool(keys);
        let mut done = 0;
        while let Some(k) = pop_min_pivot(&mut pool, |k|
            self.complex.contains_key(k).then(|| self.complex.elim_cost(k))
        ) {
            // `pop_min_pivot` returns the cheapest pivot; once it exceeds the cap, so do all the
            // rest — stop and defer them (with the whole remaining frontier) to the matrix reduction.
            if let Some(max) = self.config.max_elim_cost {
                let cost = self.complex.elim_cost(&k);
                if cost > max {
                    debug!("{}   deferred {} pivots to matrix (min cost 2^{} > cap {max})", self.current_step(), pool.len() + 1, cost.ilog2());
                    break;
                }
            }
            if self.try_eliminate_at(&k).is_some() {
                done += 1;
                if targets > PROGRESS_LOG_MIN && done % PROGRESS_LOG_STEP == 0 {
                    debug!("{}   ... eliminated {done}/{targets} in C[{i}] (rank: {})", self.current_step(), self.complex.rank(i));
                }
            }
        }

        let after = self.complex.rank(i) as isize;

        debug!("{}   eliminated C[{i}]: {} (diff: {})", self.current_step(), after, after - before);
        debug!("{}   neighbors: C[{}] {} / C[{}] {}",
            self.current_step(), i - 1, self.complex.rank(i - 1), i + 1, self.complex.rank(i + 1));
    }

    // Eliminate at `k` via an invertible in- or out-edge; returns the fill cost paid
    // (`edge_weight` of the chosen pivot = Schur block size), or `None` if nothing to eliminate.
    pub fn try_eliminate_at(&mut self, k: &TngComplexKey) -> Option<usize> {
        if let Some(&j) = self.choose_inv_edge_into(&k) {
            let cost = self.complex.edge_weight(&j, k);
            self.eliminate(&j, &k);
            Some(cost)
        } else if let Some(&l) = self.choose_inv_edge_from(&k) {
            let cost = self.complex.edge_weight(k, &l);
            self.eliminate(&k, &l);
            Some(cost)
        } else {
            None
        }
    }

    pub fn eliminate(&mut self, i: &TngComplexKey, j: &TngComplexKey) {
        trace!("{} eliminate {}: {} -> {}", self.stat(), self.complex.edge(i, j), self.complex.vertex(i), self.complex.vertex(j));
        
        self.elements.eliminate(&self.complex, i, j);
        self.complex.eliminate(i, j);
    }

    // Cheapest invertible in-/out-edge within the cost cap (over-cap pivots are left for the matrix
    // pass — this is what gates greedy's inline elim as well as the `eliminate_in` sweep).
    fn choose_inv_edge_into(&self, k: &TngComplexKey) -> Option<&TngComplexKey> {
        let cap = self.config.max_elim_cost;
        self.complex.vertex(k).in_edges().filter_map(|j|
            self.complex.edge(j, k).is_invertible().then_some(j)
        )
        .filter(|j| cap.map_or(true, |max| self.complex.edge_weight(j, k) <= max))
        .min_by_key(|j| (self.complex.edge_weight(j, k), **j))
    }

    fn choose_inv_edge_from(&self, k: &TngComplexKey) -> Option<&TngComplexKey> {
        let cap = self.config.max_elim_cost;
        self.complex.vertex(k).out_edges().filter_map(|l|
            self.complex.edge(k, l).is_invertible().then_some(l)
        )
        .filter(|l| cap.map_or(true, |max| self.complex.edge_weight(k, l) <= max))
        .min_by_key(|l| (self.complex.edge_weight(k, l), **l))
    }

    pub fn process_free_loops(&mut self) {
        while !self.loops.is_empty() { 
            let c = self.loops.remove(0);

            self.elements.insert_loop(c);

            let (h, t) = self.complex.ht();
            let marked = self.complex.base_pt() == Some(c);
            let c = TngComplex::from_loop(h, t, c, marked);
            self.merge(c, vec![]);

            if self.config.mode.auto_deloop() {
                self.deloop_all();
            }
        }
    }

    fn finalize(&mut self) {
        if !self.should_deloop() {
            info!("{} skip finalize (deloop deferred): {}", self.current_step(), self.stat());
            return;
        }

        if self.complex.is_completely_delooped() {
            info!("{} completely delooped: {}", self.current_step(), self.stat());
            return;
        }

        info!("{} finalize: {}", self.current_step(), self.stat());

        self.deloop_all();

        // Deloop marked circles only when there are no other unmarked components left. 
        if self.complex.is_closed() {
            for i in self.complex.h_range() {
                self.deloop_in_with(i, true);
            }
        }

        info!("{} finalized: {}", self.current_step(), self.stat());
    }

    pub fn into_tng_complex(self) -> TngComplex<R> {
        self.complex
    }

    /// Convert to the raw complex, applying `config.q_range` (matrix-level deloop filter). This is the
    /// only q-filter on the `no_full_deloop` path, where circles expand into generators here.
    pub fn into_raw_complex(self) -> ChainComplex1<KhGen, R> {
        match self.config.q_range.clone() {
            Some(range) => self.complex.into_raw_complex_filtered(range),
            None => self.complex.into_raw_complex(),
        }
    }

    pub fn eval_elements(&self) -> Vec<KhChain<R>> {
        let (h, t) = self.complex.ht();
        // an un-delooped complex (`no_full_deloop`) needs the circle-expanding eval.
        if self.complex.is_completely_delooped() {
            self.elements.eval(h, t)
        } else {
            self.elements.eval_with(&self.complex, h, t)
        }
    }

    pub(crate) fn stat(&self) -> String {
        self.complex.stat()
    }

    /// Boundary-cutwidth profile of this builder's MinCut crossing order — a Cob-free pre-build
    /// dry-run (ties broken by index). The peak width predicts the dense-slice cost (~`2^peak`).
    pub(crate) fn profile(&self) -> BuildProfile {
        let nodes = self.nodes();
        let n = nodes.len();
        let mut remaining: Vec<usize> = (0..n).collect();
        let mut open: FxHashSet<Edge> = FxHashSet::default();

        let (order, widths): (Vec<usize>, Vec<usize>) = std::iter::from_fn(|| {
            let pos = remaining.iter()
                .position_min_by_key(|&&i| (cutwidth_after(&open, &[&nodes[i]]), i))?;
            let idx = remaining.swap_remove(pos);
            toggle_boundary(&mut open, &[&nodes[idx]]);
            Some((idx, open.len()))
        }).unzip();

        let peak = widths.iter().copied().max().unwrap_or(0);
        BuildProfile { n, order, widths, peak }
    }
}

/// Boundary-cutwidth profile of a crossing order: dense-slice cost peaks at ~`2^peak`.
pub(crate) struct BuildProfile {
    pub n: usize,
    #[allow(dead_code)] // replayed only by the faithfulness test
    pub order: Vec<usize>,  // node indices, in MinCut order
    pub widths: Vec<usize>, // boundary cutwidth after each step
    pub peak: usize,
}

impl fmt::Display for BuildProfile {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "n:    {}", self.n)?;
        writeln!(f, "peak: {}", self.peak)?;
        write!(f, "{}", sparkline(&self.widths, self.peak))
    }
}

/// Divide-and-conquer chunked build for a [`TngComplexBuilder`]: plan k chunks up front at the
/// deepest cutwidth valleys, build each via a child builder, and merge the reduced chunk into the
/// parent — so the parent never materializes the full dense slice.
struct ChunkBuilder<'a, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    builder: &'a TngComplexBuilder<R>,
}

impl<'a, R> ChunkBuilder<'a, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    // Plan k chunks up front and build each into a reduced sub-complex, paired with the crossings
    // it covers. Building is independent of the parent, so the caller merges them afterwards.
    fn build_chunks(&self) -> Vec<(Vec<Node>, (TngComplex<R>, Vec<TngComplexElem<R>>))> {
        let plan = self.plan();
        info!("{} chunk plan: {} pieces {:?}", self.builder.current_step(), plan.len(),
            plan.iter().map(|c| c.len()).collect_vec());

        plan.into_iter().map(|chunk| {
            let built = self.build_chunk(&chunk);
            (chunk, built)
        }).collect()
    }

    // Partition the crossings into pieces: `Auto(k)` at the deepest cutwidth valleys of the MinCut
    // order (contiguous, thin interface), or `Manual` by severing the given edge-cut(s).
    fn plan(&self) -> Vec<Vec<Node>> {
        let prof = self.builder.profile();
        let k = match &self.builder.config.cut {
            CutOption::Manual(cuts) => return self.manual_plan(cuts),
            // one node per unit here, so "cut after `c` crossings" = position `c - 1`.
            CutOption::AtCrossings(counts) => {
                let cuts = counts.iter()
                    .map(|&c| c.saturating_sub(1).min(prof.order.len().saturating_sub(2)))
                    .sorted().dedup().collect_vec();
                return self.segment_plan(&prof, &cuts);
            }
            CutOption::Auto(k) => (*k).max(1),
            CutOption::None => 1,
        };
        let cuts = select_cuts(&prof.widths, k - 1);
        self.segment_plan(&prof, &cuts)
    }

    // Segment `prof.order` after each cut position; map each index to its node.
    fn segment_plan(&self, prof: &BuildProfile, cuts: &[usize]) -> Vec<Vec<Node>> {
        let nodes = self.builder.nodes();
        let starts = std::iter::once(0).chain(cuts.iter().map(|&v| v + 1));
        let ends = cuts.iter().map(|&v| v + 1).chain(std::iter::once(prof.order.len()));
        starts.zip(ends)
            .map(|(s, e)| prof.order[s..e].iter().map(|&i| nodes[i].clone()).collect())
            .filter(|c: &Vec<Node>| !c.is_empty())
            .collect()
    }

    // Manual cut (no symmetry constraint): sever the cut edges (union), order the resulting pieces.
    fn manual_plan(&self, cuts: &[Vec<Edge>]) -> Vec<Vec<Node>> {
        let cut: FxHashSet<Edge> = cuts.iter().flatten().copied().collect();
        let nodes = self.builder.nodes();
        let pieces = cut_components(nodes, &cut);
        assert!(pieces.len() >= 2, "cut does not separate the link into ≥2 pieces");
        merge_order(nodes, pieces).into_iter()
            .map(|piece| piece.into_iter().map(|i| nodes[i].clone()).collect())
            .collect()
    }

    // Build `chunk` into a reduced sub-complex via a child builder, carrying its elements out too.
    fn build_chunk(&self, chunk: &[Node]) -> (TngComplex<R>, Vec<TngComplexElem<R>>) {
        let step = self.builder.current_step();
        let ends = boundary_edges(&chunk.iter().collect::<Vec<_>>()).into_iter().sorted().collect_vec();
        info!("{step} build chunk (n: {}, nb: {} {:?}): {}", chunk.len(), ends.len(), ends, chunk.iter().join(", "));

        let mut child = self.child_builder(chunk).run();
        let elems = child.elements_mut().take();
        let c = child.into_tng_complex();
        info!("{step} chunk built: {}", c.stat());
        (c, elems)
    }

    // A child builder over `chunk` (a sub-tangle), inheriting the parent's simplify mode;
    // chunking always uses the MinCut order, never recursing.
    fn child_builder(&self, chunk: &[Node]) -> TngComplexBuilder<R> {
        let (h, t) = self.builder.complex.ht();
        let base_pt = self.builder.complex.base_pt();
        let mut child = TngComplexBuilder::init(h, t, (0, 0), base_pt);
        child.set_nodes(chunk.iter().cloned());
        child.elements_mut().set(self.builder.elements().content().to_vec());

        // cap the child to the chunk's reachable band: a chunk vertex of weight
        // > b - deg_shift.0 can never reach the window (weight only grows).
        let h_range = self.builder.config.h_range.as_ref().map(|r| {
            let s = self.builder.complex.deg_shift().0;
            0 ..= (*r.end() - s).max(0)
        });
        let config = BuildConfig { mode: self.builder.config.mode, node_order: NodeOrder::MinCut, h_range, ..Default::default() };
        child.with_config(config)
    }
}

// Predicted raw size of the degree-`d` product slice, computable before building it:
// vertices `Σ_{d₁+d₂=d} |L[d₁]|·|R[d₂]|`, edges from the degree-`d-1` sources weighted by their
// out-edge counts (an upper bound: zero-reduced cobordisms and out-of-window targets drop out).
fn slice_estimate<R>(left: &TngComplex<R>, right: &TngComplex<R>, d: isize) -> (usize, usize)
where R: Ring, for<'x> &'x R: RingOps<R> {
    let out_edges = |c: &TngComplex<R>, i: isize| -> usize {
        c.keys_of_deg(i).map(|k| c.vertex(k).out_edges().count()).sum()
    };
    let verts = left.h_range().map(|d1| {
        left.rank(d1) * right.rank(d - d1)
    }).sum();
    let edges = left.h_range().map(|d1| {
        let d2 = d - 1 - d1;
        out_edges(left, d1) * right.rank(d2) + left.rank(d1) * out_edges(right, d2)
    }).sum();
    (verts, edges)
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    
    use super::*;

    // `profile`'s dry-run open-edge set must equal the real complex's `boundary_ends` at every step.
    #[test]
    fn dry_run_matches_real_boundary() {
        for name in ["6_2", "7_3", "8_19"] {
            let l = Link::test_data(name);
            let prof = TngComplexBuilder::<i32>::from_link(&l, &0, &0, false).profile();
            let nodes: Vec<Node> = l.nodes().cloned().collect();

            // raw merge (no deloop / eliminate) so we read the pure tangle boundary
            let mut b = TngComplexBuilder::<i32>::init(&0, &0, (0, 0), None)
                .with_config(BuildConfig { mode: BuildMode::None, ..Default::default() });

            let mut open: FxHashSet<Edge> = FxHashSet::default();
            for (step, &idx) in prof.order.iter().enumerate() {
                b.append_node(&nodes[idx]);
                toggle_boundary(&mut open, &[&nodes[idx]]);
                let real: FxHashSet<Edge> = b.complex().boundary_ends().collect();
                assert_eq!(open, real, "{name} step {step}: open-set vs boundary_ends");
                assert_eq!(prof.widths[step], open.len(), "{name} step {step}: width");
            }
            assert_eq!(*prof.widths.last().unwrap(), 0, "{name} should close up");
        }
    }

    #[test]
    fn test_unknot() {
        let l = Link::unknot();
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).run();
        let c = b.into_tng_complex().into_raw_complex();

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_unknot_rm1() {
        let l = Link::test_data("unknot_l_twist");
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).run();
        let c = b.into_tng_complex().into_raw_complex();

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_unknot_rm1_neg() {
        let l = Link::test_data("unknot_r_twist");
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).run();
        let c = b.into_tng_complex().into_raw_complex();

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
    }

    #[test]
    fn test_unknot_rm2() {
        let l = Link::test_data("unknot_lr_twist");
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).run();
        let c = b.into_tng_complex().into_raw_complex();

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
        assert_eq!(c[ 1].rank(), 0);
    }

    #[test]
    fn test_unlink_2() {
        // the R2 diagram: "unlink2" itself loads unoriented (its over-component has no under-anchor)
        // and Kh needs the orientation for its grading.
        let l = Link::test_data("unlink2_r2");
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).run();
        let c = b.into_tng_complex().into_raw_complex();

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 4);
        assert_eq!(c[ 1].rank(), 0);
    }

    #[test]
    fn test_tangle() { 
        let mut c = TngComplexBuilder::init(&0, &0, (0, 0), None);
        c.set_nodes([
            Node::from_pd_code([4,2,5,1]),
            Node::from_pd_code([3,6,4,1])
        ]);

        c.process_nodes();
        
        assert!(!c.complex.is_completely_delooped());
    }

    #[test]
    fn test_hopf_link() {
        let l = Link::test_data("L2a1");
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).run();
        let c = b.into_tng_complex().into_raw_complex();

        c.check_d_all();

        assert_eq!(c[-2].rank(), 2);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
    }

    #[test]
    fn test_build_modes_agree() {
        // all build modes must produce identical homology (incl. torsion).
        let l = Link::test_data("8_19");
        let build = |mode| {
            let config = BuildConfig { mode, ..Default::default() };
            TngComplexBuilder::from_link(&l, &0, &0, false).with_config(config).run()
                .into_tng_complex().into_raw_complex()
        };

        let ref_h = build(BuildMode::Greedy).homology();
        for mode in [BuildMode::MinFill, BuildMode::NoElim, BuildMode::None] {
            let c = build(mode);
            c.check_d_all();
            let h = c.homology();
            for i in 0..=8 {
                assert_eq!(h[i].rank(), ref_h[i].rank(), "rank at {i}, {mode:?}");
                assert_eq!(h[i].tors(), ref_h[i].tors(), "tors at {i}, {mode:?}");
            }
        }
    }

    #[test]
    fn test_no_full_deloop_agrees() {
        // no_full_deloop must not change homology: the deferred deloop is redone by into_raw_complex.
        let l = Link::test_data("8_19");
        let build = |skip| {
            let config = BuildConfig { no_full_deloop: skip, ..Default::default() };
            TngComplexBuilder::from_link(&l, &0, &0, false).with_config(config).run()
                .into_tng_complex().into_raw_complex()
        };

        let ref_h = build(false).homology();
        let c = build(true);
        c.check_d_all();
        let h = c.homology();
        for i in 0..=8 {
            assert_eq!(h[i].rank(), ref_h[i].rank(), "rank at {i}");
            assert_eq!(h[i].tors(), ref_h[i].tors(), "tors at {i}");
        }
    }

    #[test]
    fn test_chunk_build_matches() {
        // chunked builds must reproduce the non-chunked homology (incl. torsion).
        let l = Link::test_data("8_19");
        let build = |chunks: Option<usize>| {
            let config = BuildConfig { cut: chunks.map_or(CutOption::None, CutOption::Auto), ..Default::default() };
            TngComplexBuilder::from_link(&l, &0, &0, false).with_config(config).run()
                .into_tng_complex().into_raw_complex()
        };

        let ref_h = build(None).homology();
        for k in [Some(2), Some(3), Some(4)] {
            let c = build(k);
            c.check_d_all();
            let h = c.homology();
            for i in 0..=8 {
                assert_eq!(h[i].rank(), ref_h[i].rank(), "rank at {i}, chunks {k:?}");
                assert_eq!(h[i].tors(), ref_h[i].tors(), "tors at {i}, chunks {k:?}");
            }
        }
    }

    // Over Khovanov (d preserves q), a q-window keeps exactly the in-window generators, and the
    // bigraded homology at each kept (i, q) is unchanged. Runs both filter paths (greedy deloop and
    // no_full_deloop / matrix-level).
    #[test]
    fn q_filter_matches_full() {
        let l = Link::test_data("8_19");
        let full = KhComplex::new(&l, &0, &0, false);
        let (h_range, q_range) = (full.h_range(), full.q_range());
        let full_h = full.homology();

        // an interior window: drop the outermost occupied q on each side.
        let lo = *q_range.start() + 2;
        let hi = *q_range.end() - 2;

        for no_full_deloop in [false, true] {
            let config = BuildConfig { q_range: Some(lo..=hi), no_full_deloop, ..Default::default() };
            let win = KhComplex::new_with_config(&l, &0, &0, false, config);

            for i in win.h_range() {
                for x in win[i].raw_generators() {
                    let q = win.q_deg_of(x);
                    assert!((lo..=hi).contains(&q), "gen out of window: ({i}, {q}), no_full_deloop={no_full_deloop}");
                }
            }

            let win_h = win.homology();
            for i in h_range.clone() {
                for q in q_range.clone().step_by(2) {
                    let expected = if (lo..=hi).contains(&q) { full_h[(i, q)].rank() } else { 0 };
                    assert_eq!(win_h[(i, q)].rank(), expected, "rank ({i}, {q}), no_full_deloop={no_full_deloop}");
                    if (lo..=hi).contains(&q) {
                        assert_eq!(win_h[(i, q)].tors(), full_h[(i, q)].tors(), "tors ({i}, {q}), no_full_deloop={no_full_deloop}");
                    }
                }
            }
        }
    }

    // Over 𝔽₂[H] (cross-q edges via H), a window covering the whole complex reproduces the
    // unfiltered homology exactly — checks the edge-skipping path is a faithful no-op. Compared
    // singly-graded (per h), since `deg H = −2` makes the bigraded split ill-defined here.
    #[test]
    fn q_filter_full_window_identity() {
        use yui_core::poly::Poly;
        use yui_core::num::FF2;
        type P = Poly<'H', FF2>;

        let l = Link::test_data("6_2");
        let (h, t) = (P::variable(), P::zero());
        let full = KhComplex::new(&l, &h, &t, false);
        let q = full.q_range();
        let wide = (*q.start() - 4) ..= (*q.end() + 4);

        let config = BuildConfig { q_range: Some(wide), ..Default::default() };
        let win = KhComplex::new_with_config(&l, &h, &t, false, config);

        let (fh, wh) = (full.homology(), win.homology());
        for i in full.h_range() {
            assert_eq!(wh[i].rank(), fh[i].rank(), "rank at {i}");
            assert_eq!(wh[i].tors(), fh[i].tors(), "tors at {i}");
        }
    }

    // Over 𝔽₂[H], `d` raises generator-q (via H), so an upper-unbounded window `{q ≥ lo}` is a
    // genuine subcomplex: a proper truncation that must still satisfy `d² = 0`.
    #[test]
    fn q_filter_subcomplex_valid() {
        use yui_core::poly::Poly;
        use yui_core::num::FF2;
        type P = Poly<'H', FF2>;

        let l = Link::test_data("6_2");
        let (h, t) = (P::variable(), P::zero());
        let full = KhComplex::new(&l, &h, &t, false);
        let lo = *full.q_range().start() + 2; // drop the bottom q-degree(s)

        let config = BuildConfig { q_range: Some(lo ..= isize::MAX), ..Default::default() };
        let win = KhComplex::new_with_config(&l, &h, &t, false, config);

        win.inner().check_d_all(); // `{q ≥ lo}` is a subcomplex: valid d² = 0
        for i in win.h_range() {
            for x in win[i].raw_generators() {
                assert!(win.q_deg_of(x) >= lo, "gen below window at ({i}, {})", win.q_deg_of(x));
            }
        }

        let gens = |c: &KhComplex<P>| c.h_range().map(|i| c[i].rank()).sum::<usize>();
        assert!(gens(&win) < gens(&full), "window dropped no generator");
    }

    // chunked builds must track the canon cycles too: the Lee-class divisibility (the ss
    // ingredient) is computed from each chunked homology and must match the non-chunked one.
    #[test]
    fn test_chunk_elements_match() {
        use crate::kh::KhHomology;
        use crate::util::calc::div_vec;

        let l = Link::test_data("8_19");
        let c = 2;
        let div = |chunks: Option<usize>| {
            let config = BuildConfig { cut: chunks.map_or(CutOption::None, CutOption::Auto), ..Default::default() };
            let kh = KhHomology::new_with_config(&l, &c, &0, false, config);
            kh.canon_cycles().iter()
                .map(|z| div_vec(&kh[0].vectorize_euc(z).subvec(0..2), &c).unwrap())
                .collect_vec()
        };

        let ref_d = div(None);
        for k in [Some(2), Some(3), Some(4)] {
            assert_eq!(div(k), ref_d, "divisibility, chunks {k:?}");
        }
    }

    #[test]
    fn test_8_19() {
        let l = Link::test_data("8_19");
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).run();
        let c = b.into_tng_complex().into_raw_complex();

        c.check_d_all();

        let h = c.homology();

        for i in [1,6,7,8] {
            assert_eq!(h[i].rank(), 0);
            assert!(h[i].is_free());
        }

        for i in [0,4,5] {
            assert_eq!(h[i].rank(), 2);
            assert!(h[i].is_free());
        }

        assert_eq!(h[2].rank(), 1);
        assert!(h[2].is_free());

        assert_eq!(h[3].rank(), 1);
        assert_eq!(h[3].tors(), &vec![2]);
    }

    #[test]
    fn test_8_19_h_range() {
        let l = Link::test_data("8_19");
        let config = BuildConfig { h_range: Some(2..=6), ..Default::default() };
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).with_config(config).run();
        let c = b.into_tng_complex().into_raw_complex();

        // literal truncation: chain groups vanish outside [2, 6].
        for i in [0, 1, 7, 8] {
            assert_eq!(c[i].rank(), 0);
        }

        c.check_d_all();

        let h = c.homology();

        // interior degrees are correct (the endpoints 2 and 6 are not).
        assert_eq!(h[3].rank(), 1);
        assert_eq!(h[3].tors(), &vec![2]);

        assert_eq!(h[4].rank(), 2);
        assert!(h[4].is_free());

        assert_eq!(h[5].rank(), 2);
        assert!(h[5].is_free());
    }

    #[test]
    fn test_8_19_h_range_full() {
        // A range covering the whole complex must reproduce the full homology.
        let l = Link::test_data("8_19");
        let config = BuildConfig { h_range: Some(0..=8), ..Default::default() };
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).with_config(config).run();
        let c = b.into_tng_complex().into_raw_complex();

        c.check_d_all();

        let h = c.homology();

        for i in [1, 6, 7, 8] {
            assert_eq!(h[i].rank(), 0);
        }
        for i in [0, 4, 5] {
            assert_eq!(h[i].rank(), 2);
            assert!(h[i].is_free());
        }
        assert_eq!(h[2].rank(), 1);
        assert!(h[2].is_free());
        assert_eq!(h[3].rank(), 1);
        assert_eq!(h[3].tors(), &vec![2]);
    }

    #[test]
    fn test_8_19_h_range_empty() {
        // A range disjoint from the complex's degrees yields an empty complex.
        let l = Link::test_data("8_19");
        let config = BuildConfig { h_range: Some(20..=20), ..Default::default() };
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).with_config(config).run();
        let c = b.into_tng_complex().into_raw_complex();

        for i in 0..=8 {
            assert_eq!(c[i].rank(), 0);
        }
    }

    #[test]
    fn canon_cycle_trefoil() {
        let l = Link::test_data("3_1");
        let b = TngComplexBuilder::from_link(&l, &1, &0, false).run();
        let zs = b.eval_elements();
        let c = b.into_tng_complex().into_raw_complex();

        assert_eq!(zs.len(), 2);
        assert_ne!(zs[0], zs[1]);
        
        for z in zs {
            assert!(c.d(0, &z).is_zero());
        }
    }
}