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

use rustc_hash::{FxHashSet, FxHashMap};
use itertools::Itertools;
use log::{debug, info, trace};
use yui_core::{Ring, RingOps};
use yui_link::{Node, Edge, Link};

use crate::kh::{KhChain, KhComplex};
use crate::tng::{End, TngComplexElem, LcCobTrait, TngComplex, TngComplexKey};
use super::{reachable_range, pop_min_pivot, sparkline, cutwidth_after, toggle_boundary, TngElemBuilder};

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
    Selective, // deloop only productive circles, eliminate immediately (full deloop at merge end)
    MinFill,   // deloop a whole degree, then eliminate by global min-fill (Markowitz)
    NoElim,    // deloop every circle but don't eliminate (delooped, unreduced complex)
    None,      // don't deloop, don't eliminate (raw merge; finalize still deloops to a valid complex)
}

impl BuildMode {
    // whether any simplification happens during the build (None = raw merge only).
    pub fn is_active(&self) -> bool {
        *self != BuildMode::None
    }

    pub fn is_selective(&self) -> bool {
        *self == BuildMode::Selective
    }

    pub fn is_min_fill(&self) -> bool {
        *self == BuildMode::MinFill
    }

    // whether each newly-delooped vertex is eliminated inline (vs deferred / not at all).
    pub fn immediate_elim(&self) -> bool {
        matches!(self, BuildMode::Greedy | BuildMode::Selective)
    }
}

/// Toggles for the automatic simplification done while building.
#[derive(Clone, Debug)]
pub struct BuildConfig {
    pub node_order: NodeOrder,
    pub mode: BuildMode,
    pub h_range: Option<RangeInclusive<isize>>,
}

impl Default for BuildConfig {
    fn default() -> Self {
        Self { node_order: NodeOrder::default(), mode: BuildMode::default(), h_range: None }
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
        let base_pt = if reduced { l.base_pt() } else { None };
        let deg_shift = KhComplex::deg_shift_for(l, reduced);

        let mut b = Self::init(h, t, deg_shift, base_pt);
        b.set_nodes(l.nodes().cloned());
        b.set_loops(l.loops().iter().cloned());

        if t.is_zero() && l.is_knot() {
            let canon = TngComplexElem::canon_cycles(l, base_pt);
            b.set_elements(canon);
        }

        b
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
        // canon cycles live in h-degree 0; drop them if the range excludes it.
        if let Some(range) = &config.h_range {
            if !range.contains(&0) {
                self.elements.clear();
            }
        }
        self.config = config;
        self
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

    pub fn set_elements<I>(&mut self, elements: I)
    where I: IntoIterator<Item = TngComplexElem<R>> {
        self.elements.set(elements);
    }

    pub(crate) fn take_elements(&mut self) -> Vec<TngComplexElem<R>> {
        self.elements.take()
    }

    pub fn run(mut self) -> Self {
        info!("build config:\n{:#?}", self.config);
        info!("cutwidth profile:\n{}", self.profile());
        self.process_nodes();
        self.process_free_loops();
        self.finalize();
        self
    }

    // See [BN07, §7] (scan-and-cancel algorithm).
    pub(crate) fn process_nodes(&mut self) {
        info!("{} process {} nodes", self.current_step(), self.n_nodes());

        while let Some(x) = self.choose_next_node().cloned() {
            self.append_node(&x)
        }
    }

    /// Pick the next node_order: maximize `score_node`, ties broken by earliest crossing order
    /// (`self.nodes` keeps PD order — `prepare_append` removes via order-preserving `Vec::remove`).
    pub(crate) fn choose_next_node(&self) -> Option<&Node> {
        self.nodes.iter().enumerate()
            .min_by_key(|(i, x)| (-self.score_node(x), *i))
            .map(|(_, x)| x)
    }

    /// Strategy score for appending `x` (higher is better).
    pub(crate) fn score_node(&self, x: &Node) -> isize {
        match self.config.node_order {
            NodeOrder::MinCut => -self.cutwidth(x),
            NodeOrder::Given => 0, // constant → ties broken by earliest index = given order
        }
    }

    /// Boundary cutwidth (open-edge count) after appending `x`. `boundary_ends` is cheap, so
    /// recomputing it per call is fine.
    pub(crate) fn cutwidth(&self, x: &Node) -> isize {
        self.cutwidth_of(x.edges().iter().copied())
    }

    /// Cutwidth after toggling an arbitrary edge multiset — for τ-pairs, where `x` and `τx`
    /// must be scored together (a shared axis edge toggles twice and cancels).
    pub(crate) fn cutwidth_of(&self, edges: impl IntoIterator<Item = Edge>) -> isize {
        let boundary: FxHashSet<Edge> = self.complex.boundary_ends().collect();
        let mut cnt: FxHashMap<Edge, u32> = FxHashMap::default();
        for e in edges {
            *cnt.entry(e).or_default() += 1;
        }
        let delta: isize = cnt.iter()
            .filter(|(_, c)| *c % 2 == 1)
            .map(|(e, _)| if boundary.contains(e) { -1 } else { 1 })
            .sum();
        boundary.len() as isize + delta
    }

    // "(committed/total)" crossing progress, for log prefixes.
    pub(crate) fn current_step(&self) -> String {
        format!("({}/{})", self.complex.dim(), self.complex.dim() + self.n_nodes())
    }

    pub(crate) fn append_node(&mut self, x: &Node) {
        info!("{} append: {x}", self.current_step());

        self.prepare_append(x);
        
        let (h, t) = self.complex.ht();
        let cx = TngComplex::from_node(h, t, x, self.complex.base_pt());
        self.merge(cx);
    }

    pub(crate) fn prepare_append(&mut self, x: &Node) { 
        if let Some(i) = self.nodes.iter().find_position(|&e| e == x) { 
            self.nodes.remove(i.0);
        }

        self.elements.append_node(x);
    }

    pub(crate) fn merge(&mut self, other: TngComplex<R>) {
        let (left, right) = self.complex.prepare_merge(other);
        let range = reachable_range(self.complex.h_range(), &self.config.h_range, self.n_nodes());

        debug!("{} merge {} <- {}", self.current_step(), left.stat(), right.stat());
        debug!("  merge range: {:?}", range);

        match self.config.mode {
            BuildMode::None    => self.complex.merge_with(&left, &right),
            BuildMode::MinFill => self.merge_deferred(&left, &right, range),
            _                  => self.merge_default(&left, &right, range),
        }

        self.prune_h_range();
        debug!("{} merged: {}", self.current_step(), self.stat());
    }

    // Default path (Greedy / Selective / NoElim): per degree, eliminate (if the mode does) then deloop.
    // Deloop reads `selective` from the mode (defer non-productive circles, then full-deloop + re-pass at the end).
    fn merge_default(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, range: RangeInclusive<isize>) {
        let top = *range.end();
        let selective = self.config.mode.is_selective();

        for i in range {
            debug!("{} build C[{i}]...", self.current_step());
            self.merge_slice(left, right, i);
            if self.config.mode.immediate_elim() {
                self.eliminate_in(i - 1);
            }
            self.deloop_in(i - 1);
            debug!("{} built C[{i}]: {}", self.current_step(), self.complex.rank(i));
        }

        self.prune_isolated_top(top);
        self.deloop_in(top);

        // re-run selective to a fixpoint (catch loops turned productive by later elims), then full-deloop the rest.
        if selective {
            self.deloop_all_selective();
            self.deloop_all_forced();
        }
    }

    // Repeatedly deloop productive circles until none remain — each elimination can turn a
    // previously-deferred loop productive.
    fn deloop_all_selective(&mut self) {
        for step in 1.. {
            let before = self.complex.n_verts();
            debug!("selective re-pass {step}: start ({before} verts)");
            self.deloop_all();
            let after = self.complex.n_verts();
            debug!("  selective re-pass {step}: {before} -> {after} verts (diff {})",
                after as isize - before as isize);
            if after == before { break }
        }
    }

    // MinFill: per degree, deloop then eliminate i-1, i by global min-fill (incremental Markowitz).
    // Always full-deloops — combining it with selective delooping only balloons the transient complex.
    fn merge_deferred(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, range: RangeInclusive<isize>) {
        let top = *range.end();

        for i in range {
            debug!("{} build C[{i}]...", self.current_step());
            self.merge_slice(left, right, i);
            self.deloop_in(i - 1);
            self.eliminate_in(i - 2);
            self.eliminate_in(i - 1);
            debug!("{} built C[{i}]: {}", self.current_step(), self.complex.rank(i));
        }

        self.prune_isolated_top(top);
        self.deloop_in(top);
        self.eliminate_in(top - 1);
        // no eliminate_in(top): top has no outgoing edges, so it would be a no-op.
    }

    // Build degree `i`: merge in its vertices and the edges into it.
    pub(super) fn merge_slice(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, i: isize) {
        let nv = self.complex.merge_vertices(left, right, i);
        debug!("  +{nv} verts");
        let ne = self.complex.merge_edges(left, right, i - 1);
        debug!("  +{ne} edges");
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

    /// Keys at degree `i` matching `pred`, each paired with `weight(k)`; callers pivot by least
    /// weight via [`pop_min_pivot`].
    pub(crate) fn collect_keys<F, W>(&self, i: isize, pred: F, weight: W) -> Vec<(TngComplexKey, usize)>
    where F: Fn(&TngComplexKey) -> bool, W: Fn(&TngComplexKey) -> usize {
        self.complex.keys_of_deg(i)
            .filter(|k| pred(k))
            .map(|k| (*k, weight(k)))
            .collect_vec()
    }

    // Selective mode deloops only "productive" circles — those whose no-dot cap yields an
    // invertible edge (deloop+elim fires, no doubling); otherwise any circle.
    pub(crate) fn find_loop_in(&self, k: &TngComplexKey, allow_based: bool, selective: bool) -> Option<usize> {
        let v = self.complex.vertex(k);
        v.tng().comps().enumerate()
            .filter(|(_, c)| c.is_circle() && (allow_based || !c.is_marked()))
            .find(|(_, c)| !selective
             || v.out_edges().any(|l| self.complex.edge(k, l).is_invertible_after_cap(End::Src, c))
             || v.in_edges().any(|j| self.complex.edge(j, k).is_invertible_after_cap(End::Tgt, c)))
            .map(|(r, _)| r)
    }

    // Deloop unmarked loops over all degrees, selective per `config.mode`.
    pub(crate) fn deloop_all(&mut self) {
        for i in self.complex.h_range() {
            self.deloop_in(i);
        }
    }

    // Deloop all unmarked loops, ignoring the `selective` config — for the final cleanup.
    fn deloop_all_forced(&mut self) {
        for i in self.complex.h_range() {
            self.deloop_in_with(i, false, false);
        }
    }

    // Deloop the marked (based) loops into a single summand. All unmarked loops must already
    // be gone — else delooping only the marked ones would break the complex.
    fn deloop_all_marked(&mut self) {
        debug_assert!(
            self.complex.keys().all(|k| self.find_loop_in(k, false, false).is_none()),
            "deloop_all_marked: unmarked loops remain"
        );
        for i in self.complex.h_range() {
            self.deloop_in_with(i, true, false);
        }
    }

    fn deloop_in(&mut self, i: isize) {
        let selective = self.config.mode.is_selective();
        self.deloop_in_with(i, false, selective);
    }

    fn deloop_in_with(&mut self, i: isize, allow_based: bool, selective: bool) {
        let mut keys = self.collect_keys(i,
            |k| self.find_loop_in(k, allow_based, selective).is_some(),
            |k| self.complex.vertex(k).c_weight(),
        );
        if keys.is_empty() { return }

        debug!("{} deloop in C[{i}], targets: {}{}.", self.current_step(), keys.len(), if selective { " (selective)" } else { "" });

        let before = self.complex.rank(i) as isize;

        while let Some(k) = pop_min_pivot(&mut keys, |k|
            self.complex.contains_key(k).then(|| self.complex.vertex(k).c_weight())
        ) {
            let Some(r) = self.find_loop_in(&k, allow_based, selective) else { continue };

            for new_key in self.deloop(&k, r) {
                if self.find_loop_in(&new_key, allow_based, selective).is_some() {
                    let w = self.complex.vertex(&new_key).c_weight();
                    keys.push((new_key, w));
                }
            }
        }

        let after = self.complex.rank(i) as isize;

        debug!("{}   delooped C[{i}]: {} (diff: {}).", self.current_step(), after, after - before);
    }

    pub(crate) fn deloop(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> {
        let c = self.complex.vertex(k).tng().comp(r);

        trace!("{} deloop {c} in {}", self.stat(), self.complex.vertex(k));

        self.elements.deloop(k, c);

        let mut added = self.complex.deloop(k, r);

        // immediate elim eliminates each new vertex now; min-fill leaves them for the post-deloop
        // global pass, None leaves them entirely.
        if self.config.mode.immediate_elim() {
            // retain only the keys that weren't eliminated
            added.retain(|k| !self.try_eliminate_at(k));
        }
        added
    }

    pub(crate) fn eliminate_in(&mut self, i: isize) {
        let mut keys = self.collect_keys(i,
            |k| self.complex.vertex(k).out_edges().any(|l|
                self.complex.edge(k, l).is_invertible()
            ),
            |k| self.complex.elim_cost(k),
        );
        if keys.is_empty() { return }

        debug!("{} eliminate in C[{i}], targets: {}", self.current_step(), keys.len());

        let before = self.complex.rank(i) as isize;

        while let Some(k) = pop_min_pivot(&mut keys, |k|
            self.complex.contains_key(k).then(|| self.complex.elim_cost(k))
        ) {
            self.try_eliminate_at(&k);
        }

        let after = self.complex.rank(i) as isize;

        debug!("{}   eliminated C[{i}]: {} (diff: {}).", self.current_step(), after, after - before);
    }

    pub(crate) fn try_eliminate_at(&mut self, k: &TngComplexKey) -> bool {
        if let Some(&j) = self.choose_inv_edge_into(&k) { 
            self.eliminate(&j, &k);
            true
        } else if let Some(&l) = self.choose_inv_edge_from(&k) { 
            self.eliminate(&k, &l);
            true
        } else { 
            false
        }
    }

    fn choose_inv_edge_into(&self, k: &TngComplexKey) -> Option<&TngComplexKey> { 
        self.complex.vertex(k).in_edges().filter_map(|j|
            self.complex.edge(j, k).is_invertible().then_some(j)
        )
        .min_by_key(|j| (self.complex.edge_weight(j, k), **j))
    }

    fn choose_inv_edge_from(&self, k: &TngComplexKey) -> Option<&TngComplexKey> {
        self.complex.vertex(k).out_edges().filter_map(|l|
            self.complex.edge(k, l).is_invertible().then_some(l)
        )
        .min_by_key(|l| (self.complex.edge_weight(k, l), **l))
    }

    pub(crate) fn eliminate(&mut self, i: &TngComplexKey, j: &TngComplexKey) {
        trace!("{} eliminate {}: {} -> {}", self.stat(), self.complex.edge(i, j), self.complex.vertex(i), self.complex.vertex(j));
        
        self.elements.eliminate(&self.complex, i, j);
        self.complex.eliminate(i, j);
    }

    pub(crate) fn process_free_loops(&mut self) {
        while !self.loops.is_empty() { 
            let c = self.loops.remove(0);

            self.elements.insert_loop(c);

            let (h, t) = self.complex.ht();
            let marked = self.complex.base_pt() == Some(c);
            let c = TngComplex::from_loop(h, t, c, marked);
            self.merge(c);

            if self.config.mode.is_active() {
                self.deloop_all_forced();
            }
        }
    }

    fn finalize(&mut self) {
        if self.complex.is_completely_delooped() {
            info!("{} completely delooped: {}", self.current_step(), self.stat());
            return;
        }

        info!("{} finalize: {}", self.current_step(), self.stat());

        self.deloop_all_forced();
        self.deloop_all_marked(); // deloop marked loops

        info!("{} finalized: {}", self.current_step(), self.stat());
    }

    pub fn into_tng_complex(self) -> TngComplex<R> {
        self.complex
    }

    pub fn eval_elements(&self) -> Vec<KhChain<R>> {
        let (h, t) = self.complex.ht();
        self.elements.eval(h, t)
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
        let l = Link::test_data("unlink2");
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
        for mode in [BuildMode::Selective, BuildMode::MinFill, BuildMode::NoElim, BuildMode::None] {
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