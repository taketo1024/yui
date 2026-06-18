//! Up-front chunk planning for [`SymTngBuilder`]: partition the crossings into `k` τ-closed pieces
//! (or at user-specified cuts), ordered so each merges behind a thin interface. This module only
//! decides the partition — the parent builds each piece and merges it back.

use itertools::Itertools;
use log::info;
use rustc_hash::{FxHashMap, FxHashSet};
use yui_core::{Ring, RingOps};
use yui_link::{Node, Edge};

use super::SymTngBuilder;
use super::boundary_edges;

/// How to partition the crossings into chunks.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub enum ChunkStrategy {
    /// Cut the MinCut crossing order at its `k-1` deepest cutwidth valleys (thin build frontier).
    #[default]
    Frontier,
    /// Recursive min-boundary bisection of the τ-unit graph into `k` balanced pieces (thin interface).
    Boundary,
    /// User-specified cut edges (each list a τ-symmetric cut); the pieces are the severed components.
    Manual(Vec<Vec<Edge>>),
}

pub(crate) struct ChunkBuilder<'a, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    builder: &'a SymTngBuilder<R>,
}

impl<'a, R> ChunkBuilder<'a, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub(crate) fn new(builder: &'a SymTngBuilder<R>) -> Self {
        Self { builder }
    }

    // Plan the partition into τ-closed pieces (node lists), ordered so each merges behind a thin interface.
    pub(crate) fn plan(&self) -> Vec<Vec<Node>> {
        let nodes = self.builder.nodes();
        let cfg = self.builder.config();

        let pieces = match &cfg.chunk_strategy {
            ChunkStrategy::Frontier => self.frontier(cfg.chunks.expect("`chunks` required for Frontier")),
            ChunkStrategy::Boundary => self.boundary(cfg.chunks.expect("`chunks` required for Boundary")),
            ChunkStrategy::Manual(cuts) => self.manual(cuts),
        };

        let ordered = self.merge_order(pieces);
        self.log_plan(&ordered);

        ordered.into_iter()
            .map(|piece| piece.iter().map(|&i| nodes[i].clone()).collect())
            .collect()
    }

    // τ-units: on-axis crossings as singletons, off-axis crossings paired at the lower index.
    fn units(&self) -> Vec<Vec<usize>> {
        let nodes = self.builder.nodes();
        let idx_of: FxHashMap<Node, usize> = nodes.iter().enumerate().map(|(i, x)| (x.clone(), i)).collect();
        (0..nodes.len()).filter_map(|i| {
            let j = idx_of[self.builder.inv_node(&nodes[i])];
            match j {
                _ if j == i => Some(vec![i]),
                _ if i < j  => Some(vec![i, j]),
                _           => None,
            }
        }).collect()
    }

    // Frontier: cut `profile_sym`'s MinCut order at its k-1 deepest cutwidth valleys → k segments.
    fn frontier(&self, k: usize) -> Vec<Vec<usize>> {
        let prof = self.builder.profile_sym();
        let cuts = deepest_valleys(&prof.widths, k - 1);
        if cuts.len() < k - 1 {
            info!("frontier: only {} valleys for {k} chunks → {} pieces", cuts.len(), cuts.len() + 1);
        }
        let mut pieces = vec![];
        let mut start = 0;
        for end in cuts.iter().copied().chain([prof.order.len() - 1]) {
            pieces.push(prof.order[start..=end].iter().flatten().copied().collect());
            start = end + 1;
        }
        pieces
    }

    // Boundary: recursive min-boundary bisection of the τ-unit graph into k balanced pieces.
    fn boundary(&self, k: usize) -> Vec<Vec<usize>> {
        let units = self.units();
        let adj = self.unit_adjacency(&units);
        let sizes: Vec<usize> = units.iter().map(|u| u.len()).collect();
        partition_k(&adj, &sizes, k).into_iter()
            .map(|piece| piece.into_iter().flat_map(|u| units[u].clone()).collect())
            .collect()
    }

    // Manual: validate the cut edges, then return the components left after severing them.
    fn manual(&self, cuts: &[Vec<Edge>]) -> Vec<Vec<usize>> {
        let cut: FxHashSet<Edge> = cuts.iter().flatten().copied().collect();
        self.validate_cut(&cut);
        self.cut_components(&cut)
    }

    // A manual cut must be τ-symmetric (closed under `inv_edge`), separate the link into ≥2 pieces,
    // and leave each piece τ-invariant (so the sym build can pair `x` with `τx` inside it).
    fn validate_cut(&self, cut: &FxHashSet<Edge>) {
        for &e in cut {
            assert!(cut.contains(&self.builder.inv_edge(e)), "cut not τ-symmetric: τ-image of edge {e} missing");
        }
        let comps = self.cut_components(cut);
        assert!(comps.len() >= 2, "cut does not separate the link into ≥2 pieces");

        let nodes = self.builder.nodes();
        let idx_of: FxHashMap<Node, usize> = nodes.iter().enumerate().map(|(i, x)| (x.clone(), i)).collect();
        for comp in &comps {
            let set: FxHashSet<usize> = comp.iter().copied().collect();
            let tau_in = comp.iter().all(|&i| set.contains(&idx_of[self.builder.inv_node(&nodes[i])]));
            assert!(tau_in, "a cut piece is not τ-invariant (τ maps it outside)");
        }
    }

    // Crossing components after severing the `cut` edges (BFS over non-cut shared edges).
    fn cut_components(&self, cut: &FxHashSet<Edge>) -> Vec<Vec<usize>> {
        let nodes = self.builder.nodes();
        let mut by_edge: FxHashMap<Edge, Vec<usize>> = FxHashMap::default();
        for (i, x) in nodes.iter().enumerate() {
            for &e in x.edges() {
                if !cut.contains(&e) { by_edge.entry(e).or_default().push(i); }
            }
        }
        let mut seen = vec![false; nodes.len()];
        let mut comps = vec![];
        for start in 0..nodes.len() {
            if seen[start] { continue }
            seen[start] = true;
            let mut stack = vec![start];
            let mut comp = vec![];
            while let Some(i) = stack.pop() {
                comp.push(i);
                for &e in nodes[i].edges() {
                    if cut.contains(&e) { continue }
                    for &j in by_edge.get(&e).into_iter().flatten() {
                        if !seen[j] { seen[j] = true; stack.push(j); }
                    }
                }
            }
            comps.push(comp);
        }
        comps
    }

    // Unit-adjacency graph: each knot-edge straddling two units is one strand crossing their boundary.
    fn unit_adjacency(&self, units: &[Vec<usize>]) -> Vec<Vec<(usize, usize)>> {
        let nodes = self.builder.nodes();
        let mut unit_of = vec![0usize; nodes.len()];
        for (u, members) in units.iter().enumerate() {
            for &i in members { unit_of[i] = u; }
        }
        let mut edge_units: FxHashMap<Edge, Vec<usize>> = FxHashMap::default();
        for (i, x) in nodes.iter().enumerate() {
            for &e in x.edges() { edge_units.entry(e).or_default().push(unit_of[i]); }
        }
        let mut adj: Vec<FxHashMap<usize, usize>> = vec![FxHashMap::default(); units.len()];
        for us in edge_units.values() {
            if let [u, v] = us[..] {
                if u != v {
                    *adj[u].entry(v).or_default() += 1;
                    *adj[v].entry(u).or_default() += 1;
                }
            }
        }
        adj.into_iter().map(|m| m.into_iter().collect()).collect()
    }

    // Boundary edges of a node-index piece — its merge interface (edges with one endpoint inside).
    fn chunk_ends(&self, piece: &[usize]) -> Vec<Edge> {
        let nodes = self.builder.nodes();
        let subset: Vec<&Node> = piece.iter().map(|&i| &nodes[i]).collect();
        boundary_edges(&subset).into_iter().sorted().collect()
    }

    // Order pieces to minimize the peak merge interface (exact for few pieces, greedy beyond).
    fn merge_order(&self, plan: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
        let k = plan.len();
        if k <= 2 {
            return plan;
        }
        let ends: Vec<FxHashSet<Edge>> = plan.iter().map(|c| self.chunk_ends(c).into_iter().collect()).collect();
        let order = if k <= 8 {
            (0..k).permutations(k).min_by_key(|o| interface_profile(o, &ends)).unwrap()
        } else {
            greedy_merge_order(&ends)
        };
        order.into_iter().map(|i| plan[i].clone()).collect()
    }

    // Sequential-merge interfaces (plan order): shared boundary edges between the merged-so-far and the next.
    fn merge_interfaces(&self, plan: &[Vec<usize>]) -> Vec<usize> {
        let mut acc: Vec<usize> = vec![];
        let mut out = vec![];
        for chunk in plan {
            if !acc.is_empty() {
                let acc_b: FxHashSet<Edge> = self.chunk_ends(&acc).into_iter().collect();
                out.push(self.chunk_ends(chunk).iter().filter(|e| acc_b.contains(e)).count());
            }
            acc.extend_from_slice(chunk);
        }
        out
    }

    fn log_plan(&self, plan: &[Vec<usize>]) {
        let interfaces = self.merge_interfaces(plan);
        let peak = interfaces.iter().max().copied().unwrap_or(0);
        let sizes = plan.iter().map(|c| c.len()).collect_vec();
        info!("chunk plan: {} pieces {:?}, merge interfaces {:?} (peak {peak})", plan.len(), sizes, interfaces);
    }
}

// The k deepest cutwidth-valley positions (a descent that turns back up; flats ignored), as indices.
fn deepest_valleys(widths: &[usize], k: usize) -> Vec<usize> {
    let mut valleys = vec![];
    let mut descended = false;
    for i in 1..widths.len() {
        if widths[i] < widths[i - 1] {
            descended = true;
        } else if widths[i] > widths[i - 1] && descended {
            valleys.push(i - 1);
            descended = false;
        }
    }
    valleys.sort_by_key(|&p| widths[p]); // deepest (smallest width) first
    valleys.truncate(k);
    valleys.sort();
    valleys
}

// Recursive min-boundary bisection into k pieces: repeatedly split the largest splittable piece.
fn partition_k(adj: &[Vec<(usize, usize)>], sizes: &[usize], k: usize) -> Vec<Vec<usize>> {
    let mut pieces: Vec<Vec<usize>> = vec![(0..adj.len()).collect()];
    while pieces.len() < k {
        let Some(idx) = pieces.iter().enumerate()
            .filter(|(_, p)| p.len() > 1)
            .max_by_key(|(_, p)| p.iter().map(|&u| sizes[u]).sum::<usize>())
            .map(|(i, _)| i)
        else { break };
        let piece = pieces.swap_remove(idx);
        let (a, b) = min_boundary_bisection(&piece, adj, sizes);
        pieces.push(a);
        pieces.push(b);
    }
    pieces
}

// Min-boundary bisection (Fiduccia–Mattheyses): thinnest edge-cut, ties to the most balanced split;
// each side kept ≥ ⌊total/3⌋ crossings — the floor that bars the degenerate peel-one-crossing cut.
fn min_boundary_bisection(piece: &[usize], adj: &[Vec<(usize, usize)>], sizes: &[usize]) -> (Vec<usize>, Vec<usize>) {
    let n = piece.len();
    let pos: FxHashMap<usize, usize> = piece.iter().enumerate().map(|(p, &u)| (u, p)).collect();
    let nbr: Vec<Vec<(usize, usize)>> = piece.iter().map(|&u| {
        adj[u].iter().filter_map(|&(v, w)| pos.get(&v).map(|&q| (q, w))).collect()
    }).collect();
    let wt: Vec<usize> = piece.iter().map(|&u| sizes[u]).collect();
    let min_side = (wt.iter().sum::<usize>() / 3).max(1);

    let mass = |s: &[bool], v: bool| (0..n).filter(|&p| s[p] == v).map(|p| wt[p]).sum::<usize>();
    let cut_of = |s: &[bool]| -> usize {
        (0..n).flat_map(|p| nbr[p].iter().map(move |&(q, w)| (p, q, w)))
            .filter(|&(p, q, _)| p < q && s[p] != s[q]).map(|(_, _, w)| w).sum()
    };
    // width first, then most balanced (smallest crossing-mass gap).
    let score = |s: &[bool]| -> (usize, usize) { (cut_of(s), mass(s, true).abs_diff(mass(s, false))) };

    let mut side: Vec<bool> = (0..n).map(|p| p >= n / 2).collect();

    loop {
        let mut cur = side.clone();
        let mut locked = vec![false; n];
        let (mut best, mut best_score) = (cur.clone(), score(&cur));
        for _ in 0..n {
            let pick = (0..n)
                .filter(|&p| !locked[p] && mass(&cur, cur[p]) - wt[p] >= min_side)
                .max_by_key(|&p| nbr[p].iter()
                    .map(|&(q, w)| if cur[q] != cur[p] { w as isize } else { -(w as isize) })
                    .sum::<isize>());
            let Some(p) = pick else { break };
            cur[p] = !cur[p];
            locked[p] = true;
            let s = score(&cur);
            if s < best_score { (best_score, best) = (s, cur.clone()); }
        }
        if best_score < score(&side) { side = best; } else { break; }
    }

    let (mut a, mut b) = (vec![], vec![]);
    for (p, &u) in piece.iter().enumerate() {
        if side[p] { a.push(u); } else { b.push(u); }
    }
    (a, b)
}

// Interfaces of a merge order, sorted worst-first (0-share → ∞ penalty) — the key we minimize.
fn interface_profile(order: &[usize], ends: &[FxHashSet<Edge>]) -> Vec<usize> {
    let mut acc = ends[order[0]].clone();
    let mut ifs: Vec<usize> = order[1..].iter().map(|&i| {
        let shared = ends[i].iter().filter(|e| acc.contains(e)).count();
        ends[i].iter().for_each(|&e| if !acc.remove(&e) { acc.insert(e); });
        if shared == 0 { usize::MAX } else { shared }
    }).collect();
    ifs.sort_unstable_by(|a, b| b.cmp(a));
    ifs
}

// Connected-greedy merge order (many pieces): seed thinnest boundary, append fewest-nonzero-shared.
fn greedy_merge_order(ends: &[FxHashSet<Edge>]) -> Vec<usize> {
    let seed = (0..ends.len()).min_by_key(|&i| ends[i].len()).unwrap();
    let mut remaining: Vec<usize> = (0..ends.len()).filter(|&i| i != seed).collect();
    let mut acc = ends[seed].clone();
    let mut order = vec![seed];
    while !remaining.is_empty() {
        let pick = remaining.iter().copied()
            .min_by_key(|&i| match ends[i].iter().filter(|e| acc.contains(e)).count() { 0 => usize::MAX, s => s })
            .unwrap();
        remaining.retain(|&i| i != pick);
        ends[pick].iter().for_each(|&e| if !acc.remove(&e) { acc.insert(e); });
        order.push(pick);
    }
    order
}
