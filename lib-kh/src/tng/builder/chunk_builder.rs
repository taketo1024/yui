//! Divide-and-conquer chunked build, shared by [`TngComplexBuilder`](super::TngComplexBuilder)
//! and [`SymTngBuilder`](super::SymTngBuilder): plan k chunks up front at thin cutwidth
//! interfaces, build each via a child builder, and merge the reduced chunks into the parent —
//! so the parent never materializes the full dense slice.

use itertools::Itertools;
use log::info;
use rustc_hash::{FxHashMap, FxHashSet};
use yui_link::{Edge, Node};

use super::{boundary_edges, CutOption};

/// The planning view of a boundary-cutwidth profile: node-index units in build order (1 node
/// plain / on-axis, 2 for a τ-pair) and the boundary cutwidth after each unit. Cuts land
/// between units, so units stay whole.
pub(crate) struct PlanProfile {
    pub units: Vec<Vec<usize>>,
    pub widths: Vec<usize>,
}

/// A builder that can be chunk-built: planning inputs plus the child build.
pub(crate) trait Chunkable: Sized {
    type Elem;

    fn nodes(&self) -> &[Node];
    fn cut_option(&self) -> &CutOption;
    fn profile(&self) -> PlanProfile;
    fn stat(&self) -> String;

    /// A child builder over `chunk` (a sub-tangle).
    fn init_child(&self, chunk: &[Node]) -> Self;
    fn run(self) -> Self;
    fn take_elements(&mut self) -> Vec<Self::Elem>;

    /// Manual-plan hook: fuse the cut pieces (τ-orbit fusion for the sym builder).
    fn process_pieces(&self, pieces: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
        pieces
    }

    /// Manual-plan hook: validate the cut and its pieces (τ-symmetry for the sym builder).
    fn validate_cut(&self, _cut: &FxHashSet<Edge>, _pieces: &[Vec<usize>]) {}
}

pub(crate) struct ChunkBuilder<'a, B> {
    builder: &'a B,
}

impl<'a, B: Chunkable> ChunkBuilder<'a, B> {
    pub(crate) fn new(builder: &'a B) -> Self {
        Self { builder }
    }

    // Plan the chunks and build each into a reduced sub-complex, paired with the crossings it
    // covers. Building is independent of the parent, so the caller merges them afterwards.
    pub(crate) fn build_chunks(&self) -> Vec<(Vec<Node>, B, Vec<B::Elem>)> {
        let plan = self.plan();
        info!("chunk plan: {} pieces {:?}", plan.len(),
            plan.iter().map(|c| c.len()).collect_vec());

        plan.into_iter().map(|chunk| {
            let (child, elems) = self.build_chunk(&chunk);
            (chunk, child, elems)
        }).collect()
    }

    // Partition the crossings into pieces: `Auto(k)` at the deepest cutwidth valleys of the
    // MinCut order (contiguous, thin interface), `AtCrossings` at the requested crossing counts,
    // or `Manual` by severing the given edge-cut(s).
    fn plan(&self) -> Vec<Vec<Node>> {
        let k = match self.builder.cut_option() {
            CutOption::Manual(cuts) => return self.manual_plan(cuts),
            CutOption::AtCrossings(counts) => return self.at_crossings_plan(counts),
            CutOption::Auto(k) => (*k).max(1),
            CutOption::None => 1,
        };
        let prof = self.builder.profile();
        let cuts = select_cuts(&prof.widths, k - 1);
        self.segment_plan(&prof, &cuts)
    }

    // Cut after the unit positions whose cumulative crossing count is closest to each requested
    // count — direct control over chunk balance (units stay whole, so counts land within ±1).
    fn at_crossings_plan(&self, counts: &[usize]) -> Vec<Vec<Node>> {
        let prof = self.builder.profile();
        let cum: Vec<usize> = prof.units.iter()
            .scan(0, |acc, unit| {
                *acc += unit.len();
                Some(*acc)
            })
            .collect();

        // the final position is excluded — cutting after everything leaves an empty tail chunk.
        let cuts = counts.iter()
            .filter_map(|&c| (0..cum.len().saturating_sub(1)).min_by_key(|&p| cum[p].abs_diff(c)))
            .sorted().dedup().collect_vec();
        for (&c, &p) in counts.iter().sorted().zip(cuts.iter()) {
            info!("cut requested at {c} crossings -> position {p} ({} crossings, width {})", cum[p], prof.widths[p]);
        }
        self.segment_plan(&prof, &cuts)
    }

    // Segment the profile order after each cut position; expand each unit to its nodes.
    fn segment_plan(&self, prof: &PlanProfile, cuts: &[usize]) -> Vec<Vec<Node>> {
        let nodes = self.builder.nodes();
        let starts = std::iter::once(0).chain(cuts.iter().map(|&v| v + 1));
        let ends = cuts.iter().map(|&v| v + 1).chain(std::iter::once(prof.units.len()));
        starts.zip(ends)
            .map(|(s, e)| prof.units[s..e].iter().flatten().map(|&i| nodes[i].clone()).collect())
            .filter(|c: &Vec<Node>| !c.is_empty())
            .collect()
    }

    // Manual cut: sever the cut edges (union of all cut-lines), let the builder fuse and
    // validate the pieces, then order them for merging.
    fn manual_plan(&self, cuts: &[Vec<Edge>]) -> Vec<Vec<Node>> {
        let cut: FxHashSet<Edge> = cuts.iter().flatten().copied().collect();
        let nodes = self.builder.nodes();
        let comps = cut_components(nodes, &cut);
        let pieces = self.builder.process_pieces(comps);
        self.builder.validate_cut(&cut, &pieces);
        assert!(pieces.len() >= 2, "cut does not separate the link into ≥2 pieces");
        merge_order(nodes, pieces).into_iter()
            .map(|piece| piece.into_iter().map(|i| nodes[i].clone()).collect())
            .collect()
    }

    // Build `chunk` into a reduced sub-complex via a child builder.
    fn build_chunk(&self, chunk: &[Node]) -> (B, Vec<B::Elem>) {
        let ends = boundary_edges(&chunk.iter().collect::<Vec<_>>()).into_iter().sorted().collect_vec();
        info!("build chunk (n: {}, nb: {} {:?}): {}", chunk.len(), ends.len(), ends, chunk.iter().join(", "));

        let mut child = self.builder.init_child(chunk).run();
        info!("chunk built: {}", child.stat());
        let elems = child.take_elements();
        (child, elems)
    }
}

// ---- Auto chunk planning: pick the cut positions in a cutwidth profile ----

// Valley positions (a descent that turns back up) in `widths`, ascending. The `scan` carries a
// `descended` flag so flats are ignored — a mid-descent plateau isn't mistaken for the bottom.
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

// `n_cuts` cut positions in `widths`: take the deepest valleys first; if more cuts than valleys
// are needed, keep every valley and split the widest pieces evenly (greedy, which minimizes the
// largest piece).
fn select_cuts(widths: &[usize], n_cuts: usize) -> Vec<usize> {
    let n = widths.len();
    let n_cuts = n_cuts.min(n.saturating_sub(1));
    if n_cuts == 0 {
        return vec![];
    }

    let vs = valleys(widths); // ascending positions
    if vs.len() >= n_cuts {
        // enough valleys: keep the `n_cuts` deepest, back in position order.
        return vs.into_iter()
            .sorted_by_key(|&p| widths[p])
            .take(n_cuts)
            .sorted()
            .collect();
    }

    // too few valleys: every valley is a cut; spend the rest splitting the widest pieces evenly.
    // `alloc[i]` is the number of pieces segment `i` is divided into.
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

// ---- Manual edge-cut chunking ----

// Crossing components after severing the `cut` edges (BFS over non-cut shared edges).
fn cut_components(nodes: &[Node], cut: &FxHashSet<Edge>) -> Vec<Vec<usize>> {
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

// Boundary edges of a node-index piece — its merge interface (edges with one endpoint inside).
fn chunk_ends(nodes: &[Node], piece: &[usize]) -> FxHashSet<Edge> {
    let subset: Vec<&Node> = piece.iter().map(|&i| &nodes[i]).collect();
    boundary_edges(&subset).into_iter().collect()
}

// Order pieces so each merges behind a thin interface. Exhaustive over all `k!` orders while that's
// cheap (`8! ≈ 40k`); past 8 pieces `k!` explodes, so fall back to a greedy order.
fn merge_order(nodes: &[Node], pieces: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    let k = pieces.len();
    if k <= 2 {
        return pieces;
    }
    let ends: Vec<FxHashSet<Edge>> = pieces.iter().map(|c| chunk_ends(nodes, c)).collect();
    let order = if k <= 8 {
        (0..k).permutations(k).min_by_key(|o| interface_profile(o, &ends)).unwrap()
    } else {
        greedy_merge_order(&ends)
    };
    order.into_iter().map(|i| pieces[i].clone()).collect()
}

// Sorted-descending shared-edge counts as each piece merges in `order` (0 = disjoint → worst).
// Compared lexicographically to pick the order whose merges share the most boundary throughout.
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

// Greedy merge order for >8 pieces: seed at the thinnest boundary, then always append the piece
// sharing the most edges with the accumulated frontier.
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
