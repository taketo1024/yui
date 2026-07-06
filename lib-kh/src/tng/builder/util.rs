use std::ops::RangeInclusive;
use itertools::Itertools;
use rustc_hash::FxHashSet;
use yui_link::{Node, Edge};
use crate::tng::TngComplexKey;

// Pop the least-weight live key from a `(key, weight)` pool. Cached weights go stale as fill grows
// degrees, so re-check the popped key via `live_weight` and re-store it if it got pricier (drop if gone).
// Ties break by key (not Vec/hash-table order) so the elimination order is identical across platforms.
pub(crate) fn pop_min_pivot<F>(keys: &mut Vec<(TngComplexKey, usize)>, mut live_weight: F) -> Option<TngComplexKey>
where F: FnMut(&TngComplexKey) -> Option<usize> {
    while !keys.is_empty() {
        let idx = keys.iter().enumerate().min_by_key(|(_, (k, w))| (*w, *k)).map(|(i, _)| i).unwrap();
        let (key, w) = keys[idx];
        match live_weight(&key) {
            None => { keys.swap_remove(idx); }
            Some(w_now) if w_now > w => keys[idx].1 = w_now, // pricier than cached → re-store, reselect
            Some(_) => { keys.swap_remove(idx); return Some(key); } // ≤ cached → already the best, take it
        }
    }
    None
}

// Indices in `base` that can still reach `window` with `r` pending crossings: an index
// `i` ends up in `[i, i+r]`, so keep `i ∈ [a-r, b] ∩ base` (all of `base` if no window).
pub(crate) fn reachable_range(base: RangeInclusive<isize>, window: &Option<RangeInclusive<isize>>, r: usize) -> RangeInclusive<isize> {
    match window {
        Some(w) => (*base.start()).max(*w.start() - r as isize) ..= (*base.end()).min(*w.end()),
        None => base,
    }
}

// ---- Cutwidth-profile helpers (for TngComplexBuilder::profile / SymTngBuilder::profile_sym) ----
// A `node_unit` is one crossing (on-axis) or a τ-pair; `open` is the set of boundary (open) arc-ends.

// Unicode bar chart of a width sequence, scaled to `peak`.
pub(crate) fn sparkline(widths: &[usize], peak: usize) -> String {
    const BARS: [char; 8] = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];
    widths.iter().map(|&w| {
        let i = if peak == 0 { 0 } else { w * 7 / peak };
        BARS[i.min(7)]
    }).collect()
}

// Bucket `costs` into 32 log2 bins: bin 0 = cost 0 (free), bin b = cost in [2^(b-1), 2^b);
// costs ≥ 2^31 saturate into the top bin.
pub(crate) fn fill_cost_histogram(costs: impl Iterator<Item = usize>) -> [usize; 32] {
    let mut hist = [0usize; 32];
    for c in costs {
        let b = if c == 0 { 0 } else { (usize::BITS - c.leading_zeros()) as usize };
        hist[b.min(31)] += 1;
    }
    hist
}

// Sparkline of a target set's fill-cost distribution (log2 buckets). Kept out of the `debug!` call
// site so the histogram is built only when debug logging is on. A `|` marks the cap — but only when
// it actually splits the distribution (some edges above it defer to the matrix); if the cap sits at
// or above the top non-empty bucket, everything is kept and no marker is drawn.
pub(crate) fn fill_cost_sparkline(keys: &[(TngComplexKey, usize)], max: Option<usize>) -> String {
    let hist = fill_cost_histogram(keys.iter().map(|(_, c)| *c));
    let hi = hist.iter().rposition(|&c| c > 0).unwrap_or(0);
    let peak = hist.iter().copied().max().unwrap_or(0);
    let bars = sparkline(&hist[..=hi], peak);
    match max {
        // last fully-kept bucket = ⌊log2(cap+1)⌋ (bucket b covers [2^(b-1), 2^b), kept iff 2^b-1 ≤ cap).
        Some(m) if (m + 1).ilog2() < hi as u32 => {
            let cut = (m + 1).ilog2() as usize;
            bars.chars().enumerate()
                .flat_map(|(b, ch)| if b == cut { vec![ch, '|'] } else { vec![ch] })
                .collect()
        }
        _ => bars,
    }
}

// The node-unit's boundary arc-ends: its edges with odd incidence (one endpoint inside the unit).
pub(crate) fn boundary_edges(node_unit: &[&Node]) -> Vec<Edge> {
    node_unit.iter().flat_map(|x| x.edges().iter().copied()).counts().into_iter()
        .filter_map(|(e, c)| (c % 2 == 1).then_some(e))
        .collect()
}

// Boundary cutwidth that appending `node_unit` would yield, without mutating `open`.
pub(crate) fn cutwidth_after(open: &FxHashSet<Edge>, node_unit: &[&Node]) -> usize {
    let delta: isize = boundary_edges(node_unit).iter()
        .map(|e| if open.contains(e) { -1 } else { 1 })
        .sum();
    (open.len() as isize + delta) as usize
}

// Toggle `node_unit`'s flip-edges into/out of the open-edge set (open ↦ open △ boundary_edges).
pub(crate) fn toggle_boundary(open: &mut FxHashSet<Edge>, node_unit: &[&Node]) {
    boundary_edges(node_unit).into_iter().for_each(|e| {
        if !open.remove(&e) {
            open.insert(e);
        }
    });
}

// ---- Chunk planning: pick the cut positions in a cutwidth profile (shared by both builders) ----

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
pub(crate) fn select_cuts(widths: &[usize], n_cuts: usize) -> Vec<usize> {
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

// ---- Manual edge-cut chunking (shared by both builders; the sym builder adds its own τ-checks) ----

// Crossing components after severing the `cut` edges (BFS over non-cut shared edges).
pub(crate) fn cut_components(nodes: &[Node], cut: &FxHashSet<Edge>) -> Vec<Vec<usize>> {
    let mut by_edge: rustc_hash::FxHashMap<Edge, Vec<usize>> = Default::default();
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
pub(crate) fn merge_order(nodes: &[Node], pieces: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
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
