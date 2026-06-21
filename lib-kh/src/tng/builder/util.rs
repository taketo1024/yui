use std::ops::RangeInclusive;
use itertools::Itertools;
use rustc_hash::FxHashSet;
use yui_link::{Node, Edge};
use crate::tng::TngComplexKey;

// Pop the least-weight live key from a `(key, weight)` pool. Cached weights go stale as fill grows
// degrees, so re-check the popped key via `live_weight` and re-store it if it got pricier (drop if gone).
pub(crate) fn pop_min_pivot<F>(keys: &mut Vec<(TngComplexKey, usize)>, mut live_weight: F) -> Option<TngComplexKey>
where F: FnMut(&TngComplexKey) -> Option<usize> {
    while !keys.is_empty() {
        let idx = keys.iter().enumerate().min_by_key(|(_, (_, w))| *w).map(|(i, _)| i).unwrap();
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
