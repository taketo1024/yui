use std::ops::RangeInclusive;
use itertools::Itertools;
use rustc_hash::FxHashSet;
use yui_core::bitseq::Bit;
use yui_link::{Node, Edge};
use crate::tng::{TngComp, TngComplexKey};

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

// Arcs added when appending `x` to the partial diagram, with the base-point edge filtered out.
pub(crate) fn node_arcs(x: &Node, base_pt: Option<Edge>) -> Vec<TngComp> {
    let arcs = if x.is_resolved() {
        let (a0, a1) = x.arcs();
        vec![a0, a1]
    } else {
        let (a00, a01) = x.resolve(Bit::Bit0).arcs();
        let (a10, a11) = x.resolve(Bit::Bit1).arcs();
        vec![a00, a01, a10, a11]
    };
    arcs.into_iter()
        .filter(|a| base_pt.map(|e| !a.contains(e)).unwrap_or(true))
        .map(TngComp::from)
        .collect()
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

// Distinct edges of a node-unit (one or two crossings) with odd incidence — those whose open/closed flips.
fn flip_edges(node_unit: &[&Node]) -> Vec<Edge> {
    node_unit.iter().flat_map(|x| x.edges().iter().copied()).counts().into_iter()
        .filter_map(|(e, c)| (c % 2 == 1).then_some(e))
        .collect()
}

// Boundary cutwidth that appending `node_unit` would yield, without mutating `open`.
pub(crate) fn cutwidth_after(open: &FxHashSet<Edge>, node_unit: &[&Node]) -> usize {
    let delta: isize = flip_edges(node_unit).iter()
        .map(|e| if open.contains(e) { -1 } else { 1 })
        .sum();
    (open.len() as isize + delta) as usize
}

// Toggle `node_unit`'s flip-edges into/out of the open-edge set (open ↦ open △ flip_edges).
pub(crate) fn toggle_boundary(open: &mut FxHashSet<Edge>, node_unit: &[&Node]) {
    flip_edges(node_unit).into_iter().for_each(|e| {
        if !open.remove(&e) {
            open.insert(e);
        }
    });
}
