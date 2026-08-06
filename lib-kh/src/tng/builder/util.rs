use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::ops::RangeInclusive;
use itertools::Itertools;
use rustc_hash::FxHashSet;
use yui_link::{Node, Edge, InvLink};
use crate::tng::TngComplexKey;

// A lazy min-priority pool of pivot candidates, ordered by (cached weight, key). The key joins the
// ordering so the pop sequence is deterministic across platforms (no true ties — keys are unique).
pub(crate) type PivotPool = BinaryHeap<Reverse<(usize, TngComplexKey)>>;

pub(crate) fn pivot_pool(keys: impl IntoIterator<Item = (TngComplexKey, usize)>) -> PivotPool {
    keys.into_iter().map(|(k, w)| Reverse((w, k))).collect()
}

pub(crate) fn push_pivot(pool: &mut PivotPool, key: TngComplexKey, weight: usize) {
    pool.push(Reverse((weight, key)));
}

// Pop the least-weight live key. Cached weights go stale as fill grows degrees, so re-check the
// popped key via `live_weight`: if it got pricier, re-push at the new weight and reselect; drop it
// if gone. O(log n) per pop, vs a linear scan of the whole frontier.
pub(crate) fn pop_min_pivot<F>(pool: &mut PivotPool, mut live_weight: F) -> Option<TngComplexKey>
where F: FnMut(&TngComplexKey) -> Option<usize> {
    while let Some(Reverse((w, key))) = pool.pop() {
        match live_weight(&key) {
            None => {}                                                    // gone → drop
            Some(w_now) if w_now > w => pool.push(Reverse((w_now, key))), // pricier than cached → re-push, reselect
            Some(_) => return Some(key),                                  // ≤ cached → already the best, take it
        }
    }
    None
}

// v1.0 supports only strongly invertible links in a transvergent diagram: canon cycles need not
// be 1 + τ closed otherwise, and `partition_off_axis` needs an axis that separates the plane.
pub(crate) fn assert_supported_symmetry(l: &InvLink) {
    assert!(
        l.is_strongly_invertible(),
        "currently, only strongly invertible knots / links are supported"
    );
    assert!(
        l.is_transvergent(),
        "currently, only transvergent diagrams are supported (the axis must lie in the plane)"
    );
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
    let peak_h = hist.iter().copied().max().unwrap_or(0);
    let bars = sparkline(&hist[..=hi], peak_h);
    let bars: String = match max {
        // last fully-kept bucket = ⌊log2(cap+1)⌋ (bucket b covers [2^(b-1), 2^b), kept iff 2^b-1 ≤ cap).
        Some(m) if (m + 1).ilog2() < hi as u32 => {
            let cut = (m + 1).ilog2() as usize;
            bars.chars().enumerate()
                .flat_map(|(b, ch)| if b == cut { vec![ch, '|'] } else { vec![ch] })
                .collect()
        }
        _ => bars,
    };

    // total eliminatable edges, modal cost (tallest bucket b holds costs 2^(b-1)..2^b), and max cost.
    let total = keys.len();
    let max_cost = keys.iter().map(|(_, c)| *c).max().unwrap_or(0);
    let peak_bucket = hist.iter().enumerate().max_by_key(|&(_, &c)| c).map_or(0, |(b, _)| b);
    let peak = if peak_bucket == 0 { "0".to_string() } else { format!("2^{}", peak_bucket - 1) };
    let max_s = if max_cost == 0 { "0".to_string() } else { format!("2^{}", max_cost.ilog2()) };
    format!("{bars} (total: {total}, peak: {peak}, max: {max_s})")
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

