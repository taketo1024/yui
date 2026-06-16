use std::ops::RangeInclusive;
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
