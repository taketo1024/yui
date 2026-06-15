use std::ops::RangeInclusive;
use yui_core::bitseq::Bit;
use yui_link::{Node, Edge};
use crate::tng::TngComp;

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
