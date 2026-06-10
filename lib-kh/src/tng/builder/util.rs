use std::ops::RangeInclusive;

// Indices in `base` that can still reach `window` with `r` pending crossings: an index
// `i` ends up in `[i, i+r]`, so keep `i ∈ [a-r, b] ∩ base` (all of `base` if no window).
pub(crate) fn reachable_range(base: RangeInclusive<isize>, window: &Option<RangeInclusive<isize>>, r: usize) -> RangeInclusive<isize> {
    match window {
        Some(w) => (*base.start()).max(*w.start() - r as isize) ..= (*base.end()).min(*w.end()),
        None => base,
    }
}
