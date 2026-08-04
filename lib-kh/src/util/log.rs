use log::debug;

// Threshold below which progress lines are suppressed — short rounds don't need pacing.
const PROGRESS_LOG_MIN: usize = 50_000;

// One progress line per `step` boundary that `prev → done` crosses (increments may exceed 1,
// e.g. equivariant eliminations consume τ-pairs), for rounds larger than `PROGRESS_LOG_MIN`.
pub(crate) fn log_progress(done: usize, prev: usize, total: usize, step: usize) {
    if total > PROGRESS_LOG_MIN && done / step > prev / step {
        debug!("    ...{done}/{total} ({}%)", 100 * done / total);
    }
}
