//! Paced progress logging: one line per `step` boundary crossed.

use log::{log, Level};

/// One progress line per `step` boundary that `prev → done` crosses (increments may exceed 1,
/// e.g. equivariant eliminations consume τ-pairs). `depth` indents by two spaces per level.
pub fn log_progress(level: Level, done: usize, prev: usize, total: usize, step: usize, depth: usize) {
    if log_step_crossed(done, prev, total, step) {
        log!(level, "{}...{done}/{total} ({}%)", "  ".repeat(depth), 100 * done / total);
    }
}

/// The pacing test alone, for callers whose line carries more than `done/total`. The final step
/// always reports; a run of at most one step logs nothing.
pub fn log_step_crossed(done: usize, prev: usize, total: usize, step: usize) -> bool {
    total > step && (done / step > prev / step || done == total)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn step_crossing() {
        // a run of at most one step is silent, final step included.
        assert!(!log_step_crossed(10, 9, 10, 10));

        // one line per boundary, not per increment.
        assert!( log_step_crossed(10,  9, 100, 10));
        assert!(!log_step_crossed(11, 10, 100, 10));

        // increments may exceed 1 and must still fire once.
        assert!(log_step_crossed(22, 18, 100, 10));

        // the final step always reports, wherever it falls.
        assert!(log_step_crossed(100, 99, 100, 10));
        assert!(log_step_crossed( 95, 94,  95, 10));
    }
}
