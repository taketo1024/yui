use std::sync::atomic::{AtomicUsize, Ordering};

/// A thread-safe `usize` counter, backed by [`AtomicUsize`].
pub struct SyncCounter {
    count: AtomicUsize,
}

impl SyncCounter {
    pub fn new() -> Self {
        Self { count: AtomicUsize::new(0) }
    }

    pub fn incr(&self) -> usize {
        self.count.fetch_add(1, Ordering::Relaxed) + 1
    }

    pub fn count(&self) -> usize {
        self.count.load(Ordering::Relaxed)
    }
}
