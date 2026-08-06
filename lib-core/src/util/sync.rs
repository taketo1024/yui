//! [`SyncCounter`]: a thread-safe counter used to label objects during a build.

use std::sync::atomic::{AtomicUsize, Ordering};

/// A thread-safe `usize` counter, backed by [`AtomicUsize`].
pub struct SyncCounter {
    count: AtomicUsize,
}

impl SyncCounter {
    pub fn new(n: usize) -> Self {
        Self { count: AtomicUsize::new(n) }
    }

    pub fn count(&self) -> usize {
        self.count.load(Ordering::Relaxed)
    }

    pub fn incr(&self) -> usize {
        self.count.fetch_add(1, Ordering::Relaxed) + 1
    }

    pub fn set(&self, n: usize) {
        self.count.store(n, Ordering::Relaxed)
    }
}
