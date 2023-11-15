use core::sync::atomic::{AtomicBool, Ordering};

static MULTITHREAD_ENABLED: AtomicBool = AtomicBool::new(true);

pub fn is_multithread_enabled() -> bool {
    MULTITHREAD_ENABLED.load(Ordering::Relaxed)
}

pub fn set_multithread_enabled(val: bool) {
    MULTITHREAD_ENABLED.store(val, Ordering::Relaxed)
}