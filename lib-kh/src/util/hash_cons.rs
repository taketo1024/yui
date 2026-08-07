//! [`HashCons`]: an interning table handing out one shared `Arc` per distinct
//! value.

use std::hash::Hash;
use std::sync::Arc;

use dashmap::DashSet;
use rustc_hash::{FxBuildHasher, FxHashSet};

// Concurrent hash-cons (interning): structurally-equal values collapse to a single shared `Arc<T>`,
// so duplicates share one allocation and compare cheaply. `intern` hits the shared table; on
// read-heavy hot paths use `intern_cached` with a per-thread `Cache` to skip the shard lock+probe
// (the storage must be caller-side — `thread_local!` is a static — but the logic lives here).
pub(crate) struct HashCons<T> {
    table: DashSet<Arc<T>, FxBuildHasher>,
}

// Per-thread read cache for `intern_cached`. Holds the canonical `Arc`s a thread has seen.
pub(crate) type Cache<T> = FxHashSet<Arc<T>>;

impl<T: Eq + Hash> HashCons<T> {
    pub(crate) fn new() -> Self {
        Self { table: DashSet::default() }
    }

    // Canonical `Arc<T>` for `value`, inserting it if unseen.
    pub(crate) fn intern(&self, value: T) -> Arc<T> {
        if let Some(a) = self.table.get(&value) {
            return Arc::clone(&a);
        }
        let arc = Arc::new(value);
        // `insert` returns false if an equal value was interned concurrently; fetch that canonical
        // `Arc` so dedup stays exact under races.
        if !self.table.insert(Arc::clone(&arc)) {
            return Arc::clone(&self.table.get(&*arc).unwrap());
        }
        arc
    }

    // `intern`, served from a per-thread `cache` first. Only canonical `Arc`s enter `cache`, so
    // dedup stays exact while hot repeats avoid the table's lock entirely.
    pub(crate) fn intern_cached(&self, value: T, cache: &mut Cache<T>) -> Arc<T> {
        if let Some(arc) = cache.get(&value) {
            return Arc::clone(arc);
        }
        let arc = self.intern(value);
        cache.insert(Arc::clone(&arc));
        arc
    }
}

impl<T: Eq + Hash> Default for HashCons<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dedups_to_one_arc() {
        let cons = HashCons::<String>::new();
        let a = cons.intern("x".to_string());
        let b = cons.intern("x".to_string());
        let c = cons.intern("y".to_string());
        // equal values share the same allocation; distinct ones don't.
        assert!(Arc::ptr_eq(&a, &b));
        assert!(!Arc::ptr_eq(&a, &c));
        assert_eq!(*a, "x");
    }

    #[test]
    fn cached_returns_canonical_arc() {
        let cons = HashCons::<String>::new();
        let mut cache = Cache::default();
        let a = cons.intern("x".to_string());                       // populate the table
        let b = cons.intern_cached("x".to_string(), &mut cache);    // miss → table, then cached
        let c = cons.intern_cached("x".to_string(), &mut cache);    // hit → cache
        // cached path returns the same canonical `Arc` as the table.
        assert!(Arc::ptr_eq(&a, &b));
        assert!(Arc::ptr_eq(&b, &c));
    }
}
