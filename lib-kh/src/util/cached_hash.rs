use std::cmp::Ordering as CmpOrdering;
use std::fmt::{self, Debug, Display};
use std::hash::{Hash, Hasher};
use std::ops::Deref;
use std::sync::atomic::{AtomicU64, Ordering};

use rustc_hash::FxHasher;

// A value paired with a lazily-cached structural hash. Mutate the payload only via `inner_mut`,
// which invalidates the cache; `Eq`/`Hash` then short-circuit on the cached `u64`. Used as the
// hashed payload of `Cob`/`Tng`, whose values are hashed far more often than they're mutated.
pub(crate) struct CachedHash<T> {
    value: T,
    hash: AtomicU64, // 0 = uncomputed (a real hash of 0 is bumped to 1)
}

impl<T: Hash> CachedHash<T> {
    pub(crate) fn new(value: T) -> Self {
        Self { value, hash: AtomicU64::new(0) }
    }

    pub(crate) fn inner_mut(&mut self) -> &mut T {
        *self.hash.get_mut() = 0; // any mutation may change the hash → invalidate
        &mut self.value
    }

    pub(crate) fn into_inner(self) -> T {
        self.value
    }

    pub(crate) fn cached_hash(&self) -> u64 {
        let h = self.hash.load(Ordering::Relaxed);
        if h != 0 {
            return h;
        }
        let mut hasher = FxHasher::default();
        self.value.hash(&mut hasher);
        let h = hasher.finish().max(1);
        self.hash.store(h, Ordering::Relaxed);
        h
    }
}

impl<T: Clone> Clone for CachedHash<T> {
    fn clone(&self) -> Self {
        Self { value: self.value.clone(), hash: AtomicU64::new(self.hash.load(Ordering::Relaxed)) }
    }
}

impl<T: Hash + PartialEq> PartialEq for CachedHash<T> {
    fn eq(&self, other: &Self) -> bool {
        self.cached_hash() == other.cached_hash() && self.value == other.value
    }
}

impl<T: Hash + Eq> Eq for CachedHash<T> {}

// Order by the payload only — the cached `u64` is incidental. (Can't `#[derive]`: `AtomicU64`
// isn't `Ord`.) Consistent with `Eq`: equal payloads ⇒ `Equal`, and equal payloads hash equally.
impl<T: Hash + PartialOrd> PartialOrd for CachedHash<T> {
    fn partial_cmp(&self, other: &Self) -> Option<CmpOrdering> {
        self.value.partial_cmp(&other.value)
    }
}

impl<T: Hash + Ord> Ord for CachedHash<T> {
    fn cmp(&self, other: &Self) -> CmpOrdering {
        self.value.cmp(&other.value)
    }
}

impl<T: Hash> Hash for CachedHash<T> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(self.cached_hash());
    }
}

impl<T: Hash + Default> Default for CachedHash<T> {
    fn default() -> Self {
        Self::new(T::default())
    }
}

// Delegate to the payload — don't surface the cached `u64`.
impl<T: Debug> Debug for CachedHash<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.value.fmt(f)
    }
}

impl<T: Display> Display for CachedHash<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.value.fmt(f)
    }
}

// Read-only deref for the payload's read API. No `DerefMut` on purpose: writes go through
// `inner_mut` so the (rare) cache-resetting sites stay explicit.
impl<T> Deref for CachedHash<T> {
    type Target = T;
    fn deref(&self) -> &T {
        &self.value
    }
}
