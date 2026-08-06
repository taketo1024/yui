//! [`CloneAnd`]: clone a value, mutate the clone, and return it.

/// Adds `clone_and(f)`: clone `self`, apply the mutation `f`, and return
/// the modified clone — for building variants without mutating the original.
pub trait CloneAnd where Self: Clone {
    fn clone_and<F>(&self, f: F) -> Self
    where F: FnOnce(&mut Self) {
        let mut cloned = self.clone();
        f(&mut cloned);
        cloned
    }
}

impl<T> CloneAnd for T where T: Clone {}