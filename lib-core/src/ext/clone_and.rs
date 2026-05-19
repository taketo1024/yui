pub trait CloneAnd where Self: Clone { 
    fn clone_and<F>(&self, f: F) -> Self
    where F: FnOnce(&mut Self) { 
        let mut cloned = self.clone();
        f(&mut cloned);
        cloned
    }
}

impl<T> CloneAnd for T where T: Clone {}