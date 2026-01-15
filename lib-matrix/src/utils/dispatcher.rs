use std::{
    any::{Any, TypeId},
    collections::HashMap,
    sync::{Arc, RwLock},
};

pub trait FnSpec: 'static {
    type Input<T: 'static>: Sized + Any;
    type Output<T: 'static>: Any;
}

struct FnObj<S: FnSpec, T: 'static> {
    f: fn(&S::Input<T>) -> S::Output<T>,
}

trait FnCall: Send + Sync + 'static {
    fn call(&self, input: &dyn Any) -> Box<dyn Any>;
}

impl<S: FnSpec, T: 'static> FnCall for FnObj<S, T> {
    fn call(&self, input: &dyn Any) -> Box<dyn Any> {
        let typed = input
            .downcast_ref::<S::Input<T>>()
            .expect("engine called with wrong input type");
        Box::new((self.f)(typed)) as Box<dyn Any>
    }
}

pub struct Dispatcher<S: FnSpec> {
    map: RwLock<HashMap<TypeId, Arc<dyn FnCall>>>,
    _phantom: std::marker::PhantomData<S>,
}

impl<S: FnSpec> Dispatcher<S> {
    pub fn new() -> Self {
        Self {
            map: RwLock::new(HashMap::new()),
            _phantom: std::marker::PhantomData,
        }
    }

    pub fn is_callable<T: 'static>(&self) -> bool { 
        let map = self.map.read().expect("engine registry poisoned");
        map.contains_key(&TypeId::of::<T>())
    }

    pub fn set<T: 'static>(&self, f: fn(&S::Input<T>) -> S::Output<T>) {
        let mut map = self.map.write().expect("engine registry poisoned");
        map.insert(TypeId::of::<T>(), Arc::new(FnObj::<S, T> { f }));
    }

    pub fn remove<T: 'static>(&self) {
        let mut map = self.map.write().expect("engine registry poisoned");
        map.remove(&TypeId::of::<T>());
    }

    pub fn try_call<T: 'static>(&self, input: &S::Input<T>) -> Option<S::Output<T>> {
        // Clone the Arc while holding the read lock; then drop the guard.
        let engine: Arc<dyn FnCall> = {
            let map = self.map.read().expect("engine registry poisoned");
            map.get(&TypeId::of::<T>()).cloned()
        }?;

        let boxed = engine.call(input as &dyn Any);
        boxed.downcast::<S::Output<T>>().ok().map(|b| *b)
    }
}

use std::ptr::NonNull;

#[repr(transparent)]
#[derive(Copy, Clone)]
pub struct RefPtr<T>(NonNull<T>);

impl<T> RefPtr<T> {
    pub fn new(r: &T) -> Self {
        Self(NonNull::from(r))
    }

    /// # Safety
    /// The referenced value must still be alive.
    pub unsafe fn as_ref<'a>(&self) -> &'a T {
        self.0.as_ref()
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Add;
    use super::*;

    struct AddSpec;
    impl FnSpec for AddSpec {
        type Input<T: 'static> = (T, T);
        type Output<T: 'static> = T;
    }

    fn add<R>(input: &(R, R)) -> R
    where for<'a> &'a R: Add<Output = R> { 
        Add::add(&input.0, &input.1)
    }

    #[test]
    fn test_add_dispatch() { 
        let reg = Dispatcher::<AddSpec>::new();

        assert_eq!(reg.is_callable::<i32>(), false);
        assert_eq!(reg.is_callable::<i64>(), false);
        assert_eq!(reg.try_call(&(2, 3)), None);
        assert_eq!(reg.try_call(&(2_i64, 3)), None);

        reg.set(add::<i32>);

        assert_eq!(reg.is_callable::<i32>(), true);
        assert_eq!(reg.is_callable::<i64>(), false);
        assert_eq!(reg.try_call(&(2, 3)), Some(5));
        assert_eq!(reg.try_call(&(2_i64, 3)), None);

        reg.remove::<i32>();

        assert_eq!(reg.is_callable::<i32>(), false);
        assert_eq!(reg.is_callable::<i64>(), false);
        assert_eq!(reg.try_call(&(2, 3)), None);
        assert_eq!(reg.try_call(&(2_i64, 3)), None);
    }

    #[test]
    fn test_add_dispatch_closure() { 
        let reg = Dispatcher::<AddSpec>::new();

        assert_eq!(reg.is_callable::<i32>(), false);
        assert_eq!(reg.is_callable::<i64>(), false);
        assert_eq!(reg.try_call(&(2, 3)), None);
        assert_eq!(reg.try_call(&(2_i64, 3)), None);

        reg.set(|input: &(i32, i32)| input.0 + input.1);

        assert_eq!(reg.is_callable::<i32>(), true);
        assert_eq!(reg.is_callable::<i64>(), false);
        assert_eq!(reg.try_call(&(2, 3)), Some(5));
        assert_eq!(reg.try_call(&(2_i64, 3)), None);
    }

    struct RefAddSpec;
    impl FnSpec for RefAddSpec {
        type Input<T: 'static> = (RefPtr<T>, RefPtr<T>);
        type Output<T: 'static> = T;
    }

    fn add_ref<R>(input: &(RefPtr<R>, RefPtr<R>)) -> R
    where for<'a> &'a R: Add<Output = R> { 
        let (a, b) = unsafe { (input.0.as_ref(), input.1.as_ref()) };
        Add::add(a, b)
    }

    #[test]
    fn test_ref_add_dispatch() { 
        let reg = Dispatcher::<RefAddSpec>::new();
        let a = 2;
        let b = 3;
        let a_ptr = RefPtr::new(&a);
        let b_ptr = RefPtr::new(&b);

        reg.set(add_ref::<i32>);

        assert_eq!(reg.is_callable::<i32>(), true);
        assert_eq!(reg.try_call(&(a_ptr, b_ptr)), Some(5));

        reg.remove::<i32>();

        assert_eq!(reg.is_callable::<i32>(), false);
        assert_eq!(reg.try_call(&(a_ptr, b_ptr)), None);
    }
}