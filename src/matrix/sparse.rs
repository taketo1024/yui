use std::ops::{Deref, Add};

use num_traits::Zero;
use sprs::{CsMatBase, SpIndex, TriMat, CsMat};

pub trait CsMatElem: Clone + Default + Send + Sync + sprs::MulAcc {}

impl<T> CsMatElem for T
    where Self: Clone + Default + Send + Sync + sprs::MulAcc 
{}

pub trait CsMatExt { 
    type N;
    fn csc_from_vec(shape: (usize, usize), data: Vec<Self::N>) -> CsMat<Self::N>;
    fn is_zero(&self) -> bool;
}

impl<N, I, IptrStorage, IndStorage, DataStorage, Iptr> CsMatExt for CsMatBase<N, I, IptrStorage, IndStorage, DataStorage, Iptr>
where
    N: Zero + Add<Output = N> + Clone,
    I: SpIndex,
    Iptr: SpIndex,
    IptrStorage: Deref<Target = [Iptr]>,
    IndStorage: Deref<Target = [I]>,
    DataStorage: Deref<Target = [N]>,
{
    type N = N;

    fn csc_from_vec(shape: (usize, usize), data: Vec<Self::N>) -> CsMat<Self::N> {
        let (m, n) = shape;
        assert!(m * n >= data.len());

        let mut trip = TriMat::new(shape);
        for (k, a) in data.into_iter().enumerate() {
            if a.is_zero() { continue }
            let (i, j) = (k / n, k % n);
            trip.add_triplet(i, j, a);
        }
        
        trip.to_csc()
    }

    fn is_zero(&self) -> bool {
        self.data().iter().all(|a| a.is_zero())
    }
}