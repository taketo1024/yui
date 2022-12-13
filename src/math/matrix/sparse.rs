use std::ops::{Deref, Add};

use num_traits::{Zero, One};
use rand::Rng;
use sprs::{CsMatBase, SpIndex, TriMat, CsMat, PermView};

pub trait CsMatElem: Clone + Default + Send + Sync + sprs::MulAcc {}

impl<T> CsMatElem for T
    where Self: Clone + Default + Send + Sync + sprs::MulAcc 
{}

pub trait CsMatExt { 
    type N;
    fn csc_from_vec(shape: (usize, usize), data: Vec<Self::N>) -> CsMat<Self::N>;
    fn rand(shape: (usize, usize), density: f64) -> CsMat<Self::N>;
    fn is_zero(&self) -> bool;
    fn permute(self, p: PermView, q: PermView) -> CsMat<Self::N>;
}

impl<N, I, IptrStorage, IndStorage, DataStorage, Iptr> CsMatExt for CsMatBase<N, I, IptrStorage, IndStorage, DataStorage, Iptr>
where
    N: Clone + Zero + One + Add<Output = N>,
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

    fn rand(shape: (usize, usize), density: f64) -> CsMat<Self::N> { 
        let mut rng = rand::thread_rng();
        
        let (m, n) = shape;
        let mut trip = TriMat::new(shape);

        for i in 0..m { 
            for j in 0..n { 
                if rng.gen::<f64>() < density { 
                    trip.add_triplet(i, j, Self::N::one());
                }
            }
        }

        trip.to_csc()
    }

    fn is_zero(&self) -> bool {
        self.data().iter().all(|a| a.is_zero())
    }

    fn permute(self, p: PermView, q: PermView) -> CsMat<Self::N> {
        let mut trip = TriMat::new(self.shape());
        
        for (a, (i, j)) in self.into_iter() { 
            if a.is_zero() { continue }
            let (i, j) = (i.index(), j.index());
            trip.add_triplet(p.at(i), q.at(j), a.clone());
        }

        trip.to_csc()
    }
}