use std::ops::{Deref, Add, Range};

use num_traits::{Zero, One};
use rand::Rng;
use sprs::{CsMatBase, SpIndex, TriMat, CsMat, PermView, CsVec};

pub trait CsMatElem: Clone + Default + Send + Sync + sprs::MulAcc {}

impl<T> CsMatElem for T
    where Self: Clone + Default + Send + Sync + sprs::MulAcc 
{}

pub trait CsMatExt<R> { 
    fn csc_from_vec(shape: (usize, usize), data: Vec<R>) -> CsMat<R>;
    fn rand(shape: (usize, usize), density: f64) -> CsMat<R>;
    fn is_zero(&self) -> bool;
    fn permute(self, p: PermView, q: PermView) -> CsMat<R>;
    fn submatrix(&self, rows: Range<usize>, cols: Range<usize>) -> CsMat<R>;
}

impl<R, I, IptrStorage, IndStorage, DataStorage, Iptr> CsMatExt<R> for CsMatBase<R, I, IptrStorage, IndStorage, DataStorage, Iptr>
where
    R: Clone + Zero + One + Add<Output = R>,
    I: SpIndex,
    Iptr: SpIndex,
    IptrStorage: Deref<Target = [Iptr]>,
    IndStorage: Deref<Target = [I]>,
    DataStorage: Deref<Target = [R]>,
{
    fn csc_from_vec(shape: (usize, usize), data: Vec<R>) -> CsMat<R> {
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

    fn rand(shape: (usize, usize), density: f64) -> CsMat<R> { 
        let mut rng = rand::thread_rng();
        
        let (m, n) = shape;
        let mut trip = TriMat::new(shape);

        for i in 0..m { 
            for j in 0..n { 
                if rng.gen::<f64>() < density { 
                    trip.add_triplet(i, j, R::one());
                }
            }
        }

        trip.to_csc()
    }

    fn is_zero(&self) -> bool {
        self.data().iter().all(|a| a.is_zero())
    }

    fn permute(self, p: PermView, q: PermView) -> CsMat<R> {
        let mut trip = TriMat::new(self.shape());
        
        for (a, (i, j)) in self.into_iter() { 
            if a.is_zero() { continue }
            let (i, j) = (i.index(), j.index());
            trip.add_triplet(p.at(i), q.at(j), a.clone());
        }

        trip.to_csc()
    }

    fn submatrix(&self, rows: Range<usize>, cols: Range<usize>) -> CsMat<R> {
        todo!()
    }
}

pub trait CsVecExt<R> { 
    fn from_vec(data: Vec<R>) -> CsVec<R>;
}

impl<R> CsVecExt<R> for CsVec<R> 
where R: Zero {
    fn from_vec(data: Vec<R>) -> CsVec<R> {
        let n = data.len();
        let (indices, data): (Vec<usize>, Vec<R>) = data
            .into_iter()
            .enumerate()
            .filter_map(|(i, a)| 
                if a.is_zero() { None } else { Some((i, a)) }
            ).unzip();
        
        Self::new(n, indices, data)
    }
}