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
    fn divide4(self, i0: usize, j0: usize) -> [CsMat<R>; 4];
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
        let (i0, i1) = (rows.start, rows.end);
        let (j0, j1) = (cols.start, cols.end);

        assert!(i0 <= i1 && i1 <= self.rows());
        assert!(j0 <= j1 && j1 <= self.cols());

        let mut trip = TriMat::new((i1 - i0, j1 - j0));
        
        for (a, (i, j)) in self.into_iter() { 
            let (i, j) = (i.index(), j.index());
            if !a.is_zero() && rows.contains(&i) && cols.contains(&j) {
                trip.add_triplet(i - i0, j - j0, a.clone());
            }
        }

        trip.to_csc()
    }

    fn divide4(self, k: usize, l: usize) -> [CsMat<R>; 4] {
        let (m, n) = self.shape();
        assert!(k <= m);
        assert!(l <= n);

        let mut trip1 = TriMat::new((k, l));
        let mut trip2 = TriMat::new((k, n - l));
        let mut trip3 = TriMat::new((m - k, l));
        let mut trip4 = TriMat::new((m - k, n - l));
        
        for (a, (i, j)) in self.into_iter() { 
            if a.is_zero() { continue }

            let (a, i, j) = (a.clone(), i.index(), j.index());
            match ((0..k).contains(&i), (0..l).contains(&j)) { 
                (true , true ) => trip1.add_triplet(i,     j,     a),
                (true , false) => trip2.add_triplet(i,     j - l, a),
                (false, true ) => trip3.add_triplet(i - k, j,     a),
                (false, false) => trip4.add_triplet(i - k, j - l, a),
            }
        }

        [trip1.to_csc(), trip2.to_csc(), trip3.to_csc(), trip4.to_csc()]
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