use std::ops::{Add, AddAssign, Neg, Deref, Range, Sub, SubAssign, Mul, MulAssign};
use std::iter::zip;
use num_traits::{Zero, One};
use rand::Rng;
use sprs::{CsMatBase, SpIndex, TriMat, CsMat, PermView, CsVec, CsVecBase};
use auto_impl_ops::auto_ops;
use crate::math::traits::{Ring, RingOps, AddMonOps, AddGrpOps, MonOps};

use super::dense::MatType;

pub struct SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    cs_mat: CsMat<R>
}

impl<R> SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn shape(&self) -> (usize, usize) {
        self.cs_mat.shape()
    }
}

impl<R> From<CsMat<R>> for SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn from(cs_mat: CsMat<R>) -> Self {
        assert!(cs_mat.is_csc());
        Self { cs_mat }
    }
}

impl<R> Neg for SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<R> Neg for &SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    type Output = SpMat<R>;
    fn neg(self) -> Self::Output {
        let neg = self.cs_mat.map(|a| -a);
        SpMat::from(neg)
    }
}

macro_rules! impl_op {
    ($trait:ident, $method:ident) => {
        #[auto_ops]
        impl<'a, 'b, R> $trait<&'b SpMat<R>> for &'a SpMat<R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            type Output = SpMat<R>;
            fn $method(self, rhs: &'b SpMat<R>) -> Self::Output {
                let res = self.cs_mat.$method(&rhs.cs_mat);
                SpMat::from(res)
            }
        }
    };
}

impl_op!(Add, add);
impl_op!(Sub, sub);
impl_op!(Mul, mul);

macro_rules! impl_ops {
    ($trait:ident) => {
        impl<R> $trait<SpMat<R>> for SpMat<R>
        where R: Ring, for<'x> &'x R: RingOps<R> {}

        impl<R> $trait<SpMat<R>> for &SpMat<R>
        where R: Ring, for<'x> &'x R: RingOps<R> {}
    };
}

impl_ops!(AddMonOps);
impl_ops!(AddGrpOps);
impl_ops!(MonOps);
impl_ops!(RingOps);

impl<R> MatType for SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;
}

pub trait CsMatExt<R> { 
    fn id(n: usize) -> CsMat<R>;
    fn csc_from_vec(shape: (usize, usize), data: Vec<R>) -> CsMat<R>;
    fn rand(shape: (usize, usize), density: f64) -> CsMat<R>;
    fn is_zero_shape(&self) -> bool;
    fn is_zero(&self) -> bool;
    fn is_square(&self) -> bool;
    fn is_id(&self) -> bool;
    fn permute(&self, p: PermView, q: PermView) -> CsMat<R>;
    fn permute_rows(&self, p: PermView) -> CsMat<R>;
    fn permute_cols(&self, q: PermView) -> CsMat<R>;
    fn submatrix(&self, rows: Range<usize>, cols: Range<usize>) -> CsMat<R>;
    fn split_hor(&self, k: usize) -> [CsMat<R>; 2];
    fn split_ver(&self, k: usize) -> [CsMat<R>; 2];
    fn divide4(&self, i0: usize, j0: usize) -> [CsMat<R>; 4];
    fn combine4(blocks: [&CsMat<R>; 4]) -> CsMat<R>;
    fn concat(a: &CsMat<R>, b: &CsMat<R>) -> CsMat<R>;
    fn stack(a: &CsMat<R>, b: &CsMat<R>) -> CsMat<R>;
}

impl<R, I, IptrStorage, IndStorage, DataStorage, Iptr> CsMatExt<R> for CsMatBase<R, I, IptrStorage, IndStorage, DataStorage, Iptr>
where
    R: Clone + Zero + One + Add<Output = R> + PartialEq,
    I: SpIndex,
    Iptr: SpIndex,
    IptrStorage: Deref<Target = [Iptr]>,
    IndStorage: Deref<Target = [I]>,
    DataStorage: Deref<Target = [R]>,
{
    fn id(n: usize) -> CsMat<R> { 
        let indptr = (0..=n).collect();
        let indices = (0..n).collect();
        let data = vec![R::one(); n];
        CsMat::new_csc((n, n), indptr, indices, data)
    }

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

    fn is_zero_shape(&self) -> bool {
        self.rows() == 0 || self.cols() == 0
    }

    fn is_zero(&self) -> bool {
        self.is_zero_shape() || 
        self.data().iter().all(|a| a.is_zero())
    }

    fn is_square(&self) -> bool {
        self.rows() == self.cols()
    }

    fn is_id(&self) -> bool {
        self.is_square() && self.into_iter().all(|(a, (i, j))| 
            (i == j && a.is_one()) || (i != j && a.is_zero())
        )
    }

    fn permute(&self, p: PermView, q: PermView) -> CsMat<R> {
        let mut trip = TriMat::new(self.shape());
        
        for (a, (i, j)) in self.into_iter() { 
            if a.is_zero() { continue }
            let (i, j) = (i.index(), j.index());
            trip.add_triplet(p.at(i), q.at(j), a.clone());
        }

        trip.to_csc()
    }

    fn permute_rows(&self, p: PermView) -> CsMat<R> {
        let id = PermView::identity(self.cols());
        self.permute(p, id)
    }

    fn permute_cols(&self, q: PermView) -> CsMat<R> {
        let id = PermView::identity(self.rows());
        self.permute(id, q)
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

    fn split_hor(&self, k: usize) -> [CsMat<R>; 2] {
        let [a, b, _, _] = self.divide4(self.rows(), k);
        [a, b]
    }

    fn split_ver(&self, k: usize) -> [CsMat<R>; 2] {
        let [a, _, c, _] = self.divide4(k, self.cols());
        [a, c]
    }

    fn divide4(&self, k: usize, l: usize) -> [CsMat<R>; 4] {
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

    fn combine4(blocks: [&CsMat<R>; 4]) -> CsMat<R> {
        let [a, b, c, d] = blocks;

        assert_eq!(a.rows(), b.rows());
        assert_eq!(c.rows(), d.rows());
        assert_eq!(a.cols(), c.cols());
        assert_eq!(b.cols(), d.cols());

        let (m, n) = (a.rows() + c.rows(), a.cols() + b.cols());
        let (k, l) = a.shape();

        let mut trip = TriMat::<R>::new((m, n));

        for (x, (di, dj)) in zip([a,b,c,d], [(0,0),(0,l),(k,0),(k,l)]) { 
            for (r, (i, j)) in x.into_iter() {
                trip.add_triplet(i + di, j + dj, r.clone());
            }
        }

        trip.to_csc()
    }

    fn concat(a: &CsMat<R>, b: &CsMat<R>) -> CsMat<R> {
        let zero = |m, n| CsMat::<R>::zero((m, n));
        Self::combine4([
            a, b, 
            &zero(0, a.cols()), &zero(0, b.cols())
        ])
    }

    fn stack(a: &CsMat<R>, b: &CsMat<R>) -> CsMat<R> {
        let zero = |m, n| CsMat::<R>::zero((m, n));
        Self::combine4([
            a, &zero(a.rows(), 0), 
            b, &zero(b.rows(), 0)
        ])
    }

}

pub trait CsVecExt<R> { 
    fn from_vec(data: Vec<R>) -> CsVec<R>;
    fn is_zero(&self) -> bool;
    fn permute(&self, p: PermView) -> CsVec<R>;
    fn subvec(&self, range: Range<usize>) -> CsVec<R>;
    fn divide2(self, r: usize) -> (CsVec<R>, CsVec<R>);
}

impl<IStorage, DStorage, R, I: SpIndex> CsVecExt<R> for CsVecBase<IStorage, DStorage, R, I>
where
    R: Clone + Zero,
    IStorage: Deref<Target = [I]>,
    DStorage: Deref<Target = [R]>
{
   fn from_vec(data: Vec<R>) -> CsVec<R> {
        let n = data.len();
        let (indices, data): (Vec<usize>, Vec<R>) = data
            .into_iter()
            .enumerate()
            .filter_map(|(i, a)| 
                if a.is_zero() { None } else { Some((i, a)) }
            ).unzip();
        
        CsVec::new(n, indices, data)
    }

    fn is_zero(&self) -> bool {
        self.data().iter().all(|a| a.is_zero())
    }

    fn permute(&self, p: PermView) -> CsVec<R> { 
        let n = self.dim();
        let mut ind = vec![];
        let mut val = vec![];

        for (i, a) in self.iter() { 
            ind.push(p.at(i.index()));
            val.push(a.clone());
        }

        CsVec::new_from_unsorted(n, ind, val).ok().unwrap()
    }

    fn subvec(&self, range: Range<usize>) -> CsVec<R> {
        let (i0, i1) = (range.start, range.end);
        assert!(i0 <= i1 && i1 <= self.dim());

        let mut ind = vec![];
        let mut val = vec![];

        for (i, a) in self.iter() {
            let i = i.index();
            if !a.is_zero() && range.contains(&i) {
                ind.push(i - i0);
                val.push(a.clone());
            }
        }

        CsVec::new(i1 - i0, ind, val)
    }

    fn divide2(self, r: usize) -> (CsVec<R>, CsVec<R>) { 
        let n = self.dim();
        assert!(r <= n);

        let mut ind1: Vec<usize> = vec![];
        let mut val1: Vec<R> = vec![];

        let mut ind2: Vec<usize> = vec![];
        let mut val2: Vec<R> = vec![];

        for (i, a) in self.iter() {
            let i = i.index();
            if i < r { 
                ind1.push(i);
                val1.push(a.clone());
            } else { 
                ind2.push(i - r);
                val2.push(a.clone());
            }
        }

        let v1 = CsVec::new(r, ind1, val1);
        let v2 = CsVec::new(n -r, ind2, val2);
        
        (v1, v2)
    }
}

#[cfg(test)]
mod tests {
    use sprs::PermOwned;
    use super::*;

    #[test]
    fn divide4() { 
        let a = CsMat::csc_from_vec((6, 9), vec![
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 2, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 3, 0, 0, 0, 0,
            0, 1, 0, 1, 0, 0, 1, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 1
        ]);

        let [p,q,r,s] = a.divide4(3, 5);
        assert_eq!(p, CsMat::csc_from_vec((3, 5), vec![
            1, 0, 1, 0, 0,
            0, 1, 1, 1, 0,
            0, 0, 1, 1, 0,
        ]));
        assert_eq!(q, CsMat::csc_from_vec((3, 4), vec![
            1, 1, 0, 1,
            1, 0, 2, 0,
            0, 0, 1, 1,
        ]));
        assert_eq!(r, CsMat::csc_from_vec((3, 5), vec![
            0, 1, 1, 0, 3,
            0, 1, 0, 1, 0,
            1, 0, 1, 0, 1,
        ]));
        assert_eq!(s, CsMat::csc_from_vec((3, 4), vec![
            0, 0, 0, 0,
            0, 1, 0, 1,
            1, 0, 1, 1,
        ]));
    }

    #[test]
    fn combine4() { 
        let p = CsMat::csc_from_vec((3, 5), vec![
            1, 0, 1, 0, 0,
            0, 1, 1, 1, 0,
            0, 0, 1, 1, 0,
        ]);
        let q = CsMat::csc_from_vec((3, 4), vec![
            1, 1, 0, 1,
            1, 0, 2, 0,
            0, 0, 1, 1,
        ]);
        let r = CsMat::csc_from_vec((3, 5), vec![
            0, 1, 1, 0, 3,
            0, 1, 0, 1, 0,
            1, 0, 1, 0, 1,
        ]);
        let s = CsMat::csc_from_vec((3, 4), vec![
            0, 0, 0, 0,
            0, 1, 0, 1,
            1, 0, 1, 1,
        ]);

        assert_eq!(CsMat::combine4([&p,&q,&r,&s]), CsMat::csc_from_vec((6, 9), vec![
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 2, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 3, 0, 0, 0, 0,
            0, 1, 0, 1, 0, 0, 1, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 1
        ]));
    }

    #[test]
    fn permute_vec() { 
        let p = PermOwned::new(vec![1,3,0,2]);
        let v = CsVec::new(4, vec![0,1,2,3], vec![0,1,2,3]);
        let w = v.permute(p.view());
        assert_eq!(w, CsVec::new(4, vec![0,1,2,3], vec![2,0,3,1]))
    }

    #[test]
    fn subvec() {
        let v = CsVec::new(10, (0..10).collect(), (0..10).collect());
        let w = v.subvec(3..7);
        assert_eq!(w, CsVec::new(4, vec![0,1,2,3], vec![3,4,5,6]))
    }

    #[test]
    fn divide2_vec() {
        let v = CsVec::new(10, (0..10).collect(), (0..10).collect());
        let (v1, v2) = v.divide2(3);
        assert_eq!(v1, CsVec::new(3, vec![0,1,2], vec![0,1,2]));
        assert_eq!(v2, CsVec::new(7, vec![0,1,2,3,4,5,6], vec![3,4,5,6,7,8,9]));
    }
}