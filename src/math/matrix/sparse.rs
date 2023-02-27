use std::ops::{Add, AddAssign, Neg, Deref, Range, Sub, SubAssign, Mul, MulAssign};
use std::iter::zip;
use std::fmt::Display;
use num_traits::Zero;
use rand::Rng;
use sprs::{SpIndex, TriMat, CsMat, PermView, CsVec, CsVecBase, CsVecView};
use auto_impl_ops::auto_ops;
use crate::math::traits::{Ring, RingOps, AddMonOps, AddGrpOps, MonOps};

use super::DnsMat;

pub use super::dense::MatType;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    cs_mat: CsMat<R>
}

impl<R> SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn generate<F>(shape: (usize, usize), f: F) -> Self
    where F: FnOnce(&mut (dyn FnMut(usize, usize, R) + Send + Sync)) { 
        let mut t = TriMat::new(shape);
        f( &mut |i, j, a| { 
            if !a.is_zero() { 
                t.add_triplet(i, j, a)
            }
        });
        let cs_mat = t.to_csc();
        Self::from(cs_mat)
    }

    pub fn from_vec(shape: (usize, usize), grid: Vec<R>) -> Self { 
        let n = shape.1;
        Self::generate(shape, |init| { 
            grid.into_iter().enumerate().for_each(|(k, a)| { 
                let (i, j) = (k / n, k % n);
                init(i, j, a)
            })
        })
    }

    pub fn cs_mat(&self) -> &CsMat<R> { 
        &self.cs_mat
    }

    pub fn cs_mat_into(self) -> CsMat<R> { 
        self.cs_mat
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, usize, &R)> { 
        self.cs_mat.iter().map(|(a, (i, j))| {
            (i, j, a)
        })
    }

    pub fn to_dense(&self) -> DnsMat<R> { 
        self.into()
    }
}

impl<R> From<CsMat<R>> for SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn from(cs_mat: CsMat<R>) -> Self {
        assert!(cs_mat.is_csc());
        Self { cs_mat }
    }
}

impl<R> From<&DnsMat<R>> for SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn from(a: &DnsMat<R>) -> Self {
        let n = a.cols();
        SpMat::generate(a.shape(), |set| { 
            for (k, a) in a.array().iter().enumerate() {
                if a.is_zero() { continue }
                let (i, j) = (k / n, k % n);
                set(i, j, a.clone());
            }
        })
    }
}

impl<R> Display for SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "") // TODO
    }
}

impl<R> Default for SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn default() -> Self {
        todo!()
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

    fn shape(&self) -> (usize, usize) {
        self.cs_mat.shape()
    }

    fn zero(shape: (usize, usize)) -> Self {
        Self::from(CsMat::zero(shape).to_csc())
    }

    fn is_zero(&self) -> bool {
        self.cs_mat.data().iter().all(|a| a.is_zero())
    }

    fn id(n: usize) -> Self { 
        let indptr = (0..=n).collect();
        let indices = (0..n).collect();
        let data = vec![R::one(); n];
        let cs_mat = CsMat::new_csc((n, n), indptr, indices, data);
        Self::from(cs_mat)
    }

    fn is_id(&self) -> bool {
        self.is_square() && self.cs_mat.into_iter().all(|(a, (i, j))| 
            (i == j && a.is_one()) || (i != j && a.is_zero())
        )
    }
}

impl<R> SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn rand(shape: (usize, usize), density: f64) -> SpMat<R> { 
        let (m, n) = shape;
        let mut rng = rand::thread_rng();

        Self::generate(shape, |set| { 
            for i in 0..m { 
                for j in 0..n { 
                    if rng.gen::<f64>() < density { 
                        set(i, j, R::one());
                    }
                }
            }
        })
    }

    pub fn transpose(&self) -> Self { 
        Self::from(self.cs_mat.transpose_view().to_csc())
    }

    pub fn permute(&self, p: PermView, q: PermView) -> Self { 
        Self::generate(self.shape(), |set| { 
            for (i, j, a) in self.iter() { 
                if a.is_zero() { continue }
                let (i, j) = (i.index(), j.index());
                set(p.at(i), q.at(j), a.clone());
            }
        })
    }

    pub fn permute_rows(&self, p: PermView) -> Self { 
        let id = PermView::identity(self.cols());
        self.permute(p, id)
    }
    
    pub fn permute_cols(&self, q: PermView) -> Self { 
        let id = PermView::identity(self.rows());
        self.permute(id, q)
    }

    pub fn submatrix(&self, rows: Range<usize>, cols: Range<usize>) -> Self { 
        let (i0, i1) = (rows.start, rows.end);
        let (j0, j1) = (cols.start, cols.end);

        assert!(i0 <= i1 && i1 <= self.rows());
        assert!(j0 <= j1 && j1 <= self.cols());

        Self::generate((i1 - i0, j1 - j0), |set| { 
            for (i, j, a) in self.iter() { 
                let (i, j) = (i.index(), j.index());
                if !a.is_zero() && rows.contains(&i) && cols.contains(&j) {
                    set(i - i0, j - j0, a.clone());
                }
            }
        })
    }

    pub fn split_blocks(&self, k: usize, l: usize) -> [SpMat<R>; 4] {
        let (m, n) = self.shape();
        assert!(k <= m);
        assert!(l <= n);

        let mut trips = [
            TriMat::new((k, l)),
            TriMat::new((k, n - l)),
            TriMat::new((m - k, l)),
            TriMat::new((m - k, n - l))
        ];
        
        for (i, j, a) in self.iter() { 
            if a.is_zero() { continue }

            let (a, i, j) = (a.clone(), i.index(), j.index());
            match ((0..k).contains(&i), (0..l).contains(&j)) { 
                (true , true ) => trips[0].add_triplet(i,     j,     a),
                (true , false) => trips[1].add_triplet(i,     j - l, a),
                (false, true ) => trips[2].add_triplet(i - k, j,     a),
                (false, false) => trips[3].add_triplet(i - k, j - l, a),
            }
        }

        trips.map(|t| SpMat::from(t.to_csc()))
    }

    pub fn split_hor(&self, k: usize) -> (Self, Self) { 
        let [a, b, _, _] = self.split_blocks(self.rows(), k);
        (a, b)
    }

    pub fn split_ver(&self, k: usize) -> (Self, Self) { 
        let [a, _, c, _] = self.split_blocks(k, self.cols());
        (a, c)
    }

    pub fn combine_blocks(blocks: [&SpMat<R>; 4]) -> SpMat<R> {
        let [a, b, c, d] = blocks;

        assert_eq!(a.rows(), b.rows());
        assert_eq!(c.rows(), d.rows());
        assert_eq!(a.cols(), c.cols());
        assert_eq!(b.cols(), d.cols());

        let (m, n) = (a.rows() + c.rows(), a.cols() + b.cols());
        let (k, l) = a.shape();
        let targets = zip(
            [a, b, c, d], 
            [(0,0), (0,l), (k,0), (k,l)]
        );

        Self::generate((m, n), |set| { 
            for (x, (di, dj)) in targets { 
                for (i, j, r) in x.iter() {
                    set(i + di, j + dj, r.clone());
                }
            }    
        })
    }

    pub fn concat(&self, b: &Self) -> Self { 
        let zero = |m, n| SpMat::<R>::zero((m, n));
        Self::combine_blocks([
            self, 
            b, 
            &zero(0, self.cols()), 
            &zero(0, b.cols())
        ])
    }

    pub fn stack(&self, b: &Self) -> Self { 
        let zero = |m, n| SpMat::<R>::zero((m, n));
        Self::combine_blocks([
            self, 
            &zero(self.rows(), 0), 
            b, 
            &zero(b.rows(), 0)
        ])
    }

    pub fn col(&self, j: usize) -> CsVec<R> { 
        self.col_view(j).to_owned()
    }

    pub fn col_view(&self, j: usize) -> CsVecView<R> { 
        self.cs_mat.outer_view(j).unwrap()
    }
}

impl<R> Mul<&CsVec<R>> for &SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    type Output = CsVec<R>;

    fn mul(self, rhs: &CsVec<R>) -> Self::Output {
        self.cs_mat() * rhs
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn init() { 
        let a = SpMat::generate((2, 2), |set| { 
            set(0, 0, 1);
            set(0, 1, 2);
            set(1, 0, 3);
            set(1, 1, 4);
        });
        assert_eq!(&a.cs_mat, &CsMat::new_csc((2, 2), vec![0, 2, 4], vec![0, 1, 0, 1], vec![1, 3, 2, 4]));
    }

    #[test]
    fn from_grid() { 
        let a = SpMat::from_vec((2, 2), vec![1,2,3,4]);
        assert_eq!(&a.cs_mat, &CsMat::new_csc((2, 2), vec![0, 2, 4], vec![0, 1, 0, 1], vec![1, 3, 2, 4]));
    }

    #[test]
    fn to_dense() { 
        let a = SpMat::generate((2, 2), |set| { 
            set(0, 0, 1);
            set(0, 1, 2);
            set(1, 0, 3);
            set(1, 1, 4);
        });
        assert_eq!(a.to_dense(), DnsMat::from(ndarray::array![[1, 2], [3, 4]]));
    }
}

// -- old code -- //

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
mod tests_old {
    use sprs::PermOwned;
    use super::*;

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