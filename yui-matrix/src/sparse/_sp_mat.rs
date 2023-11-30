#![allow(unused)]
use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, MulAssign, Range};
use std::iter::zip;
use std::fmt::{Display, Debug};
use std::sync::Mutex;
use delegate::delegate;
use derive_more::Display;
use nalgebra_sparse::na::{Scalar, ClosedAdd, ClosedSub, ClosedMul};
use nalgebra_sparse::{CscMatrix, CooMatrix};
use itertools::Itertools;
use num_traits::{Zero, One};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use auto_impl_ops::auto_ops;
use sprs::PermView;
use yui::{Ring, RingOps, AddMonOps, AddGrpOps, MonOps, AddGrp};
use crate::dense::*;
use super::sp_vec::SpVec;

#[derive(Clone, PartialEq, Eq)]
pub struct SpMat<R> { 
    inner: CscMatrix<R>
}

impl<R> SpMat<R> { 
    pub fn view(&self) -> SpMatView<R> { 
        SpMatView::from(self)
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, usize, &R)> { 
        self.inner.triplet_iter()
    }

    pub fn iter_nz(&self) -> impl Iterator<Item = (usize, usize, &R)>
    where R: Zero { 
        self.iter().filter(|e| !e.2.is_zero())
    }

    pub fn col(&self, j: usize) -> SpVec<R>
    where R: Clone + Zero { 
        let col = self.inner.col(j);
        let iter = Iterator::zip(col.row_indices().iter().cloned(), col.values().iter().cloned());
        SpVec::from_entries(self.rows(), iter)
    }

    pub fn transpose(&self) -> Self
    where R: Scalar { 
        self.inner.transpose().into()
    }

    pub fn disassemble(self) -> (Vec<usize>, Vec<usize>, Vec<R>) { 
        self.inner.disassemble()
    }
}

impl<R> MatType for SpMat<R> {
    fn shape(&self) -> (usize, usize) {
        (self.inner.nrows(), self.inner.ncols())
    }
}

impl<R> From<CscMatrix<R>> for SpMat<R> {
    fn from(inner: CscMatrix<R>) -> Self {
        Self { inner }
    }
}

impl<R> From<CooMatrix<R>> for SpMat<R>
where R: Scalar + Zero + ClosedAdd {
    fn from(coo: CooMatrix<R>) -> Self {
        let csc = CscMatrix::from(&coo);
        Self::from(csc)
    }
}

impl<R> From<Mat<R>> for SpMat<R>
where R: Scalar + Zero + ClosedAdd {
    fn from(a: Mat<R>) -> Self {
        let n = a.cols();
        let entries = a.array().into_iter().enumerate().filter_map(|(k, a)| {
            if !a.is_zero() { 
                Some((k / n, k % n, a.clone()))
            } else {
                None 
            }
        });
        SpMat::from_entries(a.shape(), entries)
    }
}

impl<R> SpMat<R> 
where R: Scalar + Clone + Zero + ClosedAdd { 
    pub fn from_entries<T>(shape: (usize, usize), entries: T) -> Self
    where T: IntoIterator<Item = (usize, usize, R)> {
        let mut coo = CooMatrix::new(shape.0, shape.1);
        for (i, j, a) in entries { 
            if a.is_zero() { 
                continue;
            }
            coo.push(i, j, a)
        }
        Self::from(coo)
    }

    pub fn from_par_entries<T>(shape: (usize, usize), entries: T) -> Self
    where 
        R: Send + Sync,
        T: IntoParallelIterator<Item = (usize, usize, R)>
    {
        let t = Mutex::new(CooMatrix::new(shape.0, shape.1));
        entries.into_par_iter().for_each(|(i, j, a)| { 
            if a.is_zero() { 
                return;
            }
            t.lock().unwrap().push(i, j, a)
        });
        let coo = t.into_inner().unwrap();
        Self::from(coo)
    }

    pub fn from_dense_data<I>(shape: (usize, usize), data: I) -> Self
    where I: IntoIterator<Item = R> { 
        let n = shape.1;
        Self::from_entries(
            shape, 
            data.into_iter().enumerate().map(|(k, a)| { 
                let (i, j) = (k / n, k % n);
                (i, j, a)
            })
        )
    }
}

// impl<R> IntoIterator for SpMat<R>
// where R: Clone + Zero {
//     type Item = (usize, usize, R);
//     type IntoIter = std::vec::IntoIter<Self::Item>;

//     fn into_iter(self) -> Self::IntoIter {
//         let
//         // MEMO improve this
//         self.iter().map(|(i, j, a)| (i, j, a.clone())).collect_vec().into_iter()
//     }
// }

// impl<R> SpMat<R>
// where R: Clone {
//     pub fn col_vec(&self, j: usize) -> SpVec<R> { 
//         let cs_vec = self.col_view(j).to_owned();
//         SpVec::from(cs_vec)
//     }
// }

impl<R> Default for SpMat<R> {
    fn default() -> Self {
        Self::zero((0, 0))
    }
}

impl<R> SpMat<R> { 
    pub fn zero(shape: (usize, usize)) -> Self {
        let csc = CscMatrix::zeros(shape.0, shape.1);
        Self::from(csc)
    }

    pub fn is_zero(&self) -> bool
    where R: Zero {
        self.inner.values().iter().all(|a| a.is_zero())
    }

    pub fn id(n: usize) -> Self
    where R: Scalar + One { 
        let csc = CscMatrix::identity(n);
        Self::from(csc)
    }

    pub fn is_id(&self) -> bool
    where R: Scalar + One + Zero {
        self.is_square() && self.iter().all(|(i, j, a)| 
            (i == j && a.is_one()) || (i != j && a.is_zero())
        )
    }
}

impl<R> Neg for SpMat<R>
where R: Scalar + Neg<Output = R> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self::from(-self.inner)
    }
}

impl<R> Neg for &SpMat<R>
where R: Scalar + Neg<Output = R> {
    type Output = SpMat<R>;
    fn neg(self) -> Self::Output {
        SpMat::from(-&self.inner)
    }
}

// see: nalgebra_sparse::ops::impl_std_ops.
macro_rules! impl_binop {
    ($trait:ident, $method:ident) => {
        #[auto_ops]
        impl<'a, 'b, R> $trait<&'b SpMat<R>> for &'a SpMat<R>
        where R: Scalar + ClosedAdd + ClosedSub + ClosedMul + Zero + One + Neg<Output = R> {
            type Output = SpMat<R>;
            fn $method(self, rhs: &'b SpMat<R>) -> Self::Output {
                let res = (&self.inner).$method(&rhs.inner);
                SpMat::from(res)
            }
        }
    };
}

impl_binop!(Add, add);
impl_binop!(Sub, sub);
impl_binop!(Mul, mul);

macro_rules! impl_ops {
    ($trait:ident) => {
        impl<R> $trait<SpMat<R>> for SpMat<R>
        where R: Scalar + ClosedAdd + ClosedSub + ClosedMul + Zero + One + Neg<Output = R> {}

        impl<R> $trait<SpMat<R>> for &SpMat<R>
        where R:Scalar + ClosedAdd + ClosedSub + ClosedMul + Zero + One + Neg<Output = R> {}
    };
}

impl_ops!(AddMonOps);
impl_ops!(AddGrpOps);
impl_ops!(MonOps);
impl_ops!(RingOps);

impl<R> SpMat<R>
where R: Scalar + Clone + Zero + ClosedAdd { 
    pub fn permute(&self, p: PermView, q: PermView) -> SpMat<R> { 
        self.view().permute(p, q).collect()
    }

    pub fn permute_rows(&self, p: PermView) -> SpMat<R> { 
        self.view().permute_rows(p).collect()
    }
    
    pub fn permute_cols(&self, q: PermView) -> SpMat<R> { 
        self.view().permute_cols(q).collect()
    }

    pub fn submat(&self, rows: Range<usize>, cols: Range<usize>) -> SpMat<R> { 
        self.view().submat(rows, cols).collect()
    }

    pub fn submat_rows(&self, rows: Range<usize>) -> SpMat<R> { 
        self.view().submat_rows(rows).collect()
    }

    pub fn submat_cols(&self, cols: Range<usize>) -> SpMat<R> { 
        self.view().submat_cols(cols).collect()
    }

    pub fn combine_blocks(blocks: [&SpMat<R>; 4]) -> SpMat<R> {
        let [a, b, c, d] = blocks;

        assert_eq!(a.rows(), b.rows());
        assert_eq!(c.rows(), d.rows());
        assert_eq!(a.cols(), c.cols());
        assert_eq!(b.cols(), d.cols());

        let (m, n) = (a.rows() + c.rows(), a.cols() + b.cols());
        let (k, l) = a.shape();

        let entries = zip(
            [a, b, c, d], 
            [(0,0), (0,l), (k,0), (k,l)]
        ).flat_map(|(x, (di, dj))| 
            x.iter().map(move |(i, j, r)|
                (i + di, j + dj, r.clone())
            )
        );

        Self::from_entries((m, n), entries)
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

    // pub fn to_dense(self) -> Mat<R> { 
    //     self.into()
    // }
}

impl<R> Display for SpMat<R>
where R: Display + Debug {
    delegate! { to self.inner { 
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result;
    }}
}

impl<R> Debug for SpMat<R>
where R: Display + Debug {
    delegate! { to self.inner { 
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result;
    }}
}

pub struct SpMatView<'a, 'b, R>  {
    target: &'a SpMat<R>,
    shape: (usize, usize),
    trans: Box<dyn Fn(usize, usize) -> Option<(usize, usize)> + 'b>
}

impl<'a, R> From<&'a SpMat<R>> for SpMatView<'a, '_, R> {
    fn from(target: &'a SpMat<R>) -> Self {
        Self::new(target, target.shape(), |i, j| Some((i, j)))
    }
}

impl<'a, 'b, R> SpMatView<'a, 'b, R> {
    fn new<F>(target: &'a SpMat<R>, shape: (usize, usize), trans: F) -> Self
    where F: Fn(usize, usize) -> Option<(usize, usize)> + 'b {
        Self { target, shape, trans: Box::new(trans) }
    }

    pub fn shape(&self) -> (usize, usize) { 
        self.shape
    }

    pub fn nrows(&self) -> usize { 
        self.shape.0
    }

    pub fn ncols(&self) -> usize { 
        self.shape.1
    }

    pub fn transpose(self) -> SpMatView<'a, 'b, R> { 
        SpMatView::new(self.target, (self.shape.1, self.shape.0), move |i, j| (self.trans)(j, i))
    }

    pub fn permute(self, p: PermView<'b>, q: PermView<'b>) -> SpMatView<'a, 'b, R> { 
        assert_eq!(self.nrows(), p.dim());
        assert_eq!(self.ncols(), q.dim());
        let trans = self.trans;
        SpMatView::new(self.target, self.shape, move |i, j| (trans)(p.at(i), q.at(j)))
    }

    pub fn permute_rows(self, p: PermView<'b>) -> SpMatView<'a, 'b, R> { 
        let id = PermView::identity(self.ncols());
        self.permute(p, id)
    }
    
    pub fn permute_cols(self, q: PermView<'b>) -> SpMatView<'a, 'b, R> { 
        let id = PermView::identity(self.nrows());
        self.permute(id, q)
    }

    pub fn submat(self, rows: Range<usize>, cols: Range<usize>) -> SpMatView<'a, 'b, R> { 
        let (i0, i1) = (rows.start, rows.end);
        let (j0, j1) = (cols.start, cols.end);

        assert!(i0 <= i1 && i1 <= self.nrows());
        assert!(j0 <= j1 && j1 <= self.ncols());

        SpMatView::new(self.target, (i1 - i0, j1 - j0), move |i, j| { 
            let (i, j) = (self.trans)(i, j)?;
            if rows.contains(&i) && cols.contains(&j) {
                Some((i - i0, j - j0))
            } else { 
                None
            }
        })
    }

    pub fn submat_rows(self, rows: Range<usize>) -> SpMatView<'a, 'b, R> { 
        let n = self.ncols();
        self.submat(rows, 0 .. n)
    }

    pub fn submat_cols(self, cols: Range<usize>) -> SpMatView<'a, 'b, R> { 
        let m = self.nrows();
        self.submat(0 .. m, cols)
    }

    pub fn collect(self) -> SpMat<R> 
    where R: Scalar + Clone + Zero + ClosedAdd {
        SpMat::from_entries(self.shape(), self.iter().map(|(i, j, a)| 
            (i, j, a.clone())
        ))
    }
}

impl<'a, 'b, R> SpMatView<'a, 'b, R>
where 'a: 'b, R: Zero {
    pub fn iter(self) -> impl Iterator<Item = (usize, usize, &'a R)> + 'b {
        let trans = self.trans;
        self.target.iter().filter_map(move |(i, j, a)| { 
            (trans)(i, j).map(|(i, j)| (i, j, a))
        })
    }    
}

// impl<R> SpMat<R> where R: Clone + Zero + One { 
//     // row_perm(p) * a == a.permute_rows(p)
//     pub fn from_row_perm(p: PermView) -> Self {
//         let n = p.dim();
//         Self::from_entries((n, n), (0..n).map(|i|
//             (p.at(i), i, R::one())
//         ))
//     }

//     // a * col_perm(p) == a.permute_cols(p)
//     pub fn from_col_perm(p: PermView) -> Self {
//         let n = p.dim();
//         Self::from_entries((n, n), (0..n).map(|i|
//             (i, p.at(i), R::one())
//         ))
//     }
// }

// #[cfg(test)]
// impl<R> SpMat<R>
// where R: Ring, for<'a> &'a R: RingOps<R> { 
//     pub fn rand(shape: (usize, usize), density: f64) -> Self {
//         use cartesian::cartesian;
//         use rand::Rng;
    
//         let (m, n) = shape;
//         let range = cartesian!(0..m, 0..n);
//         let mut rng = rand::thread_rng();
    
//         Self::from_entries(shape, range.filter_map(|(i, j)|
//             if rng.gen::<f64>() < density { 
//                 Some((i, j, R::one()))
//             } else { 
//                 None
//             }
//         ))
//     }
// }

#[cfg(test)]
pub(super) mod tests { 
    use sprs::PermOwned;
    use yui::Ratio;

    use super::*;

    #[test]
    fn init() { 
        let a = SpMat::from_entries((2, 2), [
            (0, 0, 1),
            (0, 1, 2),
            (1, 0, 3),
            (1, 1, 4)
        ]);
        assert_eq!(a.disassemble(), (vec![0, 2, 4], vec![0, 1, 0, 1], vec![1, 3, 2, 4]));
    }

    #[test]
    fn init_ratio() { 
        type R = Ratio<i64>;
        let vals = (0..4).map(|i| R::new(i + 1, 5)).collect_vec();
        let a = SpMat::from_entries((2, 2), [
            (0, 0, vals[0].clone()),
            (0, 1, vals[2].clone()),
            (1, 0, vals[1].clone()),
            (1, 1, vals[3].clone())
        ]);
        assert_eq!(a.disassemble(), (vec![0, 2, 4], vec![0, 1, 0, 1], vals));
    }

    #[test]
    fn from_grid() { 
        let a = SpMat::from_dense_data((2, 2), [1,2,3,4]);
        assert_eq!(a.disassemble(), (vec![0, 2, 4], vec![0, 1, 0, 1], vec![1, 3, 2, 4]));
    }

//     #[test]
//     fn to_dense() { 
//         let a = SpMat::from_entries((2, 2), [
//             (0, 0, 1),
//             (0, 1, 2),
//             (1, 0, 3),
//             (1, 1, 4)
//         ]);
//         assert_eq!(a.to_dense(), Mat::from(ndarray::array![[1, 2], [3, 4]]));
//     }

    #[test]
    fn permute() { 
        let p = PermOwned::new(vec![1,2,3,0]);
        let q = PermOwned::new(vec![3,0,2,1]);
        let a = SpMat::from_dense_data((4,4), 0..16);
        let b = a.permute(p.view(), q.view());
        assert_eq!(b, SpMat::from_dense_data((4,4), vec![
            13, 15, 14, 12,
             1,  3,  2,  0,
             5,  7,  6,  4,
             9, 11, 10,  8,
        ]));
    }

    #[test]
    fn transpose() { 
        let a = SpMat::from_dense_data((3,4), 0..12);
        let b = a.transpose();

        assert_eq!(b, SpMat::from_dense_data((4,3), vec![
            0, 4, 8, 
            1, 5, 9, 
            2, 6, 10, 
            3, 7, 11, 
        ]));
    }

//     #[test]
//     fn row_perm() {
//         let a = SpMat::from_dense_data((3, 4), 0..12);
//         let p = PermOwned::new(vec![2,0,1]);
//         let q = SpMat::from_row_perm(p.view());
//         assert!(q * &a == a.permute_rows(p.view()))
//     }

//     #[test]
//     fn col_perm() {
//         let a = SpMat::from_dense_data((3, 4), 0..12);
//         let p = PermOwned::new(vec![2,0,1,3]);
//         let q = SpMat::from_col_perm(p.view());
//         assert!(&a * q == a.permute_cols(p.view()))
//     }
}