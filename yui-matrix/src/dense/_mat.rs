use std::ops::{Add, Neg, Sub, Mul, Index, IndexMut, AddAssign, SubAssign, MulAssign, Range};
use std::fmt::Debug;
use nalgebra::{ClosedSub, ClosedMul};
use nalgebra_sparse::CscMatrix;
use nalgebra_sparse::na::{Scalar, ClosedAdd, DMatrix};
use delegate::delegate;
use derive_more::Display;
use auto_impl_ops::auto_ops;
use num_traits::{Zero, One};
use crate::sparse::SpMat;

pub trait MatType {
    fn shape(&self) -> (usize, usize);
    fn rows(&self) -> usize { self.shape().0 }
    fn cols(&self) -> usize { self.shape().1 }
    fn is_square(&self) -> bool { 
        let (m, n) = self.shape();
        m == n
    }
}

#[derive(Clone, Debug, Display, PartialEq, Eq)]
pub struct Mat<R> {
    inner: DMatrix<R>
}

impl<R> MatType for Mat<R> {
    fn shape(&self) -> (usize, usize) {
        (self.inner.nrows(), self.inner.ncols())
    }
}

impl<R> Mat<R> {
    pub fn inner(&self) -> &DMatrix<R> {
        &self.inner
    }

    pub fn inner_mut(&mut self) -> &mut DMatrix<R> {
        &mut self.inner
    }

    pub fn into_inner(self) -> DMatrix<R> {
        self.inner
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, usize, &R)> { 
        let m = self.rows();
        self.inner.iter().enumerate().map(move |(i, a)| 
            (i % m, i / m, a)
        )
    }
}

impl<R> Mat<R>
where R: Scalar {
    pub fn from_data<I>(shape: (usize, usize), data: I) -> Self
    where I: IntoIterator<Item = R> { 
        DMatrix::from_row_iterator(shape.0, shape.1, data).into()
    }

    pub fn zero(shape: (usize, usize)) -> Self
    where R: Zero { 
        let inner = DMatrix::zeros(shape.0, shape.1);
        Self::from(inner)
    }

    pub fn is_zero(&self) -> bool
    where R: Zero { 
        self.iter().all(|e| e.2.is_zero())
    }

    pub fn id(size: usize) -> Self
    where R: Zero + One { 
        let inner = DMatrix::identity(size, size);
        Self::from(inner)
    }

    pub fn is_id(&self) -> bool
    where R: Zero + One { 
        self.is_square() && self.iter().all(|(i, j, a)| 
            i == j && a.is_one() || 
            i != j && a.is_zero()
        )
    }

    pub fn diag<I>(shape: (usize, usize), entries: I) -> Self
    where R: Zero, I: IntoIterator<Item = R> {
        let mut mat = Self::zero(shape);
        for (i, a) in entries.into_iter().enumerate() {
            mat[(i, i)] = a;
        }
        mat
    }

    pub fn is_diag(&self) -> bool
    where R: Zero { 
        self.iter().all(|(i, j, a)| 
            i == j || a.is_zero()
        )
    }

    pub fn submat(&self, rows: Range<usize>, cols: Range<usize>) -> Mat<R> { 
        let (i0, i1) = (rows.start, rows.end);
        let (j0, j1) = (cols.start, cols.end);

        assert!(i0 <= i1 && i1 <= self.rows());
        assert!(j0 <= j1 && j1 <= self.cols());

        let slice = self.inner.view((i0, j0), (i1 - i0, j1 - j0));
        Self::from(slice.clone_owned())
    }

    pub fn submat_rows(&self, rows: Range<usize>) -> Mat<R> { 
        let n = self.cols();
        self.submat(rows, 0 .. n)
    }

    pub fn submat_cols(&self, cols: Range<usize>) -> Mat<R> { 
        let m = self.rows();
        self.submat(0 .. m, cols)
    }

    pub fn to_sparse(self) -> SpMat<R>
    where R: Zero + ClosedAdd { 
        self.into()
    }
}

impl<R> From<DMatrix<R>> for Mat<R> {
    fn from(inner: DMatrix<R>) -> Self {
        Self { inner }
    }
}

// TODO move to SpMat.
impl<R> From<Mat<R>> for SpMat<R>
where R: Scalar + Zero {
    fn from(value: Mat<R>) -> Self {
        let csc = CscMatrix::from(&value.inner);
        Self::from(csc)
    }
}

impl<R> From<SpMat<R>> for Mat<R>
where R: Scalar + Zero + ClosedAdd {
    fn from(value: SpMat<R>) -> Self {
        let inner = DMatrix::from(value.inner());
        Self::from(inner)
    }
}
 
impl<R> Index<(usize, usize)> for Mat<R> {
    type Output = R;
    delegate! { 
        to self.inner { 
            fn index(&self, index: (usize, usize)) -> &R;
        }
    }
}

impl<R> IndexMut<(usize, usize)> for Mat<R> {
    delegate! { 
        to self.inner { 
            fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output;
        }
    }
}

impl<R> Default for Mat<R>
where R: Scalar + Zero {
    fn default() -> Self {
        Self::zero((0, 0))
    }
}

// TODO move to SpMat.
impl<R> SpMat<R> { 
    pub fn to_dense(self) -> Mat<R>
    where R: Scalar + Zero + ClosedAdd { 
        self.into()
    }
}

impl<R> Neg for Mat<R>
where R: Scalar + Neg<Output = R> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Mat::from(-self.inner)
    }
}

impl<R> Neg for &Mat<R>
where R: Scalar + Neg<Output = R> {
    type Output = Mat<R>;
    fn neg(self) -> Self::Output {
        Mat::from(-&self.inner)
    }
}

#[auto_ops]
impl<R> AddAssign<&Mat<R>> for Mat<R>
where R: Scalar + ClosedAdd {
    fn add_assign(&mut self, rhs: &Self) {
        self.inner += &rhs.inner;
    }
}

#[auto_ops]
impl<R> SubAssign<&Mat<R>> for Mat<R>
where R: Scalar + ClosedSub {
    fn sub_assign(&mut self, rhs: &Self) {
        self.inner -= &rhs.inner
    }
}

#[auto_ops]
impl<'a, 'b, R> Mul<&'b Mat<R>> for &'a Mat<R>
where R: Scalar + Zero + One + ClosedAdd + ClosedMul {
    type Output = Mat<R>;
    fn mul(self, rhs: &'b Mat<R>) -> Self::Output {
        let prod = &self.inner * &rhs.inner;
        Mat::from(prod)
    }
}

impl<R> Mat<R>
where R: Scalar { 
    pub fn swap_rows(&mut self, i: usize, j: usize) {
        self.inner.swap_rows(i, j);
    }

    pub fn swap_cols(&mut self, i: usize, j: usize) {
        self.inner.swap_columns(i, j);
    }

    pub fn mul_row(&mut self, i: usize, r: &R)
    where R: ClosedMul {
        self.inner.row_mut(i).mul_assign(r.clone())
    }

    pub fn mul_col(&mut self, j: usize, r: &R)
    where R: ClosedMul {
        self.inner.column_mut(j).mul_assign(r.clone())
    }

    pub fn add_row_to(&mut self, i: usize, j: usize, r: &R)
    where R: ClosedAdd + ClosedMul { 
        let row = self.inner.row(i).mul(r.clone());
        self.inner.row_mut(j).add_assign(row)
    }

    pub fn add_col_to(&mut self, i: usize, j: usize, r: &R)
    where R: ClosedAdd + ClosedMul {  
        let col = self.inner.column(i).mul(r.clone());
        self.inner.column_mut(j).add_assign(col)
    }

    // Multiply [a, b; c, d] from left. 
    pub fn left_elementary(&mut self, comps: [&R; 4], i: usize, j: usize)
    where R: ClosedAdd + ClosedMul { 
        let [a, b, c, d] = comps.map(Clone::clone);

        let r_i = self.inner.row(i);
        let r_j = self.inner.row(j);
        
        let s_i = &r_i * a + &r_j * b;
        let s_j = &r_i * c + &r_j * d;

        self.inner.set_row(i, &s_i);
        self.inner.set_row(j, &s_j);
    }

    // Multiply [a, c; b, d] from right. 
    pub fn right_elementary(&mut self, comps: [&R; 4], i: usize, j: usize) 
    where R: ClosedAdd + ClosedMul { 
        let [a, b, c, d] = comps.map(Clone::clone);

        let r_i = self.inner.column(i);
        let r_j = self.inner.column(j);
        
        let s_i = &r_i * a + &r_j * b;
        let s_j = &r_i * c + &r_j * d;

        self.inner.set_column(i, &s_i);
        self.inner.set_column(j, &s_j);
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn init() { 
        let a = Mat::from_data((2, 3), [1,2,3,4,5,6]);

        assert_eq!(a.rows(), 2);
        assert_eq!(a.cols(), 3);
        assert_eq!(a.into_inner(), DMatrix::from_row_slice(2, 3, &[1,2,3,4,5,6]));
    }

    #[test]
    fn eq() {
        let a = Mat::from_data((2, 3), [1,2,3,4,5,6]);
        let b = Mat::from_data((2, 3), [1,2,0,4,5,6]);
        let c = Mat::from_data((3, 2), [1,2,3,4,5,6]);

        assert_eq!(a, a);
        assert_ne!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn square() {
        let a: Mat<i32> = Mat::zero((3, 3));
        assert!(a.is_square());

        let a: Mat<i32> = Mat::zero((3, 2));
        assert!(!a.is_square());
    }

    #[test]
    fn zero() {
        let a: Mat<i32> = Mat::zero((3, 2));
        assert!(a.is_zero());

        let a = Mat::from_data((2, 3), [1,2,3,4,5,6]);
        assert!(!a.is_zero());
    }

    #[test]
    fn id() {
        let a: Mat<i32> = Mat::id(3);
        assert!(a.is_id());

        let a = Mat::from_data((2, 2), [1,2,3,4]);
        assert!(!a.is_id());

        let a = Mat::from_data((2, 3), [1,0,0,0,1,0]);
        assert!(!a.is_id());
    }

    #[test]
    fn swap_rows() { 
        let mut a = Mat::from_data((3, 4), 1..=12);
        a.swap_rows(0, 1);
        assert_eq!(a, Mat::from_data((3, 4), [5,6,7,8,1,2,3,4,9,10,11,12]));
    }

    #[test]
    fn swap_cols() { 
        let mut a = Mat::from_data((3, 4), 1..=12);
        a.swap_cols(0, 1);
        assert_eq!(a, Mat::from_data((3, 4), [2,1,3,4,6,5,7,8,10,9,11,12]));
    }

    #[test]
    fn mul_row() { 
        let mut a = Mat::from_data((3, 3), 1..=9);
        a.mul_row(1, &10);
        assert_eq!(a, Mat::from_data((3, 3), [1,2,3,40,50,60,7,8,9]));
    }

    #[test]
    fn mul_col() { 
        let mut a = Mat::from_data((3, 3), 1..=9);
        a.mul_col(1, &10);
        assert_eq!(a, Mat::from_data((3, 3), [1,20,3,4,50,6,7,80,9]));
    }

    #[test]
    fn add_row_to() { 
        let mut a = Mat::from_data((3, 3), 1..=9);
        a.add_row_to(0, 1, &10);
        assert_eq!(a, Mat::from_data((3, 3), [1,2,3,14,25,36,7,8,9]));
    }

    #[test]
    fn add_col_to() { 
        let mut a = Mat::from_data((3, 3), 1..=9);
        a.add_col_to(0, 1, &10);
        assert_eq!(a, Mat::from_data((3, 3), [1,12,3,4,45,6,7,78,9]));
    }

    #[test]
    fn add() { 
        let a = Mat::from_data((3, 2), [1,2,3,4,5,6]);
        let b = Mat::from_data((3, 2), [8,2,4,0,2,1]);
        let c = a + b;
        assert_eq!(c, Mat::from_data((3, 2), [9,4,7,4,7,7]));
    }

    #[test]
    fn sub() { 
        let a = Mat::from_data((3, 2), [1,2,3,4,5,6]);
        let b = Mat::from_data((3, 2), [8,2,4,0,2,1]);
        let c = a - b;
        assert_eq!(c, Mat::from_data((3, 2), [-7,0,-1,4,3,5]));
    }

    #[test]
    fn neg() { 
        let a = Mat::from_data((3, 2), [1,2,3,4,5,6]);
        assert_eq!(-a, Mat::from_data((3, 2), [-1,-2,-3,-4,-5,-6]));
    }

    #[test]
    fn mul() { 
        let a = Mat::from_data((2, 3), [1,2,3,4,5,6]);
        let b = Mat::from_data((3, 2), [1,2,1,-1,0,2]);
        let c = a * b;
        assert_eq!(c, Mat::from_data((2, 2), [3,6,9,15]));
    }

    #[test]
    fn to_sparse() { 
        let dns = Mat::from_data((2, 3), [1,2,3,4,5,6]);
        let sps = dns.to_sparse();
        assert_eq!(sps, SpMat::from_dense_data((2, 3), [1,2,3,4,5,6]));
    }

    #[test]
    fn from_sparse() { 
        let sps = SpMat::from_dense_data((2, 3), [1,2,3,4,5,6]);
        let dns = Mat::from(sps);
        assert_eq!(dns, Mat::from_data((2, 3), [1,2,3,4,5,6]));
    }

    #[test]
    fn submat() { 
        let a = Mat::from_data((3, 4), [
            1, 2, 3, 7,
            4, 5, 6, 8,
            9,10,11,12           
        ]);
        let b = a.submat(1..3, 2..4);
        assert_eq!(b, Mat::from_data((2, 2), [
             6, 8,
            11,12           
        ]));
    }
}