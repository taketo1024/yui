use std::ops::{Add, Neg, Sub, Mul, Index, IndexMut, AddAssign, SubAssign, MulAssign};
use std::cmp::min;
use std::fmt::{Debug, Display};
use ndarray::{Array2, s};
use derive_more::Display;
use auto_impl_ops::auto_ops;
use yui_core::{Ring, RingOps, AddMonOps, AddGrpOps, MonOps};
use crate::sparse::SpMat;

pub trait MatType: Clone + Debug + Display + Default + PartialEq + Eq + RingOps<Self> 
where for<'x> &'x Self: RingOps<Self>
{
    type R;

    fn shape(&self) -> (usize, usize);
    fn rows(&self) -> usize { self.shape().0 }
    fn cols(&self) -> usize { self.shape().1 }

    fn is_square(&self) -> bool { 
        let (m, n) = self.shape();
        m == n
    }

    fn zero(shape: (usize, usize)) -> Self;
    fn is_zero(&self) -> bool;

    fn id(size: usize) -> Self;
    fn is_id(&self) -> bool;
}

#[derive(Clone, Debug, Display, Default, PartialEq, Eq)]
pub struct Mat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    array: Array2<R>
}

impl<R> Mat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    pub fn array(&self) -> &Array2<R> {
        &self.array
    }

    pub fn array_mut(&mut self) -> &mut Array2<R> {
        &mut self.array
    }

    pub fn array_into(self) -> Array2<R> {
        self.array
    }

    pub fn diag(shape: (usize, usize), entries: Vec<R>) -> Self {
        assert!( entries.len() <= min(shape.0, shape.1) );
        let mut mat = Self::zero(shape);
        for (i, a) in entries.into_iter().enumerate() {
            mat[[i, i]] = a;
        }
        mat
    }

    pub fn is_diag(&self) -> bool { 
        self.array.indexed_iter().all(|((i, j), a)| 
            i == j || a.is_zero()
        )
    }

    pub fn to_sparse(&self) -> SpMat<R> { 
        self.into()
    }

    pub fn swap_rows(&mut self, i: usize, j: usize) {
        debug_assert_ne!(i, j);
        debug_assert!(self.is_valid_row_index(i));
        debug_assert!(self.is_valid_row_index(j));

        let (s_i, s_j) = (s![i, ..], s![j, ..]);
        let (row_i, row_j) = self.array.multi_slice_mut((s_i, s_j));
        ndarray::Zip::from(row_i).and(row_j).for_each(std::mem::swap);
    }

    pub fn swap_cols(&mut self, i: usize, j: usize) {
        debug_assert_ne!(i, j);
        debug_assert!(self.is_valid_col_index(i));
        debug_assert!(self.is_valid_col_index(j));

        let (s_i, s_j) = (s![.., i], s![.., j]);
        let (col_i, col_j) = self.array.multi_slice_mut((s_i, s_j));
        ndarray::Zip::from(col_i).and(col_j).for_each(std::mem::swap);
    }

    pub fn mul_row(&mut self, i: usize, r: &R) {
        debug_assert!(self.is_valid_row_index(i));
        let mut row = self.array.row_mut(i);
        for a in row.iter_mut() {
            *a *= r;
        }
    }

    pub fn mul_col(&mut self, i: usize, r: &R) {
        debug_assert!(self.is_valid_col_index(i));
        let mut col = self.array.column_mut(i);
        for a in col.iter_mut() {
            *a *= r;
        }
    }

    pub fn add_row_to(&mut self, i: usize, j: usize, r: &R) { 
        debug_assert!(self.is_valid_row_index(i));
        debug_assert!(self.is_valid_row_index(j));

        let (s_i, s_j) = (s![i, ..], s![j, ..]);
        let (row_i, row_j) = self.array.multi_slice_mut((s_i, s_j));
        ndarray::Zip::from(row_i).and(row_j).for_each(|x, y| { 
            let x = x as &R;
            *y += r * x;
        });
    }

    pub fn add_col_to(&mut self, i: usize, j: usize, r: &R) { 
        debug_assert!(self.is_valid_col_index(i));
        debug_assert!(self.is_valid_col_index(j));

        let (s_i, s_j) = (s![.., i], s![.., j]);
        let (col_i, col_j) = self.array.multi_slice_mut((s_i, s_j));
        ndarray::Zip::from(col_i).and(col_j).for_each(|x, y| { 
            let x = x as &R;
            *y += r * x;
        });
    }

    // Multiply [a, b; c, d] from left. 
    pub fn left_elementary(&mut self, comps: [&R; 4], i: usize, j: usize) { 
        debug_assert!(self.is_valid_row_index(i));
        debug_assert!(self.is_valid_row_index(j));

        let [a, b, c, d] = comps;
        let (s_i, s_j) = (s![i, ..], s![j, ..]);
        let (row_i, row_j) = self.array.multi_slice_mut((s_i, s_j));
        ndarray::Zip::from(row_i).and(row_j).for_each(|x_mut, y_mut| { 
            let (x, y) = (x_mut as &R, y_mut as &R);
            (*x_mut, *y_mut) = (
                a * x + b * y,
                c * x + d * y
            )
        });
    }

    // Multiply [a, c; b, d] from left. 
    pub fn right_elementary(&mut self, comps: [&R; 4], i: usize, j: usize) { 
        debug_assert!(self.is_valid_col_index(i));
        debug_assert!(self.is_valid_col_index(j));

        let [a, b, c, d] = comps;
        let (s_i, s_j) = (s![.., i], s![.., j]);
        let (col_i, col_j) = self.array.multi_slice_mut((s_i, s_j));
        ndarray::Zip::from(col_i).and(col_j).for_each(|x_mut, y_mut| { 
            let (x, y) = (x_mut as &R, y_mut as &R);
            (*x_mut, *y_mut) = (
                x * a + y * b,
                x * c + y * d
            )
        });
    }

    // private methods // 

    fn is_valid_row_index(&self, i: usize) -> bool { 
        (0..self.rows()).contains(&i)
    }

    fn is_valid_col_index(&self, j: usize) -> bool { 
        (0..self.cols()).contains(&j)
    }
}

impl<R> From<Array2<R>> for Mat<R> 
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn from(array: Array2<R>) -> Self {
        Self { array }
    }
}

impl<R> From<&SpMat<R>> for Mat<R> 
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn from(a: &SpMat<R>) -> Self {
        Mat::from(a.cs_mat().to_dense())
    }
}

impl<R> Index<[usize; 2]> for Mat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    type Output = R;
    fn index(&self, index: [usize; 2]) -> &R {
        &self.array[index]
    }
}

impl<R> IndexMut<[usize; 2]> for Mat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        &mut self.array[index]
    }
}

impl<R> Neg for Mat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Mat::from(-self.array)
    }
}

impl<R> Neg for &Mat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    type Output = Mat<R>;
    fn neg(self) -> Self::Output {
        Mat::from(-&self.array)
    }
}

#[auto_ops]
impl<R> AddAssign<&Mat<R>> for Mat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn add_assign(&mut self, rhs: &Self) {
        assert_eq!(self.shape(), rhs.shape());
        self.array += &rhs.array;
    }
}

#[auto_ops]
impl<R> SubAssign<&Mat<R>> for Mat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn sub_assign(&mut self, rhs: &Self) {
        assert_eq!(self.shape(), rhs.shape());
        self.array -= &rhs.array
    }
}

#[auto_ops]
impl<'a, 'b, R> Mul<&'b Mat<R>> for &'a Mat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = Mat<R>;

    fn mul(self, rhs: &'b Mat<R>) -> Self::Output {
        assert_eq!(self.cols(), rhs.rows());
        let (l, m, n) = (self.rows(), self.cols(), rhs.cols());
        let array = Array2::from_shape_fn((l, n), |(i, k)| {
            (0..m).map(|j| {
                &self[[i, j]] * &rhs[[j, k]]
            }).sum()
        });
        Mat::from(array)
    }
}

macro_rules! impl_ops {
    ($trait:ident) => {
        impl<R> $trait<Mat<R>> for Mat<R>
        where R: Ring, for<'x> &'x R: RingOps<R> {}

        impl<R> $trait<Mat<R>> for &Mat<R>
        where R: Ring, for<'x> &'x R: RingOps<R> {}
    };
}

impl_ops!(AddMonOps);
impl_ops!(AddGrpOps);
impl_ops!(MonOps);
impl_ops!(RingOps);

impl<R> MatType for Mat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn shape(&self) -> (usize, usize) {
        (self.array.nrows(), self.array.ncols())
    }

    fn zero(shape: (usize, usize)) -> Self { 
        Self::from(Array2::zeros(shape))
    }

    fn is_zero(&self) -> bool {
        self.array.iter().all(|a| a.is_zero())
    }

    fn id(size: usize) -> Self { 
        let array = Array2::from_diag_elem(size, R::one());
        Self::from(array)
    }

    fn is_id(&self) -> bool { 
        self.is_square() && self.array.indexed_iter().all(|((i, j), a)| 
            i == j && a.is_one() || i != j && a.is_zero()
        )
    }
}

#[cfg(test)]
mod tests { 
    use ndarray::array;
    use super::*;

    #[test]
    fn init() { 
        let a = Mat::from(array![[1,2,3],[4,5,6]]);

        assert_eq!(a.rows(), 2);
        assert_eq!(a.cols(), 3);
        assert_eq!(a.array(), array![[1,2,3],[4,5,6]]);
    }

    #[test]
    fn eq() {
        let a = Mat::from(array![[1,2,3],[4,5,6]]);
        let b = Mat::from(array![[1,2,3],[4,5,6]]);
        let c = Mat::from(array![[1,2],[3,4],[5,6]]);

        assert_eq!(a, b);
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

        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let a = Mat::from(array);
        assert!(!a.is_zero());
    }

    #[test]
    fn id() {
        let a: Mat<i32> = Mat::id(3);
        assert!(a.is_id());

        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let a = Mat::from(array);
        assert!(!a.is_id());
    }

    
    #[test]
    fn swap_rows() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = Mat::from(array);
        a.swap_rows(0, 1);
        assert_eq!(a.array, array![[4,5,6],[1,2,3],[7,8,9]]);
    }

    #[test]
    fn swap_cols() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = Mat::from(array);
        a.swap_cols(0, 1);
        assert_eq!(a.array, array![[2,1,3],[5,4,6],[8,7,9]]);
    }

    #[test]
    fn mul_row() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = Mat::from(array);
        a.mul_row(1, &10);
        assert_eq!(a.array, array![[1,2,3],[40,50,60],[7,8,9]]);
    }

    #[test]
    fn mul_col() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = Mat::from(array);
        a.mul_col(1, &10);
        assert_eq!(a.array, array![[1,20,3],[4,50,6],[7,80,9]]);
    }

    #[test]
    fn add_row_to() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = Mat::from(array);
        a.add_row_to(0, 1, &10);
        assert_eq!(a.array, array![[1,2,3],[14,25,36],[7,8,9]]);
    }

    #[test]
    fn add_col_to() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = Mat::from(array);
        a.add_col_to(0, 1, &10);
        assert_eq!(a.array, array![[1,12,3],[4,45,6],[7,78,9]]);
    }

    #[test]
    fn add() { 
        let a = Mat::from(array![[1,2,3],[4,5,6]]);
        let b = Mat::from(array![[8,2,4],[0,2,1]]);
        let c = a + b;
        assert_eq!(c, Mat::from(array![[9,4,7],[4,7,7]]));
    }

    #[test]
    fn sub() { 
        let a = Mat::from(array![[1,2,3],[4,5,6]]);
        let b = Mat::from(array![[8,2,4],[0,2,1]]);
        let c = a - b;
        assert_eq!(c, Mat::from(array![[-7,0,-1],[4,3,5]]));
    }

    #[test]
    fn neg() { 
        let a = Mat::from(array![[1,2,3],[4,5,6]]);
        assert_eq!(-a, Mat::from(array![[-1,-2,-3],[-4,-5,-6]]));
    }

    #[test]
    fn mul() { 
        let a = Mat::from(array![[1,2,3],[4,5,6]]);
        let b = Mat::from(array![[1,2],[1,-1],[0,2]]);
        let c = a * b;
        assert_eq!(c, Mat::from(array![[3,6],[9,15]]));
    }

    #[test]
    fn to_sparse() { 
        let dns = Mat::from(array![[1,2,3],[4,5,6]]);
        let sps = dns.to_sparse();
        assert_eq!(sps, SpMat::from_vec((2, 3), vec![1,2,3,4,5,6]));
    }

    #[test]
    fn from_sparse() { 
        let sps = SpMat::from_vec((2, 3), vec![1,2,3,4,5,6]);
        let dns = Mat::from(&sps);
        assert_eq!(dns, Mat::from(array![[1,2,3],[4,5,6]]));
    }
}