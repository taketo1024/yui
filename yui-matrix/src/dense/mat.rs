use std::ops::{Add, Neg, Sub, Mul, Index, IndexMut, AddAssign, SubAssign, MulAssign, Range};
use std::fmt::Debug;
use ndarray::{Array2, Axis, s, concatenate, array};
use derive_more::Display;
use auto_impl_ops::auto_ops;
use num_traits::{Zero, One};
use yui::{AddMon, AddMonOps, AddGrp, AddGrpOps, MonOps, Ring, RingOps};
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

#[derive(Clone, Debug, Display, Default, PartialEq, Eq)]
pub struct Mat<R> {
    array: Array2<R>
}

impl<R> MatType for Mat<R> {
    fn shape(&self) -> (usize, usize) {
        (self.array.nrows(), self.array.ncols())
    }
}

impl<R> Mat<R> {
    pub fn array(&self) -> &Array2<R> {
        &self.array
    }

    pub fn array_mut(&mut self) -> &mut Array2<R> {
        &mut self.array
    }

    pub fn array_into(self) -> Array2<R> {
        self.array
    }

    fn is_valid_row_index(&self, i: usize) -> bool { 
        i < self.rows()
    }

    fn is_valid_col_index(&self, j: usize) -> bool { 
        j < self.cols()
    }
}

impl<R> From<Array2<R>> for Mat<R> {
    fn from(array: Array2<R>) -> Self {
        Self { array }
    }
}

impl<R> From<SpMat<R>> for Mat<R>
where R: Clone + Zero {
    fn from(a: SpMat<R>) -> Self {
        Mat::from(a.cs_mat().to_dense())
    }
}

impl<R> Index<[usize; 2]> for Mat<R> {
    type Output = R;
    fn index(&self, index: [usize; 2]) -> &R {
        &self.array[index]
    }
}

impl<R> IndexMut<[usize; 2]> for Mat<R> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        &mut self.array[index]
    }
}

impl<R> Mat<R>
where R: Clone {
    pub fn clone_row(&self, i: usize) -> Vec<R> {
        assert!(self.is_valid_row_index(i));
        self.array.slice(s![i, ..]).iter().cloned().collect()
    }

    pub fn clone_col(&self, j: usize) -> Vec<R> {
        assert!(self.is_valid_col_index(j));
        self.array.slice(s![.., j]).iter().cloned().collect()
    }

    pub fn combine_blocks(blocks: [&Mat<R>; 4]) -> Mat<R> {
        let [a, b, c, d] = blocks;

        assert_eq!(a.rows(), b.rows());
        assert_eq!(c.rows(), d.rows());
        assert_eq!(a.cols(), c.cols());
        assert_eq!(b.cols(), d.cols());

        let ab = concatenate![Axis(1), a.array, b.array];
        let cd = concatenate![Axis(1), c.array, d.array];
        let abcd  = concatenate![Axis(0), ab, cd];

        Mat::from(abcd)
    }

    pub fn concat(&self, b: &Self) -> Self { 
        let a = self;

        assert_eq!(a.rows(), b.rows());

        let ab = concatenate![Axis(1), a.array, b.array];
        Mat::from(ab)
    }

    pub fn stack(&self, c: &Self) -> Self {
        let a = self;

        assert_eq!(a.cols(), c.cols());

        let ac = concatenate![Axis(0), a.array, c.array];
        Mat::from(ac)
    }

    pub fn submat(&self, rows: Range<usize>, cols: Range<usize>) -> Mat<R> { 
        let (i0, i1) = (rows.start, rows.end);
        let (j0, j1) = (cols.start, cols.end);

        assert!(i0 <= i1 && i1 <= self.rows());
        assert!(j0 <= j1 && j1 <= self.cols());

        let slice = self.array.slice(s![rows, cols]);
        Self::from(slice.to_owned())
    }

    pub fn submat_rows(&self, rows: Range<usize>) -> Mat<R> { 
        let n = self.cols();
        self.submat(rows, 0 .. n)
    }

    pub fn submat_cols(&self, cols: Range<usize>) -> Mat<R> { 
        let m = self.rows();
        self.submat(0 .. m, cols)
    }
}

impl<R> Mat<R>
where R: Clone + Zero {
    pub fn zero(shape: (usize, usize)) -> Self { 
        Self::from(Array2::zeros(shape))
    }

    pub fn is_zero(&self) -> bool {
        self.array.iter().all(|a| a.is_zero())
    }

    pub fn diag<I>(shape: (usize, usize), entries: I) -> Self
    where I: Iterator<Item = R> {
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

    pub fn iter(&self) -> impl Iterator<Item = (usize, usize, &R)> {
        self.array.indexed_iter().map(|((i, j), a)| (i, j, a))
    }

    pub fn to_sparse(self) -> SpMat<R> { 
        self.into()
    }
}

impl<R> Mat<R>
where R: Clone + Zero + One {
    pub fn id(size: usize) -> Self { 
        let array = Array2::from_diag_elem(size, R::one());
        Self::from(array)
    }

    pub fn is_id(&self) -> bool
    where R: PartialEq { 
        self.is_square() && self.array.indexed_iter().all(|((i, j), a)| 
            i == j && a.is_one() || i != j && a.is_zero()
        )
    }
}

impl<R> Neg for Mat<R>
where R: AddGrp, for<'x> &'x R: AddGrpOps<R> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Mat::from(-self.array)
    }
}

impl<R> Neg for &Mat<R>
where R: AddGrp, for<'x> &'x R: AddGrpOps<R> {
    type Output = Mat<R>;
    fn neg(self) -> Self::Output {
        Mat::from(-&self.array)
    }
}

#[auto_ops]
impl<R> AddAssign<&Mat<R>> for Mat<R>
where R: AddMon, for<'x> &'x R: AddMonOps<R> {
    fn add_assign(&mut self, rhs: &Self) {
        assert_eq!(self.shape(), rhs.shape());
        self.array += &rhs.array;
    }
}

#[auto_ops]
impl<R> SubAssign<&Mat<R>> for Mat<R>
where R: AddGrp, for<'x> &'x R: AddGrpOps<R> {
    fn sub_assign(&mut self, rhs: &Self) {
        assert_eq!(self.shape(), rhs.shape());
        self.array -= &rhs.array
    }
}

#[auto_ops]
impl<'a, 'b, R> Mul<&'b Mat<R>> for &'a Mat<R>
where R: AddGrp, for<'x> &'x R: AddGrpOps<R> + Mul<Output = R> {
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
    ($op_trait:ident, $r_trait:ident, $r_op_trait:ident) => {
        impl<R> $op_trait<Mat<R>> for Mat<R>
        where R: $r_trait, for<'x> &'x R: $r_op_trait<R> {}

        impl<R> $op_trait<Mat<R>> for &Mat<R>
        where R: $r_trait, for<'x> &'x R: $r_op_trait<R> {}
    };
}

impl_ops!(AddMonOps, AddMon, AddMonOps);
impl_ops!(AddGrpOps, AddGrp, AddGrpOps);
impl_ops!(MonOps, Ring, RingOps);
impl_ops!(RingOps, Ring, RingOps);

impl<R> Mat<R> where R: Clone { 
    pub fn insert_row(&mut self, row: Vec<R>, i: usize) {
        // TODO no cloning
        assert!(i <= self.rows());
        assert!(row.len() == self.cols());
        let above = self.array.slice(s![0..i, ..]);
        let below = self.array.slice(s![i..,  ..]);
        let row = Array2::from_shape_vec((1, row.len()), row).unwrap();
        self.array = concatenate![Axis(0), above, row, below];
    }

    pub fn insert_col(&mut self, col: Vec<R>, j: usize) {
        // TODO no cloning
        assert!(j <= self.cols());
        assert!(col.len() == self.rows());
        let left  = self.array.slice(s![.., 0..j]);
        let right = self.array.slice(s![.., j.. ]);
        let col = Array2::from_shape_vec((col.len(), 1), col).unwrap();
        self.array = concatenate![Axis(1), left, col, right];
    }
}

impl<R> Mat<R> where R: Clone + Zero { 
    pub fn insert_zero_row(&mut self, i: usize) {
        let row = vec![R::zero(); self.cols()];
        self.insert_row(row, i);
    }

    pub fn insert_zero_col(&mut self, j: usize) {
        let col = vec![R::zero(); self.rows()];
        self.insert_col(col, j);
    }
}

impl<R> Mat<R> { 
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

    pub fn mut_row<F>(&mut self, i: usize, f: F)
    where F: Fn(&mut R) {
        debug_assert!(self.is_valid_row_index(i));
        let mut row = self.array.row_mut(i);
        for a in row.iter_mut() {
            f(a)
        }
    }

    pub fn mut_col<F>(&mut self, j: usize, f: F)
    where F: Fn(&mut R) {
        debug_assert!(self.is_valid_col_index(j));
        let mut col = self.array.column_mut(j);
        for a in col.iter_mut() {
            f(a)
        }
    }

    pub fn del_row(&mut self, i: usize) {
        assert!(self.is_valid_row_index(i));
        let (m, n) = self.shape();
        let arr = std::mem::replace(&mut self.array, array![[]]);
        let vec = arr.into_iter().enumerate().filter(|(idx, _)| idx / n != i).map(|e| e.1).collect();
        self.array = Array2::from_shape_vec((m-1, n), vec).unwrap();
    }

    pub fn del_col(&mut self, j: usize) {
        assert!(self.is_valid_col_index(j));
        let (m, n) = self.shape();
        let arr = std::mem::replace(&mut self.array, array![[]]);
        let vec = arr.into_iter().enumerate().filter(|(idx, _)| idx % n != j).map(|e| e.1).collect();
        self.array = Array2::from_shape_vec((m, n-1), vec).unwrap();
    }
}

impl<R> Mat<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn mul_row(&mut self, i: usize, r: &R) {
        self.mut_row(i, |a| *a *= r);
    }

    pub fn mul_col(&mut self, j: usize, r: &R) {
        self.mut_col(j, |a| *a *= r);
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
        let dns = Mat::from(sps);
        assert_eq!(dns, Mat::from(array![[1,2,3],[4,5,6]]));
    }

    #[test]
    fn combine_blocks() { 
        let a = Mat::from(array![
            [1,2,3],
            [4,5,6]
        ]);
        let b = Mat::from(array![
            [7],
            [8]
        ]);
        let c = Mat::from(array![
            [9, 10, 11]
        ]);
        let d = Mat::from(array![
            [12]
        ]);
        let x = Mat::combine_blocks([&a,&b,&c,&d]);
        assert_eq!(x, Mat::from(array![
            [1, 2, 3, 7],
            [4, 5, 6, 8],
            [9,10,11,12]           
        ]))
    }

    #[test]
    fn concat() { 
        let a = Mat::from(array![
            [1,2,3],
            [4,5,6]
        ]);
        let b = Mat::from(array![
            [7],
            [8]
        ]);
        let x = a.concat(&b);
        assert_eq!(x, Mat::from(array![
            [1, 2, 3, 7],
            [4, 5, 6, 8]
        ]))
    }

    #[test]
    fn stack() { 
        let a = Mat::from(array![
            [1,2,3],
            [4,5,6]
        ]);
        let c = Mat::from(array![
            [9, 10, 11]
        ]);
        let x = a.stack(&c);
        assert_eq!(x, Mat::from(array![
            [1, 2, 3],
            [4, 5, 6],
            [9,10,11]           
        ]))
    }

    #[test]
    fn submat() { 
        let a = Mat::from(array![
            [1, 2, 3, 7],
            [4, 5, 6, 8],
            [9,10,11,12]           
        ]);
        let b = a.submat(1..3, 2..4);
        assert_eq!(b, Mat::from(array![
            [ 6, 8],
            [11,12]           
        ]));
    }

    #[test]
    fn insert_row() {
        let mut a = Mat::from(array![
            [1,2,3],
            [4,5,6],
            [7,8,9],
        ]);
        let v = vec![10,20,30];
        a.insert_row(v, 2);
        
        assert_eq!(a, Mat::from(array![
            [1,2,3],
            [4,5,6],
            [10,20,30],
            [7,8,9],
        ]))
    }

    #[test]
    fn insert_col() {
        let mut a = Mat::from(array![
            [1,2,3],
            [4,5,6],
            [7,8,9],
        ]);
        let v = vec![10,20,30];
        a.insert_col(v, 2);
        
        assert_eq!(a, Mat::from(array![
            [1,2,10,3],
            [4,5,20,6],
            [7,8,30,9],
        ]))
    }

    #[test]
    fn del_row() {
        let mut a = Mat::from(array![
            [1,2,3],
            [4,5,6],
            [7,8,9],
        ]);
        a.del_row(1);
        
        assert_eq!(a, Mat::from(array![
            [1,2,3],
            [7,8,9],
        ]))
    }

    #[test]
    fn del_col() {
        let mut a = Mat::from(array![
            [1,2,3],
            [4,5,6],
            [7,8,9],
        ]);
        a.del_col(1);
        
        assert_eq!(a, Mat::from(array![
            [1,3],
            [4,6],
            [7,9],
        ]))
    }
}