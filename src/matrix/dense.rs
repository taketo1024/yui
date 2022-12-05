use std::ops::{Add, Neg, Sub, Mul, Index, IndexMut, Deref};
use std::cmp::min;
use ndarray::{Array2, s};
use sprs::{CsMat, TriMat, CsMatBase, SpIndex};
use crate::math::traits::{Ring, RingOps};

#[derive(Clone, Debug, PartialEq)]
pub struct DnsMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    array: Array2<R>
}

impl<R> From<Array2<R>> for DnsMat<R> 
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn from(array: Array2<R>) -> Self {
        Self { array }
    }
}

impl<R, I, Iptr, IptrStorage, IndStorage, DataStorage> 
    From<&CsMatBase<R, I, IptrStorage, IndStorage, DataStorage, Iptr>> 
    for DnsMat<R>
where 
    R: Ring, for<'a> &'a R: RingOps<R> ,
    I: SpIndex,
    Iptr: SpIndex,
    IptrStorage: Deref<Target = [Iptr]>,
    IndStorage: Deref<Target = [I]>,
    DataStorage: Deref<Target = [R]>,
{
    fn from(sp: &CsMatBase<R, I, IptrStorage, IndStorage, DataStorage, Iptr>) -> Self {
        DnsMat{ array: sp.to_dense() }
    }
}


impl<R> DnsMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    pub fn shape(&self) -> (usize, usize) { 
        (self.nrows(), self.ncols())
    }

    pub fn nrows(&self) -> usize { 
        self.array.nrows()
    }

    pub fn ncols(&self) -> usize { 
        self.array.ncols()
    }

    pub fn array(&self) -> &Array2<R> {
        &self.array
    }

    pub fn is_square(&self) -> bool { 
        self.nrows() == self.ncols()
    }

    pub fn eye(size: usize) -> Self { 
        let array = Array2::from_diag_elem(size, R::one());
        Self::from(array)
    }

    pub fn is_eye(&self) -> bool { 
        self.is_square() && self.array.indexed_iter().all(|((i, j), a)| 
            i == j && a.is_one() || i != j && a.is_zero()
        )
    }

    pub fn zero(shape: (usize, usize)) -> Self { 
        Self::from(Array2::zeros(shape))
    }

    pub fn is_zero(&self) -> bool {
        self.array.iter().all(|a| a.is_zero())
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
        for a_mut in row.iter_mut() {
            let a = a_mut as &R;
            *a_mut = a * r;
        }
    }

    pub fn mul_col(&mut self, i: usize, r: &R) {
        debug_assert!(self.is_valid_col_index(i));
        let mut col = self.array.column_mut(i);
        for a_mut in col.iter_mut() {
            let a = a_mut as &R;
            *a_mut = a * r;
        }
    }

    pub fn add_row_to(&mut self, i: usize, j: usize, r: &R) { 
        debug_assert!(self.is_valid_row_index(i));
        debug_assert!(self.is_valid_row_index(j));

        let (s_i, s_j) = (s![i, ..], s![j, ..]);
        let (row_i, row_j) = self.array.multi_slice_mut((s_i, s_j));
        ndarray::Zip::from(row_i).and(row_j).for_each(|x_mut, y_mut| { 
            let (x, y) = (x_mut as &R, y_mut as &R);
            *y_mut = y + &(r * x);
        });
    }

    pub fn add_col_to(&mut self, i: usize, j: usize, r: &R) { 
        debug_assert!(self.is_valid_col_index(i));
        debug_assert!(self.is_valid_col_index(j));

        let (s_i, s_j) = (s![.., i], s![.., j]);
        let (col_i, col_j) = self.array.multi_slice_mut((s_i, s_j));
        ndarray::Zip::from(col_i).and(col_j).for_each(|x_mut, y_mut| { 
            let (x, y) = (x_mut as &R, y_mut as &R);
            *y_mut = y + &(r * x);
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

    pub fn to_sparse(&self) -> CsMat<R> {
        let (m, n) = (self.array.nrows(), self.array.ncols());
        let mut sp = TriMat::new((m, n));

        for (k, a) in self.array.iter().enumerate() {
            if a.is_zero() { continue }
            let (i, j) = (k / n, k % n);
            sp.add_triplet(i, j, a.clone());
        }

        sp.to_csc()
    }

    // private methods // 

    fn is_valid_row_index(&self, i: usize) -> bool { 
        (0..self.nrows()).contains(&i)
    }

    fn is_valid_col_index(&self, j: usize) -> bool { 
        (0..self.ncols()).contains(&j)
    }
}

impl<R> Add<Self> for DnsMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.shape(), rhs.shape());
        DnsMat::from(self.array + rhs.array)
    }
}

impl<R> Neg for DnsMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        DnsMat::from(-self.array)
    }
}

impl<R> Sub<Self> for DnsMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        assert_eq!(self.shape(), rhs.shape());
        DnsMat::from(self.array - rhs.array)
    }
}

impl<R> Mul<Self> for DnsMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        assert_eq!(self.ncols(), rhs.nrows());
        let (l, m, n) = (self.nrows(), self.ncols(), rhs.ncols());
        let array = Array2::from_shape_fn((l, n), |(i, k)| {
            (0..m).map(|j| {
                &self[[i, j]] * &rhs[[j, k]]
            }).sum()
        });
        DnsMat::from(array)
    }
}

impl<R> Index<[usize; 2]> for DnsMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    type Output = R;
    fn index(&self, index: [usize; 2]) -> &R {
        &self.array[index]
    }
}

impl<R> IndexMut<[usize; 2]> for DnsMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        &mut self.array[index]
    }
}

#[cfg(test)]
mod tests { 
    use ndarray::array;
    use super::*;

    #[test]
    fn init() { 
        let a = DnsMat::from(array![[1,2,3],[4,5,6]]);

        assert_eq!(a.nrows(), 2);
        assert_eq!(a.ncols(), 3);
        assert_eq!(a.array(), array![[1,2,3],[4,5,6]]);
    }

    #[test]
    fn eq() {
        let a = DnsMat::from(array![[1,2,3],[4,5,6]]);
        let b = DnsMat::from(array![[1,2,3],[4,5,6]]);
        let c = DnsMat::from(array![[1,2],[3,4],[5,6]]);

        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn square() {
        let a: DnsMat<i32> = DnsMat::zero((3, 3));
        assert!(a.is_square());

        let a: DnsMat<i32> = DnsMat::zero((3, 2));
        assert!(!a.is_square());
    }

    #[test]
    fn zero() {
        let a: DnsMat<i32> = DnsMat::zero((3, 2));
        assert!(a.is_zero());

        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let a = DnsMat::from(array);
        assert!(!a.is_zero());
    }

    #[test]
    fn eye() {
        let a: DnsMat<i32> = DnsMat::eye(3);
        assert!(a.is_eye());

        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let a = DnsMat::from(array);
        assert!(!a.is_eye());
    }

    
    #[test]
    fn swap_rows() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = DnsMat::from(array);
        a.swap_rows(0, 1);
        assert_eq!(a.array, array![[4,5,6],[1,2,3],[7,8,9]]);
    }

    #[test]
    fn swap_cols() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = DnsMat::from(array);
        a.swap_cols(0, 1);
        assert_eq!(a.array, array![[2,1,3],[5,4,6],[8,7,9]]);
    }

    #[test]
    fn mul_row() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = DnsMat::from(array);
        a.mul_row(1, &10);
        assert_eq!(a.array, array![[1,2,3],[40,50,60],[7,8,9]]);
    }

    #[test]
    fn mul_col() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = DnsMat::from(array);
        a.mul_col(1, &10);
        assert_eq!(a.array, array![[1,20,3],[4,50,6],[7,80,9]]);
    }

    #[test]
    fn add_row_to() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = DnsMat::from(array);
        a.add_row_to(0, 1, &10);
        assert_eq!(a.array, array![[1,2,3],[14,25,36],[7,8,9]]);
    }

    #[test]
    fn add_col_to() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = DnsMat::from(array);
        a.add_col_to(0, 1, &10);
        assert_eq!(a.array, array![[1,12,3],[4,45,6],[7,78,9]]);
    }

    #[test]
    fn add() { 
        let a = DnsMat::from(array![[1,2,3],[4,5,6]]);
        let b = DnsMat::from(array![[8,2,4],[0,2,1]]);
        let c = a + b;
        assert_eq!(c, DnsMat::from(array![[9,4,7],[4,7,7]]));
    }

    #[test]
    fn sub() { 
        let a = DnsMat::from(array![[1,2,3],[4,5,6]]);
        let b = DnsMat::from(array![[8,2,4],[0,2,1]]);
        let c = a - b;
        assert_eq!(c, DnsMat::from(array![[-7,0,-1],[4,3,5]]));
    }

    #[test]
    fn neg() { 
        let a = DnsMat::from(array![[1,2,3],[4,5,6]]);
        assert_eq!(-a, DnsMat::from(array![[-1,-2,-3],[-4,-5,-6]]));
    }

    #[test]
    fn mul() { 
        let a = DnsMat::from(array![[1,2,3],[4,5,6]]);
        let b = DnsMat::from(array![[1,2],[1,-1],[0,2]]);
        let c = a * b;
        assert_eq!(c, DnsMat::from(array![[3,6],[9,15]]));
    }

    #[test]
    fn to_sparse() { 
        let dns = DnsMat::from(array![[1,2,3],[4,5,6]]);
        let sps = dns.to_sparse();
        assert_eq!(sps, CsMat::new((2, 3), vec![0,3,6], vec![0,1,2,0,1,2], vec![1,2,3,4,5,6]).into_csc());
    }

    #[test]
    fn from_sparse() { 
        let sps = CsMat::new((2, 3), vec![0,3,6], vec![0,1,2,0,1,2], vec![1,2,3,4,5,6]);
        let dns = DnsMat::from(&sps);
        assert_eq!(dns, DnsMat::from(array![[1,2,3],[4,5,6]]));
    }
}