use ndarray::{Array2, s};
use sprs::{CsMat, TriMat};
use crate::math::traits::Ring;

pub struct DnsMat<R>
where R: Ring {
    array: Array2<R>
}

impl<R> DnsMat<R> where R: Ring {
    pub fn nrows(&self) -> usize { 
        self.array.nrows()
    }

    pub fn ncols(&self) -> usize { 
        self.array.ncols()
    }

    pub fn eye(size: usize) -> Self { 
        Self::from(CsMat::eye(size))
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

    pub fn mul_row(&mut self, i: usize, r: R) {
        debug_assert!(self.is_valid_row_index(i));
        let mut row = self.array.row_mut(i);
        for a in row.iter_mut() {
            *a = a.clone() * r.clone();
        }
    }

    pub fn mul_col(&mut self, i: usize, r: R) {
        debug_assert!(self.is_valid_col_index(i));
        let mut col = self.array.column_mut(i);
        for a in col.iter_mut() {
            *a = a.clone() * r.clone();
        }
    }

    pub fn add_row_to(&mut self, i: usize, j: usize, r: R) { 
        debug_assert!(self.is_valid_row_index(i));
        debug_assert!(self.is_valid_row_index(j));

        let (s_i, s_j) = (s![i, ..], s![j, ..]);
        let (row_i, row_j) = self.array.multi_slice_mut((s_i, s_j));
        ndarray::Zip::from(row_i).and(row_j).for_each(|x, y| { 
            *y = y.clone() + r.clone() * x.clone();
        });
    }

    pub fn add_col_to(&mut self, i: usize, j: usize, r: R) { 
        debug_assert!(self.is_valid_col_index(i));
        debug_assert!(self.is_valid_col_index(j));

        let (s_i, s_j) = (s![.., i], s![.., j]);
        let (col_i, col_j) = self.array.multi_slice_mut((s_i, s_j));
        ndarray::Zip::from(col_i).and(col_j).for_each(|x, y| { 
            *y = y.clone() + r.clone() * x.clone();
        });
    }

    // private methods // 

    fn is_valid_row_index(&self, i: usize) -> bool { 
        (0..self.nrows()).contains(&i)
    }

    fn is_valid_col_index(&self, j: usize) -> bool { 
        (0..self.ncols()).contains(&j)
    }

    fn is_valid_index(&self, i: usize, j: usize) -> bool { 
        self.is_valid_row_index(i) && self.is_valid_col_index(j)
    }
}

impl<R> From<CsMat<R>> for DnsMat<R> where R: Ring {
    fn from(sp: CsMat<R>) -> Self {
        DnsMat{ array: sp.to_dense() }
    }
}

#[cfg(test)]
mod tests { 
    use ndarray::array;
    use super::DnsMat;

    #[test]
    fn swap_rows() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = DnsMat{ array };
        a.swap_rows(0, 1);
        assert_eq!(a.array, array![[4,5,6],[1,2,3],[7,8,9]]);
    }

    #[test]
    fn swap_cols() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = DnsMat{ array };
        a.swap_cols(0, 1);
        assert_eq!(a.array, array![[2,1,3],[5,4,6],[8,7,9]]);
    }

    #[test]
    fn mul_row() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = DnsMat{ array };
        a.mul_row(1, 10);
        assert_eq!(a.array, array![[1,2,3],[40,50,60],[7,8,9]]);
    }

    #[test]
    fn mul_col() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = DnsMat{ array };
        a.mul_col(1, 10);
        assert_eq!(a.array, array![[1,20,3],[4,50,6],[7,80,9]]);
    }

    #[test]
    fn add_row_to() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = DnsMat{ array };
        a.add_row_to(0, 1, 10);
        assert_eq!(a.array, array![[1,2,3],[14,25,36],[7,8,9]]);
    }

    #[test]
    fn add_col_to() { 
        let array = array![[1,2,3],[4,5,6],[7,8,9]];
        let mut a = DnsMat{ array };
        a.add_col_to(0, 1, 10);
        assert_eq!(a.array, array![[1,12,3],[4,45,6],[7,78,9]]);
    }
}