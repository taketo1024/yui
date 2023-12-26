use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, MulAssign, Range};
use std::iter::zip;
use std::fmt::{Display, Debug};
use delegate::delegate;
use nalgebra_sparse::na::{Scalar, ClosedAdd, ClosedSub, ClosedMul};
use nalgebra_sparse::{CscMatrix, CooMatrix};
use num_traits::{Zero, One, ToPrimitive};
use auto_impl_ops::auto_ops;
use sprs::PermView;
use yui::{Ring, RingOps};
use crate::dense::*;
use super::sp_vec::SpVec;
use super::triang::TriangularType;

#[derive(Clone, PartialEq, Eq)]
pub struct SpMat<R> { 
    inner: CscMatrix<R>
}

impl<R> MatTrait for SpMat<R> {
    fn shape(&self) -> (usize, usize) {
        (self.inner.nrows(), self.inner.ncols())
    }
}

impl<R> SpMat<R> { 
    pub(crate) fn inner(&self) -> &CscMatrix<R> { 
        &self.inner
    }

    pub(crate) fn into_inner(self) -> CscMatrix<R> { 
        self.inner
    }

    pub fn data(&self) -> (&[usize], &[usize], &[R]) { 
        self.inner.csc_data()
    }

    pub fn disassemble(self) -> (Vec<usize>, Vec<usize>, Vec<R>) { 
        self.inner.disassemble()
    }

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

    pub fn is_triang(&self, t: TriangularType) -> bool
    where R: Zero {
        if self.nrows() != self.ncols() { 
            return false
        }

        if t.is_upper() { 
            self.iter_nz().all(|(i, j, _)| i <= j )
        } else { 
            self.iter_nz().all(|(i, j, _)| i >= j )
        }
    }
    
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize, &R)> { 
        self.inner.triplet_iter()
    }

    pub fn iter_nz(&self) -> impl Iterator<Item = (usize, usize, &R)>
    where R: Zero { 
        self.iter().filter(|e| !e.2.is_zero())
    }

    pub fn into_dense(self) -> Mat<R>
    where R: Scalar + Zero + ClosedAdd { 
        self.into()
    }

    pub fn nnz(&self) -> usize { 
        self.inner.nnz()
    }

    pub fn density(&self) -> f64 { 
        let (m, n) = self.shape();
        if m == 0 || n == 0 { 
            return 0.0
        }

        let nnz = self.nnz().to_f64().unwrap();
        let total = (m * n).to_f64().unwrap();

        nnz / total
    }

    pub fn redundancy(&self) -> f64
    where R: Zero { 
        let nnz = self.nnz().to_f64().unwrap();
        let red = self.iter().filter(|(_, _, a)| a.is_zero()).count().to_f64().unwrap();
        red / nnz
    }

    pub fn mean_weight(&self) -> f64
    where R: Ring, for<'x> &'x R: RingOps<R> { 
        let nnz = self.nnz().to_f64().unwrap();
        let w = self.iter().map(|(_, _, a)| a.c_weight()).sum::<f64>(); 
        w / nnz
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
        let csc = CscMatrix::from(&coo);
        Self::from(csc)
    }

    pub fn from_col_vecs<I>(nrows: usize, vecs: I) -> Self 
    where I: IntoIterator<Item = SpVec<R>> { 
        let mut col_offsets = vec![0];
        let mut row_indices = vec![];
        let mut values = vec![];

        for v in vecs.into_iter() { 
            assert_eq!(nrows, v.dim());
            let (_, mut v_rows, mut v_values) = v.into_inner().disassemble();

            row_indices.append(&mut v_rows);
            values.append(&mut v_values);
            col_offsets.push(row_indices.len());
        }

        let ncols = col_offsets.len() - 1;
        let csc = CscMatrix::try_from_csc_data(nrows, ncols, col_offsets, row_indices, values).unwrap();
        Self::from(csc)
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

    pub fn col_vec(&self, j: usize) -> SpVec<R>
    where R: Scalar + Zero + ClosedAdd { 
        let col = self.inner.col(j);
        let iter = Iterator::zip(
            col.row_indices().iter().cloned(), 
            col.values().iter().cloned()
        );
        SpVec::from_entries(self.nrows(), iter)
    }

    pub fn transpose(&self) -> Self { 
        self.inner.transpose().into()
    }

    pub fn extract<F>(&self, shape: (usize, usize), f: F) -> SpMat<R>
    where F: Fn(usize, usize) -> Option<(usize, usize)> { 
        SpMat::from_entries(shape, self.iter().filter_map(|(i, j, a)|
            f(i, j).map(|(i, j)| (i, j, a.clone()))
        ))
    }

    pub fn permute(&self, p: PermView, q: PermView) -> SpMat<R> { 
        self.extract(self.shape(), |i, j| Some((p.at(i), q.at(j))))
    }

    pub fn permute_rows(&self, p: PermView) -> SpMat<R> { 
        let id = PermView::identity(self.ncols());
        self.permute(p, id)
    }
    
    pub fn permute_cols(&self, q: PermView) -> SpMat<R> { 
        let id = PermView::identity(self.nrows());
        self.permute(id, q)
    }

    pub fn submat(&self, rows: Range<usize>, cols: Range<usize>) -> SpMat<R> { 
        let (i0, i1) = (rows.start, rows.end);
        let (j0, j1) = (cols.start, cols.end);

        assert!(i0 <= i1 && i1 <= self.nrows());
        assert!(j0 <= j1 && j1 <= self.ncols());

        let shape = (i1 - i0, j1 - j0);
        self.extract(shape, |i, j|
            (rows.contains(&i) && cols.contains(&j)).then( ||
                (i - i0, j - j0)
            )
        )
    }

    pub fn submat_rows(&self, rows: Range<usize>) -> SpMat<R> { 
        let n = self.ncols();
        self.submat(rows, 0 .. n)
    }

    pub fn submat_cols(&self, cols: Range<usize>) -> SpMat<R> { 
        let m = self.nrows();
        self.submat(0 .. m, cols)
    }

    pub fn combine_blocks(blocks: [&SpMat<R>; 4]) -> SpMat<R> {
        let [a, b, c, d] = blocks;

        assert_eq!(a.nrows(), b.nrows());
        assert_eq!(c.nrows(), d.nrows());
        assert_eq!(a.ncols(), c.ncols());
        assert_eq!(b.ncols(), d.ncols());

        let (m, n) = (a.nrows() + c.nrows(), a.ncols() + b.ncols());
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
            &zero(0, self.ncols()), 
            &zero(0, b.ncols())
        ])
    }

    pub fn stack(&self, b: &Self) -> Self { 
        let zero = |m, n| SpMat::<R>::zero((m, n));
        Self::combine_blocks([
            self, 
            &zero(self.nrows(), 0), 
            b, 
            &zero(b.nrows(), 0)
        ])
    }

    pub fn extend_cols(&mut self, b: Self) { 
        assert_eq!(self.nrows(), b.nrows());

        if b.ncols() == 0 { 
            return
        }

        let shape = (self.nrows(), self.ncols() + b.ncols());
        let l = std::mem::replace(&mut self.inner, CscMatrix::zeros(0, 0));
        let r = b.inner;

        let (mut col_offsets, mut row_indices, mut values) = l.disassemble();
        let (c, mut r, mut v) = r.disassemble();
        
        let offset = col_offsets.pop().unwrap(); // pop last element.
        col_offsets.extend(c.into_iter().map(|i| offset + i));
        row_indices.append(&mut r);
        values.append(&mut v);

        self.inner = CscMatrix::try_from_csc_data(
            shape.0, shape.1, 
            col_offsets, 
            row_indices, 
            values
        ).unwrap();
    }

    // row_perm(p) * a == a.permute_rows(p)
    pub fn from_row_perm(p: PermView) -> Self
    where R: One {
        let n = p.dim();
        Self::from_entries((n, n), (0..n).map(|i|
            (p.at(i), i, R::one())
        ))
    }

    // a * col_perm(p) == a.permute_cols(p)
    pub fn from_col_perm(p: PermView) -> Self
    where R: One {
        let n = p.dim();
        Self::from_entries((n, n), (0..n).map(|i|
            (i, p.at(i), R::one())
        ))
    }
}

impl<R> From<CscMatrix<R>> for SpMat<R> {
    fn from(inner: CscMatrix<R>) -> Self {
        Self { inner }
    }
}

impl<R> From<Mat<R>> for SpMat<R>
where R: Scalar + Zero {
    fn from(value: Mat<R>) -> Self {
        let csc = CscMatrix::from(value.inner());
        Self::from(csc)
    }
}

impl<R> Default for SpMat<R> {
    fn default() -> Self {
        Self::zero((0, 0))
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

#[cfg(feature = "serde")]
impl<R> serde::Serialize for SpMat<R>
where R: Clone + serde::Serialize {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where S: serde::Serializer {
        self.inner.serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, R> serde::Deserialize<'de> for SpMat<R>
where R: Clone + serde::Deserialize<'de> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where D: serde::Deserializer<'de> {
        let inner = CscMatrix::deserialize(deserializer)?;
        let res = Self::from(inner);
        Ok(res)
    }
}

#[cfg(test)]
impl<R> SpMat<R>
where R: Scalar + Zero + One + ClosedAdd { 
    pub fn rand(shape: (usize, usize), density: f64) -> Self {
        use cartesian::cartesian;
        use rand::Rng;
    
        let (m, n) = shape;
        let range = cartesian!(0..m, 0..n);
        let mut rng = rand::thread_rng();
    
        Self::from_entries(shape, range.filter_map(|(i, j)|
            if rng.gen::<f64>() < density { 
                Some((i, j, R::one()))
            } else { 
                None
            }
        ))
    }
}

#[cfg(test)]
pub(super) mod tests { 
    use itertools::Itertools;
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

    #[test]
    fn to_dense() { 
        let a = SpMat::from_entries((2, 2), [
            (0, 0, 1),
            (0, 1, 2),
            (1, 0, 3),
            (1, 1, 4)
        ]);
        assert_eq!(a.into_dense(), Mat::from_data((2, 2), [1,2,3,4]));
    }

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
    fn submat() { 
        let a = SpMat::from_dense_data((5, 6), 0..30);
        let b = a.submat(1..3, 2..5);
        assert_eq!(b, SpMat::from_dense_data((2,3), vec![
             8,  9, 10,
            14, 15, 16
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

    #[test]
    fn extend_cols() {
        let mut a = SpMat::from_dense_data((4, 3), 0..12);
        let b = SpMat::from_dense_data((4, 2), 12..20);
        a.extend_cols(b);

        assert_eq!(a, SpMat::from_dense_data((4,5), vec![
            0,  1,  2, 12, 13,
            3,  4,  5, 14, 15,
            6,  7,  8, 16, 17,
            9, 10, 11, 18, 19,
        ]));
    }

    #[test]
    fn row_perm() {
        let a = SpMat::from_dense_data((3, 4), 0..12);
        let p = PermOwned::new(vec![2,0,1]);
        let q = SpMat::from_row_perm(p.view());
        assert!(q * &a == a.permute_rows(p.view()))
    }

    #[test]
    fn col_perm() {
        let a = SpMat::from_dense_data((3, 4), 0..12);
        let p = PermOwned::new(vec![2,0,1,3]);
        let q = SpMat::from_col_perm(p.view());
        assert!(&a * q == a.permute_cols(p.view()))
    }

    #[test]
    #[cfg(feature = "serde")]
    fn serialize() { 
        let a = SpMat::from_dense_data((3, 4), (0..12).map(|x| x % 5));
        let ser = serde_json::to_string(&a).unwrap();
        let des = serde_json::from_str(&ser).unwrap();
        assert_eq!(a, des);
    }
}