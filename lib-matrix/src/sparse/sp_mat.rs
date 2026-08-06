//! [`SpMat<R>`]: a sparse matrix in CSC form, over `nalgebra_sparse::CscMatrix`.
//! The carrier of every differential in this workspace.

use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, MulAssign, Range};
use std::fmt::{Display, Debug};
use delegate::delegate;
use itertools::{Itertools, repeat_n};
use nalgebra_sparse::na::{Scalar, ClosedAddAssign, ClosedSubAssign, ClosedMulAssign};
use nalgebra_sparse::{CscMatrix, CooMatrix};
use num_traits::{Zero, One, ToPrimitive};
use auto_impl_ops::auto_ops;
use yui_core::abst::{Ring, RingOps};
use crate::Perm;
use crate::dense::*;
use super::sp_vec::SpVec;
use super::triang::TriangularType;

/// Sparse matrix in compressed sparse column (CSC) format, backed by
/// `nalgebra_sparse::CscMatrix`.
///
/// The fundamental matrix type for all differential maps in the homology
/// pipeline. Generic over the element ring `R`; bounds are applied per-method.
///
/// Note: the underlying CSC may carry explicit zero entries (e.g. after an
/// in-place subtraction). [`iter`](Self::iter) walks all stored triplets;
/// [`iter_nz`](Self::iter_nz) filters them.
#[derive(Clone)]
pub struct SpMat<R> {
    inner: CscMatrix<R>
}

impl<R> SpMat<R> {
    pub fn try_from_csc_data(
        num_rows: usize,
        num_cols: usize,
        col_offsets: Vec<usize>,
        row_indices: Vec<usize>,
        values: Vec<R>,
    ) -> Option<Self> {
        let csc = CscMatrix::try_from_csc_data(num_rows, num_cols, col_offsets, row_indices, values);
        csc.ok().map(SpMat::from)
    }

    pub(crate) fn inner(&self) -> &CscMatrix<R> {
        &self.inner
    }

    pub(crate) fn into_inner(self) -> CscMatrix<R> {
        self.inner
    }

    pub fn csc_data(&self) -> (&[usize], &[usize], &[R]) {
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
        if self.n_rows() != self.n_cols() {
            return false
        }

        if t.is_upper() {
            self.iter_nz().all(|(i, j, _)| i <= j )
        } else {
            self.iter_nz().all(|(i, j, _)| i >= j )
        }
    }

    /// Iterates the stored `(row, col, value)` triplets — may include explicit zeros.
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize, &R)> {
        self.inner.triplet_iter()
    }

    /// [`iter`](Self::iter), filtered to non-zero values.
    pub fn iter_nz(&self) -> impl Iterator<Item = (usize, usize, &R)>
    where R: Zero {
        self.iter().filter(|e| !e.2.is_zero())
    }

    pub fn into_dense(self) -> Mat<R>
    where R: Scalar + Zero + ClosedAddAssign {
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

    pub fn block_diag<I>(blocks: I) -> SpMat<R>
    where I: IntoIterator<Item = SpMat<R>> {
        let mut shape = (0, 0);
        let mut col_offsets: Vec<usize> = vec![];
        let mut row_indices: Vec<usize> = vec![];
        let mut values: Vec<R> = vec![];

        for a in blocks {
            let a_shape = a.shape();
            let (a_cols, a_rows, mut a_vals) = a.disassemble();

            col_offsets.extend(a_cols.iter().map(|i| i + values.len()));
            col_offsets.pop(); // remove last offset

            row_indices.extend(a_rows.iter().map(|i| i + shape.0));
            values.append(&mut a_vals);

            shape.0 += a_shape.0;
            shape.1 += a_shape.1;
        }
        col_offsets.push(values.len());

        SpMat::try_from_csc_data(shape.0, shape.1, col_offsets, row_indices, values).unwrap()
    }

    pub fn map<F, S>(self, f: F) -> SpMat<S>
    where F: Fn(R) -> S {
        let (m, n) = self.shape();
        let (cols, rows, vals) = self.disassemble();
        let vals = vals.into_iter().map(f).collect_vec();
        SpMat::<S>::try_from_csc_data(m, n, cols, rows, vals).unwrap()
    }

    /// Returns the raw `(row_indices, values)` slices of column `j`.
    /// Borrow-only — no allocation, no value clones.
    pub(crate) fn col_data(&self, j: usize) -> (&[usize], &[R]) {
        let (col_offsets, row_indices, values) = self.inner.csc_data();
        let range = col_offsets[j]..col_offsets[j + 1];
        (&row_indices[range.clone()], &values[range])
    }
}

impl<R> SpMat<R>
where R: Scalar + Clone + Zero + ClosedAddAssign {
    /// Builds an `SpMat` of `shape` from `(row, col, value)` triplets. Zero
    /// values are skipped; duplicates at the same position are summed.
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

    // note: explicitly-stored zeros are dropped (CSC arithmetic keeps cancellation zeros).
    pub fn from_col_vecs<I>(nrows: usize, vecs: I) -> Self
    where I: IntoIterator<Item = SpVec<R>> {
        let mut col_offsets = vec![0];
        let mut row_indices = vec![];
        let mut values = vec![];

        for v in vecs.into_iter() {
            assert_eq!(nrows, v.dim());
            let (_, mut v_rows, mut v_values) = v.drop_zeros().into_inner().disassemble();

            row_indices.append(&mut v_rows);
            values.append(&mut v_values);
            col_offsets.push(row_indices.len());
        }

        let ncols = col_offsets.len() - 1;
        SpMat::try_from_csc_data(nrows, ncols, col_offsets, row_indices, values).unwrap()
    }

    pub fn from_row_major<I>(shape: (usize, usize), data: I) -> Self
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

    pub fn scalar(n: usize, a: &R) -> Self {
        Self::from_entries((n, n), (0..n).map(|i| (i, i, a.clone())))
    }

    pub fn diag<I>(shape: (usize, usize), entries: I) -> Self
    where I: IntoIterator<Item = R> {
        Self::from_entries(shape, entries.into_iter().enumerate().map(|(i, a)| (i, i, a)))
    }

    pub fn is_diag(&self) -> bool {
        self.iter_nz().all(|(i, j, _)| i == j)
    }

    pub fn col_vec(&self, j: usize) -> SpVec<R>
    where R: Scalar + Zero + ClosedAddAssign {
        let col = self.inner.col(j);
        let row_indices = col.row_indices().to_vec();
        let values = col.values().to_vec();
        SpVec::try_from_csc_data(self.n_rows(), row_indices, values).unwrap()
    }

    pub fn transpose(&self) -> Self {
        self.inner.transpose().into()
    }

    /// New `shape`-d matrix whose entry at `f(i, j)` (if `Some`) is `self[(i, j)]`.
    /// Entries where `f` returns `None` are dropped.
    pub fn extract<F>(&self, shape: (usize, usize), f: F) -> SpMat<R>
    where F: Fn(usize, usize) -> Option<(usize, usize)> {
        SpMat::from_entries(shape, self.iter().filter_map(|(i, j, a)|
            f(i, j).map(|(i, j)| (i, j, a.clone()))
        ))
    }

    pub fn permute(&self, p: &Perm, q: &Perm) -> SpMat<R> {
        self.extract(self.shape(), |i, j| Some((p.at(i), q.at(j))))
    }

    pub fn permute_rows(&self, p: &Perm) -> SpMat<R> {
        let id = Perm::id(self.n_cols());
        self.permute(p, &id)
    }

    pub fn permute_cols(&self, q: &Perm) -> SpMat<R> {
        let id = Perm::id(self.n_rows());
        self.permute(&id, q)
    }

    /// Applies the permutations `(p, q)` to `self` and partitions the result
    /// into four blocks at row/col `r`:
    ///
    /// ```text
    ///   paq = [[a0 | a1],   a0: r×r,     a1: r×(n-r)
    ///          [a2 | a3]]   a2: (m-r)×r, a3: (m-r)×(n-r)
    /// ```
    pub fn permute_and_split(&self, p: &Perm, q: &Perm, r: usize) -> [SpMat<R>; 4] {
        use std::cmp::Ordering::Less;

        let (m, n) = self.shape();
        assert!(r <= m && r <= n);

        let [mut a0, mut a1, mut a2, mut a3] = [vec![], vec![], vec![], vec![]];

        for (i, j, v) in self.iter() {
            let (pi, qj) = (p.at(i), q.at(j));
            let v = v.clone();
            match (pi.cmp(&r), qj.cmp(&r)) {
                (Less, Less) => a0.push((pi,     qj,     v)),
                (Less, _   ) => a1.push((pi,     qj - r, v)),
                (_   , Less) => a2.push((pi - r, qj,     v)),
                (_   , _   ) => a3.push((pi - r, qj - r, v)),
            }
        }
        [
            SpMat::from_entries((r,     r    ), a0),
            SpMat::from_entries((r,     n - r), a1),
            SpMat::from_entries((m - r, r    ), a2),
            SpMat::from_entries((m - r, n - r), a3),
        ]
    }

    pub fn submat(&self, rows: Range<usize>, cols: Range<usize>) -> SpMat<R> {
        let (i0, i1) = (rows.start, rows.end);
        let (j0, j1) = (cols.start, cols.end);

        assert!(i0 <= i1 && i1 <= self.n_rows());
        assert!(j0 <= j1 && j1 <= self.n_cols());

        let shape = (i1 - i0, j1 - j0);
        self.extract(shape, |i, j|
            (rows.contains(&i) && cols.contains(&j)).then( ||
                (i - i0, j - j0)
            )
        )
    }

    pub fn submat_rows(&self, rows: Range<usize>) -> SpMat<R> {
        let n = self.n_cols();
        self.submat(rows, 0 .. n)
    }

    pub fn submat_cols(&self, cols: Range<usize>) -> SpMat<R> {
        let m = self.n_rows();
        self.submat(0 .. m, cols)
    }

    pub fn block_split(self, point: (usize, usize)) -> [SpMat<R>; 4] {
        let (m, n) = self.shape();
        let (k, l) = point;
        assert!(k <= m);
        assert!(l <= n);

        let (offsets, rows, vals) = self.disassemble();

        let (mut a_rows, mut a_vals, mut a_offs) = (vec![], vec![], vec![0]);
        let (mut b_rows, mut b_vals, mut b_offs) = (vec![], vec![], vec![0]);
        let (mut c_rows, mut c_vals, mut c_offs) = (vec![], vec![], vec![0]);
        let (mut d_rows, mut d_vals, mut d_offs) = (vec![], vec![], vec![0]);

        let mut vals_iter = vals.into_iter();

        for j in 0..n {
            let range = offsets[j]..offsets[j + 1];
            let col_rows = &rows[range];
            let split = col_rows.partition_point(|&i| i < k);
            let (top_rows_src, bot_rows_src) = col_rows.split_at(split);

            let (top_rows, top_vals, top_offs, bot_rows, bot_vals, bot_offs) = if j < l {
                (&mut a_rows, &mut a_vals, &mut a_offs, &mut c_rows, &mut c_vals, &mut c_offs)
            } else {
                (&mut b_rows, &mut b_vals, &mut b_offs, &mut d_rows, &mut d_vals, &mut d_offs)
            };

            top_rows.extend_from_slice(top_rows_src);
            top_vals.extend(vals_iter.by_ref().take(top_rows_src.len()));
            top_offs.push(top_rows.len());

            bot_rows.extend(bot_rows_src.iter().map(|&i| i - k));
            bot_vals.extend(vals_iter.by_ref().take(bot_rows_src.len()));
            bot_offs.push(bot_rows.len());
        }

        [
            SpMat::try_from_csc_data(k,     l,     a_offs, a_rows, a_vals).unwrap(),
            SpMat::try_from_csc_data(k,     n - l, b_offs, b_rows, b_vals).unwrap(),
            SpMat::try_from_csc_data(m - k, l,     c_offs, c_rows, c_vals).unwrap(),
            SpMat::try_from_csc_data(m - k, n - l, d_offs, d_rows, d_vals).unwrap(),
        ]
    }

    pub fn h_split(self, k: usize) -> [SpMat<R>; 2] {
        let (m, n) = self.shape();
        assert!(k <= n);

        let [a, b, ..] = self.block_split((m, k));
        [a, b]
    }

    pub fn v_split(self, k: usize) -> [SpMat<R>; 2] {
        let (m, n) = self.shape();
        assert!(k <= m);

        let [a, _, b, _] = self.block_split((k, n));
        [a, b]
    }

    pub fn block_combine(blocks: [SpMat<R>; 4]) -> SpMat<R> {
        let [a, b, c, d] = blocks;

        assert_eq!(a.n_rows(), b.n_rows());
        assert_eq!(c.n_rows(), d.n_rows());
        assert_eq!(a.n_cols(), c.n_cols());
        assert_eq!(b.n_cols(), d.n_cols());

        let (m0, m1) = (a.n_rows(), c.n_rows());
        let m = m0 + m1;
        let (n0, n1) = (a.n_cols(), b.n_cols());
        let n = n0 + n1;
        let nnz = a.nnz() + b.nnz() + c.nnz() + d.nnz();

        let mut a = ColSource::from(a);
        let mut b = ColSource::from(b);
        let mut c = ColSource::from(c);
        let mut d = ColSource::from(d);

        let mut col_offsets = Vec::with_capacity(n + 1);
        let mut row_indices = Vec::with_capacity(nnz);
        let mut values = Vec::with_capacity(nnz);
        col_offsets.push(0);

        let mut push_col = |top: &mut ColSource<R>, bot: &mut ColSource<R>, j: usize| {
            let (top_rows, top_vals) = top.take_col(j);
            row_indices.extend_from_slice(top_rows);
            values.extend(top_vals);

            let (bot_rows, bot_vals) = bot.take_col(j);
            row_indices.extend(bot_rows.iter().map(|i| i + m0));
            values.extend(bot_vals);

            col_offsets.push(row_indices.len());
        };

        for j in 0..n0 { push_col(&mut a, &mut c, j); }
        for j in 0..n1 { push_col(&mut b, &mut d, j); }

        SpMat::try_from_csc_data(m, n, col_offsets, row_indices, values).unwrap()
    }

    pub fn h_stack(left: Self, right: Self) -> Self {
        assert_eq!(left.n_rows(), right.n_rows());
        let (l_cols, r_cols) = (left.n_cols(), right.n_cols());
        Self::block_combine([
            left,
            right,
            SpMat::zero((0, l_cols)),
            SpMat::zero((0, r_cols)),
        ])
    }

    pub fn v_stack(top: Self, bot: Self) -> Self {
        assert_eq!(top.n_cols(), bot.n_cols());
        let (t_rows, b_rows) = (top.n_rows(), bot.n_rows());
        Self::block_combine([
            top,
            SpMat::zero((t_rows, 0)),
            bot,
            SpMat::zero((b_rows, 0)),
        ])
    }

    pub fn extend_by_zero(&mut self, add_rows: usize, add_cols: usize) {
        let (m, n) = self.shape();
        let l = std::mem::take(&mut self.inner);
        let (mut col_offsets, row_indices, values) = l.disassemble();
        let last = *col_offsets.last().unwrap();
        col_offsets.extend(repeat_n(last, add_cols));

        self.inner = CscMatrix::try_from_csc_data(
            m + add_rows, n + add_cols,
            col_offsets, row_indices, values
        ).unwrap();
    }

    // row_perm(p) * a == a.permute_rows(p)
    pub fn row_perm_mat(p: &Perm) -> Self
    where R: One {
        let n = p.len();
        Self::from_entries((n, n), (0..n).map(|i|
            (p.at(i), i, R::one())
        ))
    }

    // a * col_perm(p) == a.permute_cols(p)
    pub fn col_perm_mat(p: &Perm) -> Self
    where R: One {
        let n = p.len();
        Self::from_entries((n, n), (0..n).map(|i|
            (i, p.at(i), R::one())
        ))
    }
}

// A column-major view of a disassembled matrix that yields one column at a
// time, moving values out without cloning. Used by `block_combine`.
struct ColSource<R> {
    offsets: Vec<usize>,
    rows: Vec<usize>,
    vals: std::vec::IntoIter<R>,
    pos: usize,
}

impl<R> ColSource<R> {
    fn from(m: SpMat<R>) -> Self {
        let (offsets, rows, vals) = m.disassemble();
        Self { offsets, rows, vals: vals.into_iter(), pos: 0 }
    }

    // Returns `(row_indices, values)` for column `j`. Must be called with
    // monotonically increasing `j` since values are moved out lazily.
    fn take_col(&mut self, j: usize) -> (&[usize], impl Iterator<Item = R> + '_) {
        debug_assert_eq!(self.pos, self.offsets[j]);
        let range = self.offsets[j]..self.offsets[j + 1];
        let count = range.len();
        self.pos += count;
        (&self.rows[range], self.vals.by_ref().take(count))
    }
}

impl<R> MatTrait for SpMat<R> {
    fn shape(&self) -> (usize, usize) {
        (self.inner.nrows(), self.inner.ncols())
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

impl<R: PartialEq + Zero> PartialEq for SpMat<R> {
    fn eq(&self, other: &Self) -> bool {
        self.shape() == other.shape() && self.iter_nz().eq(other.iter_nz())
    }
}

impl<R: Eq + Zero> Eq for SpMat<R> {}

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
        impl<R> $trait<&SpMat<R>> for &SpMat<R>
        where R: Scalar + ClosedAddAssign + ClosedSubAssign + ClosedMulAssign + Zero + One + Neg<Output = R> {
            type Output = SpMat<R>;
            fn $method(self, rhs: &SpMat<R>) -> Self::Output {
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
where R: Scalar + Zero + One + ClosedAddAssign {
    pub fn rand(shape: (usize, usize), density: f64) -> Self {
        use itertools::iproduct;
        use rand::RngExt;

        let (m, n) = shape;
        let range = iproduct!(0..m, 0..n);
        let mut rng = rand::rng();

        Self::from_entries(shape, range.filter_map(|(i, j)|
            if rng.random::<f64>() < density {
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
    use yui_core::num::Ratio;
    use crate::Perm;

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
            (0, 0, vals[0]),
            (0, 1, vals[2]),
            (1, 0, vals[1]),
            (1, 1, vals[3])
        ]);
        assert_eq!(a.disassemble(), (vec![0, 2, 4], vec![0, 1, 0, 1], vals));
    }

    #[test]
    fn from_row_major() {
        let a = SpMat::from_row_major((2, 2), [1,2,3,4]);
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
        assert_eq!(a.into_dense(), Mat::from_row_major((2, 2), [1,2,3,4]));
    }

    #[test]
    fn permute() {
        let p = Perm::new(vec![1,2,3,0]);
        let q = Perm::new(vec![3,0,2,1]);
        let a = SpMat::from_row_major((4,4), 0..16);
        let b = a.permute(&p, &q);
        assert_eq!(b, SpMat::from_row_major((4,4), vec![
            13, 15, 14, 12,
             1,  3,  2,  0,
             5,  7,  6,  4,
             9, 11, 10,  8,
        ]));
    }

    #[test]
    fn permute_and_split_identity() {
        // a = [[1,2],[3,4]], r=1, identity perms → paq = a, partition at row/col 1:
        // a0=[[1]], a1=[[2]], a2=[[3]], a3=[[4]]
        let a = SpMat::from_row_major((2, 2), [1, 2, 3, 4]);
        let id = Perm::id(2);
        let [a0, a1, a2, a3] = a.permute_and_split(&id, &id, 1);
        assert_eq!(a0, SpMat::from_row_major((1, 1), [1]));
        assert_eq!(a1, SpMat::from_row_major((1, 1), [2]));
        assert_eq!(a2, SpMat::from_row_major((1, 1), [3]));
        assert_eq!(a3, SpMat::from_row_major((1, 1), [4]));
    }

    #[test]
    fn permute_and_split_with_perm() {
        // a (3×3) with row perm p = [2,0,1] (sends row 0→2, 1→0, 2→1)
        // and col perm q = [1,2,0] (sends col 0→1, 1→2, 2→0).
        // PAQ⁻¹[p(i), q(j)] = A[i, j]. Splitting at r=2 yields the 2×2 top-left,
        // 2×1 top-right, 1×2 bottom-left, 1×1 bottom-right blocks.
        let a = SpMat::from_row_major((3, 3), [
            1, 2, 3,
            4, 5, 6,
            7, 8, 9,
        ]);
        let p = Perm::new(vec![2, 0, 1]);
        let q = Perm::new(vec![1, 2, 0]);
        let [a0, a1, a2, a3] = a.permute_and_split(&p, &q, 2);
        // PAQ⁻¹ rows in image order = [row 1, row 2, row 0], cols = [col 2, col 0, col 1]:
        //   [[6, 4, 5],
        //    [9, 7, 8],
        //    [3, 1, 2]]
        assert_eq!(a0, SpMat::from_row_major((2, 2), [6, 4, 9, 7]));
        assert_eq!(a1, SpMat::from_row_major((2, 1), [5, 8]));
        assert_eq!(a2, SpMat::from_row_major((1, 2), [3, 1]));
        assert_eq!(a3, SpMat::from_row_major((1, 1), [2]));
    }

    #[test]
    fn submat() {
        let a = SpMat::from_row_major((5, 6), 0..30);
        let b = a.submat(1..3, 2..5);
        assert_eq!(b, SpMat::from_row_major((2,3), vec![
             8,  9, 10,
            14, 15, 16
        ]));
    }

    #[test]
    fn transpose() {
        let a = SpMat::from_row_major((3,4), 0..12);
        let b = a.transpose();

        assert_eq!(b, SpMat::from_row_major((4,3), vec![
            0, 4, 8,
            1, 5, 9,
            2, 6, 10,
            3, 7, 11,
        ]));
    }

    #[test]
    fn h_stack() {
        let a = SpMat::from_row_major((4, 3), 0..12);
        let b = SpMat::from_row_major((4, 2), 12..20);
        let c = SpMat::h_stack(a, b);

        assert_eq!(c, SpMat::from_row_major((4,5), vec![
            0,  1,  2, 12, 13,
            3,  4,  5, 14, 15,
            6,  7,  8, 16, 17,
            9, 10, 11, 18, 19,
        ]));
    }

    #[test]
    fn v_stack() {
        let a = SpMat::from_row_major((2, 3), 0..6);
        let b = SpMat::from_row_major((3, 3), 6..15);
        let c = SpMat::v_stack(a, b);

        assert_eq!(c, SpMat::from_row_major((5, 3), vec![
            0,  1,  2,
            3,  4,  5,
            6,  7,  8,
            9, 10, 11,
           12, 13, 14,
        ]));
    }

    #[test]
    fn extend_by_zero() {
        // [[1,2],[3,4]] extended by 1 row and 2 cols → [[1,2,0,0],[3,4,0,0],[0,0,0,0]]
        let mut a = SpMat::from_row_major((2, 2), [1,2,3,4]);
        a.extend_by_zero(1, 2);
        assert_eq!(a.shape(), (3, 4));
        assert_eq!(a, SpMat::from_row_major((3, 4), [1,2,0,0, 3,4,0,0, 0,0,0,0]));
    }

    #[test]
    fn row_perm() {
        let a = SpMat::from_row_major((3, 4), 0..12);
        let p = Perm::new(vec![2,0,1]);
        let q = SpMat::row_perm_mat(&p);
        assert!(q * &a == a.permute_rows(&p))
    }

    #[test]
    fn col_perm() {
        let a = SpMat::from_row_major((3, 4), 0..12);
        let p = Perm::new(vec![2,0,1,3]);
        let q = SpMat::col_perm_mat(&p);
        assert!(&a * q == a.permute_cols(&p))
    }

    #[test]
    fn diag() {
        let d = SpMat::diag((3, 4), [1, 2, 3]);
        assert_eq!(d, SpMat::from_row_major((3, 4), [
            1, 0, 0, 0,
            0, 2, 0, 0,
            0, 0, 3, 0,
        ]));
        assert!(d.is_diag());
    }

    #[test]
    fn is_diag_false() {
        let a = SpMat::from_row_major((2, 2), [1, 2, 0, 3]);
        assert!(!a.is_diag());
    }

    #[test]
    fn block_diag() {
        let a = SpMat::from_row_major((2, 2), 1..=4);
        let b = SpMat::from_row_major((1, 3), 5..=7);
        let c = SpMat::from_row_major((2, 1), 8..=9);
        let d = SpMat::block_diag([a, b, c]);
        assert_eq!(d, SpMat::from_row_major((5, 6), [
            1,2,0,0,0,0,
            3,4,0,0,0,0,
            0,0,5,6,7,0,
            0,0,0,0,0,8,
            0,0,0,0,0,9
        ]))
    }

    #[test]
    #[cfg(feature = "serde")]
    fn serialize() {
        let a = SpMat::from_row_major((3, 4), (0..12).map(|x| x % 5));
        let ser = serde_json::to_string(&a).unwrap();
        let des = serde_json::from_str(&ser).unwrap();
        assert_eq!(a, des);
    }
}