use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, Range};
use std::fmt::{Display, Debug};
use nalgebra_sparse::CscMatrix;
use nalgebra_sparse::na::{Scalar, ClosedAddAssign, ClosedSubAssign, ClosedMulAssign};
use num_traits::{Zero, One};
use auto_impl_ops::auto_ops;
use yui_core::{Ring, RingOps, AddGrpOps,  AddGrp};
use super::sp_mat::SpMat;
use crate::Perm;

#[derive(Clone, Debug)]
pub struct SpVec<R> { 
    inner: CscMatrix<R> // ncols == 1
}

impl<R> SpVec<R> { 
    fn new(inner: CscMatrix<R>) -> Self { 
        assert_eq!(inner.ncols(), 1);
        Self { inner }
    }

    #[allow(unused)]
    pub(crate) fn inner(&self) -> &CscMatrix<R> { 
        &self.inner
    }

    pub(crate) fn into_inner(self) -> CscMatrix<R> { 
        self.inner
    }

    pub fn data(&self) -> (&[usize], &[R]) { 
        let (_, indices, values) = self.inner.csc_data();
        (indices, values)
    }

    pub fn zero(dim: usize) -> Self {
        let inner = CscMatrix::zeros(dim, 1);
        Self::new(inner)
    }

    pub fn is_zero(&self) -> bool
    where R: Zero {
        self.inner.values().iter().all(|a| a.is_zero())
    }

    pub fn unit(n: usize, i: usize) -> Self
    where R: One {
        let inner = CscMatrix::try_from_csc_data(
            n, 1, 
            vec![0, 1], 
            vec![i], 
            vec![R::one()]
        ).unwrap();

        Self::new(inner)
    }

    pub fn dim(&self) -> usize { 
        self.inner.nrows()
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, &R)> { 
        self.inner.triplet_iter().map(|(i, _, a)| (i, a))
    }

    pub fn iter_nz(&self) -> impl Iterator<Item = (usize, &R)>
    where R: Zero { 
        self.iter().filter(|(_, a)| !a.is_zero())
    }

    pub fn into_dense(self) -> Vec<R>
    where R: Clone + Zero {
        self.into()
    }

    pub fn into_mat(self) -> SpMat<R> { 
        self.into()
    }
}

impl<R> From<Vec<R>> for SpVec<R>
where R: Scalar + Zero + ClosedAddAssign {
    fn from(vec: Vec<R>) -> Self {
        Self::from_entries(vec.len(), vec.into_iter().enumerate())
    }
}

impl<R> From<SpVec<R>> for Vec<R>
where R: Clone + Zero {
    fn from(value: SpVec<R>) -> Self {
        let mut res = vec![R::zero(); value.dim()];
        let (_, rows, vals) = value.inner.disassemble();
        for (i, a) in rows.into_iter().zip(vals) {
            res[i] = a;
        }
        res
    }
}

// SpVec(n) as SpMat(n, 1)
impl<R> From<SpVec<R>> for SpMat<R> { 
    fn from(vec: SpVec<R>) -> Self {
        SpMat::from(vec.into_inner())
    }
}

impl<R> SpMat<R> {
    fn into_spvec(self) -> SpVec<R> { 
        assert_eq!(self.inner().ncols(), 1);
        SpVec::new(self.into_inner())
    }
}

impl<R> SpVec<R> 
where R: Scalar + Zero + ClosedAddAssign { 
    pub fn try_from_csc_data(dim: usize, row_indices: Vec<usize>, values: Vec<R>) -> Option<SpVec<R>> {
        let col_offsets = vec![0, row_indices.len()];
        let csc = CscMatrix::try_from_csc_data(dim, 1, col_offsets, row_indices, values).ok()?;
        Some(SpMat::from(csc).into_spvec())
    }
    
    pub fn from_entries<T>(dim: usize, entries: T) -> Self
    where T: IntoIterator<Item = (usize, R)> {
        SpMat::from_entries(
            (dim, 1), 
            entries.into_iter().map(|(i, a)| (i, 0, a))
        ).into_spvec()
    }

    pub fn from_sorted_entries<T>(dim: usize, entries: T) -> Self
    where T: IntoIterator<Item = (usize, R)> {
        let init = (vec![], vec![]);
        let (row_indices, values) = entries.into_iter().fold(init, |mut res, (i, a)| { 
            assert!(i < dim);
            res.0.push(i);
            res.1.push(a);
            res
        });
        Self::try_from_csc_data(dim, row_indices, values).unwrap()
    }

    pub fn extract<F>(&self, dim: usize, f: F) -> SpVec<R>
    where F: Fn(usize) -> Option<usize> { 
        SpVec::from_entries(dim, self.iter().filter_map(|(i, a)|
            f(i).map(|i| (i, a.clone()))
        ))
    }

    pub fn permute(&self, p: &Perm) -> SpVec<R> {
        self.extract(self.dim(), |i| Some(p.at(i)))
    }

    pub fn subvec(&self, range: Range<usize>) -> SpVec<R> { 
        self.extract(
            range.end - range.start, 
            |i| range.contains(&i).then(|| i - range.start)
        )
    }

    pub fn stack(top: Self, bot: Self) -> SpVec<R> {
        let (n1, n2) = (top.dim(), bot.dim());
        let (_, mut rows, mut vals) = top.inner.disassemble();
        let (_, bot_rows, bot_vals) = bot.inner.disassemble();

        rows.extend(bot_rows.into_iter().map(|i| i + n1));
        vals.extend(bot_vals);

        Self::try_from_csc_data(n1 + n2, rows, vals).unwrap()
    }

    pub fn split(self, at: usize) -> (SpVec<R>, SpVec<R>) {
        let n = self.dim();
        assert!(at <= n);

        let (_, mut rows, mut vals) = self.inner.disassemble();
        let split_idx = rows.partition_point(|&i| i < at);

        let bot_rows: Vec<usize> = rows.split_off(split_idx).into_iter().map(|i| i - at).collect();
        let bot_vals = vals.split_off(split_idx);

        let top = Self::try_from_csc_data(at, rows, vals).unwrap();
        let bot = Self::try_from_csc_data(n - at, bot_rows, bot_vals).unwrap();
        (top, bot)
    }

}

impl<R> Default for SpVec<R> {
    fn default() -> Self {
        Self::zero(0)
    }
}

impl<R: PartialEq + Zero> PartialEq for SpVec<R> {
    fn eq(&self, other: &Self) -> bool {
        self.dim() == other.dim() && self.iter_nz().eq(other.iter_nz())
    }
}

impl<R: Eq + Zero> Eq for SpVec<R> {}

impl<R> Neg for SpVec<R>
where R: AddGrp, for<'a> &'a R: AddGrpOps<R> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        SpVec { inner: -self.inner }
    }
}

impl<R> Neg for &SpVec<R>
where R: Scalar + Neg<Output = R> {
    type Output = SpVec<R>;
    fn neg(self) -> Self::Output {
        SpVec { inner: -&self.inner }
    }
}

macro_rules! impl_binop {
    ($trait:ident, $method:ident) => {
        #[auto_ops]
        impl<'a, 'b, R> $trait<&'b SpVec<R>> for &'a SpVec<R>
        where R: Scalar + ClosedAddAssign + ClosedSubAssign + ClosedMulAssign + Zero + One + Neg<Output = R> {
            type Output = SpVec<R>;
            fn $method(self, rhs: &'b SpVec<R>) -> Self::Output {
                let res = (&self.inner).$method(&rhs.inner);
                SpVec::new(res)
            }
        }
    };
}

impl_binop!(Add, add);
impl_binop!(Sub, sub);

// SpMat * SpVec
#[auto_ops(val_val, val_ref, ref_val)]
impl<'b, R> Mul<&'b SpVec<R>> for &SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = SpVec<R>;
    fn mul(self, rhs: &'b SpVec<R>) -> Self::Output {
        let res = self.inner() * &rhs.inner;
        SpVec::new(res)
    }
}

impl<R> Display for SpVec<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.inner.fmt(f)
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use super::*;

    #[test]
    fn from_vec() {
        let v = SpVec::from(vec![1,0,3,5,0]);
        assert_eq!(v.inner.disassemble(), (vec![0, 3], vec![0, 2, 3], vec![1, 3, 5]));
    }

    #[test]
    fn from_entries() {
        let v = SpVec::from_entries(5, vec![(0, 1), (4, 5), (2, 3)]);
        assert_eq!(v.inner.disassemble(), (vec![0, 3], vec![0, 2, 4], vec![1, 3, 5]));
    }

    #[test]
    fn to_dense() {
        let v = SpVec::from(vec![1,0,3,5,0]);
        assert_eq!(v.into_dense(), vec![1,0,3,5,0]);
    }

    #[test]
    fn add() {
        let v = SpVec::from(vec![1,0,3,5,0]);
        let w = SpVec::from(vec![2,1,-1,3,2]);
        assert_eq!(v + w, SpVec::from(vec![3,1,2,8,2]));
    }

    #[test]
    fn sub() {
        let v = SpVec::from(vec![1,0,3,5,0]);
        let w = SpVec::from(vec![2,1,-1,3,2]);
        assert_eq!(v - w, SpVec::from(vec![-1,-1,4,2,-2]));
    }

    #[test]
    fn neg() {
        let v = SpVec::from(vec![1,0,3,5,0]);
        assert_eq!(-v, SpVec::from(vec![-1,0,-3,-5,0]));
    }

    #[test]
    fn subvec() {
        let v = SpVec::from((0..10).collect_vec());
        let w = v.subvec(3..7);
        assert_eq!(w, SpVec::from(vec![3,4,5,6]))
    }

    #[test]
    fn subvec2() {
        let v = SpVec::from((0..10).collect_vec());
        let w = v.subvec(1..9);
        let w = w.subvec(1..4);
        assert_eq!(w, SpVec::from(vec![2,3,4]))
    }

    #[test]
    fn permute() {
        let p = Perm::new(vec![1,3,0,2]);
        let v = SpVec::from(vec![0,1,2,3]);
        let w = v.permute(&p);
        assert_eq!(w, SpVec::from(vec![2,0,3,1]));
    }

    #[test]
    fn stack() {
        let v1 = SpVec::from((0..3).collect_vec());
        let v2 = SpVec::from((5..8).collect_vec());
        let w = SpVec::stack(v1, v2);
        assert_eq!(w, SpVec::from(vec![0,1,2,5,6,7]));
    }

    #[test]
    fn split() { 
        let v = SpVec::from((0..10).collect_vec());
        let (x, y) = v.split(4);
        assert_eq!(x, SpVec::from((0..4).collect_vec()));
        assert_eq!(y, SpVec::from((4..10).collect_vec()));
    }
}