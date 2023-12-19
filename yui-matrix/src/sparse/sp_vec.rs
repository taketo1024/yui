use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, Range};
use std::fmt::{Display, Debug};
use nalgebra_sparse::CscMatrix;
use nalgebra_sparse::na::{Scalar, ClosedAdd, ClosedSub, ClosedMul};
use num_traits::{Zero, One};
use sprs::PermView;
use auto_impl_ops::auto_ops;
use yui::{Ring, RingOps, AddGrpOps,  AddGrp};
use super::sp_mat::SpMat;

#[derive(Clone, Debug, PartialEq, Eq)]
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

    pub fn view(&self) -> SpVecView<R> { 
        SpVecView::new(self, self.dim(), |i| Some(i))
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, &R)> { 
        self.inner.triplet_iter().map(|(i, _, a)| (i, a))
    }

    pub fn iter_nz(&self) -> impl Iterator<Item = (usize, &R)>
    where R: Zero { 
        self.iter().filter(|(_, a)| !a.is_zero())
    }

    pub fn into_vec(self) -> Vec<R>
    where R: Clone + Zero { 
        self.into()
    }

    pub fn into_mat(self) -> SpMat<R> { 
        self.into()
    }
}

impl<R> From<Vec<R>> for SpVec<R>
where R: Scalar + Zero + ClosedAdd {
    fn from(vec: Vec<R>) -> Self {
        Self::from_entries(vec.len(), vec.into_iter().enumerate())
    }
}

impl<R> From<SpVec<R>> for Vec<R>
where R: Clone + Zero {
    fn from(value: SpVec<R>) -> Self {
        let mut res = vec![R::zero(); value.dim()];
        for (i, a) in value.iter_nz() { 
            res[i] = a.clone();
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
where R: Scalar + Zero + ClosedAdd { 
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
        Self::from_raw_data(dim, row_indices, values)
    }

    fn from_raw_data(dim: usize, row_indices: Vec<usize>, values: Vec<R>) -> SpVec<R> { 
        let col_offsets = vec![0, row_indices.len()];
        let csc = CscMatrix::try_from_csc_data(dim, 1, col_offsets, row_indices, values).unwrap();
        SpMat::from(csc).into_spvec()
    }
    
    pub fn stack_vecs<I>(vecs: I) -> Self 
    where I: IntoIterator<Item = SpVec<R>> { 
        let init = (0, vec![], vec![]);
        let (dim, row_indices, values) = vecs.into_iter().fold(init, |mut res, v| { 
            let n1 = res.0;
            let n2 = v.dim();
            
            let (_, mut rows, mut vals) = v.inner.disassemble();
            rows.iter_mut().for_each(|i| *i += n1);

            res.0 += n2;
            res.1.append(&mut rows);
            res.2.append(&mut vals);
            res
        });
        Self::from_raw_data(dim, row_indices, values)
    }

    pub fn permute(&self, p: PermView<'_>) -> SpVec<R> { 
        self.view().permute(p).collect()
    }

    pub fn subvec(&self, range: Range<usize>) -> SpVec<R> { 
        self.view().subvec(range).collect()
    }

    pub fn stack(&self, other: &SpVec<R>) -> SpVec<R> {
        let (n1, n2) = (self.dim(), other.dim());
        Self::from_entries(n1 + n2, Iterator::chain(
            self.iter_nz().map(|(i, a)| (i, a.clone())),
            other.iter_nz().map(|(i, a)| (n1 + i, a.clone()))
        ))
    }

    pub fn to_dense(&self) -> Vec<R> { 
        let mut vec = vec![R::zero(); self.dim()];
        for (i, a) in self.iter_nz() { 
            vec[i] = a.clone();
        }
        vec
    }
}

impl<R> Default for SpVec<R> {
    fn default() -> Self {
        Self::zero(0)
    }
}

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
        where R: Scalar + ClosedAdd + ClosedSub + ClosedMul + Zero + One + Neg<Output = R> {
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
impl<'a, 'b, R> Mul<&'b SpVec<R>> for &'a SpMat<R>
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

pub struct SpVecView<'a, 'b, R> {
    target: &'a SpVec<R>,
    dim: usize,
    trans: Box<dyn Fn(usize) -> Option<usize> + 'b>
}

impl<'a, 'b, R> SpVecView<'a, 'b, R> {
    fn new<F>(target: &'a SpVec<R>, dim: usize, trans: F) -> Self
    where F: Fn(usize) -> Option<usize> + 'b {
        Self { target, dim, trans: Box::new(trans) }
    }

    pub fn dim(&self) -> usize { 
        self.dim
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, &R)> {
        self.target.iter().filter_map(|(i, a)| { 
            (self.trans)(i).map(|i| (i, a))
        })
    }

    pub fn iter_nz(&self) -> impl Iterator<Item = (usize, &R)>
    where R: Zero {
        self.target.iter_nz().filter_map(|(i, a)| { 
            (self.trans)(i).map(|i| (i, a))
        })
    }

    pub fn permute(&self, p: PermView<'b>) -> SpVecView<R> { 
        SpVecView::new(self.target, self.dim, move |i| Some(p.at(i)))
    }

    pub fn subvec(&self, range: Range<usize>) -> SpVecView<R> { 
        let (i0, i1) = (range.start, range.end);

        assert!(i0 <= i1 && i1 <= self.dim());

        SpVecView::new(self.target, i1 - i0, move |i| { 
            let i = (self.trans)(i)?;
            if range.contains(&i) {
                Some(i - i0)
            } else { 
                None
            }
        })
    }

    pub fn collect(&self) -> SpVec<R>
    where R: Scalar + Zero + ClosedAdd + Clone {
        SpVec::from_entries(self.dim(), self.iter().map(|(i, a)|
            (i, a.clone())
        ))
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use sprs::PermOwned;
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
        assert_eq!(v.to_dense(), vec![1,0,3,5,0]);
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
    fn view() {
        let v = SpVec::from(vec![0,1,2,3]);
        let w = v.view();
        assert_eq!(w.collect(), v);
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
        let p = PermOwned::new(vec![1,3,0,2]);
        let v = SpVec::from(vec![0,1,2,3]);
        let w = v.permute(p.view());
        assert_eq!(w, SpVec::from(vec![2,0,3,1]));
    }

    #[test]
    fn stack() {
        let v1 = SpVec::from((0..3).collect_vec());
        let v2 = SpVec::from((5..8).collect_vec());
        let w = v1.stack(&v2);
        assert_eq!(w, SpVec::from(vec![0,1,2,5,6,7]));
    }
}