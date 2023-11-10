use core::panic;
use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, Range};
use std::fmt::Display;
use itertools::Itertools;
use num_traits::{Zero, One};
use sprs::{CsVec, PermView};
use auto_impl_ops::auto_ops;
use yui_core::{Ring, RingOps, AddMonOps, AddGrpOps, AddMon, AddGrp};
use super::SpMat;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SpVec<R> { 
    cs_vec: CsVec<R>
}

impl<R> SpVec<R> { 
    pub fn cs_vec(&self) -> &CsVec<R> { 
        &self.cs_vec
    }

    pub fn dim(&self) -> usize { 
        self.cs_vec.dim()
    }
}

impl<R> From<CsVec<R>> for SpVec<R> {
    fn from(cs_vec: CsVec<R>) -> Self {
        Self { cs_vec }
    }
}

impl<R> From<Vec<R>> for SpVec<R>
where R: Clone + Zero {
    fn from(vec: Vec<R>) -> Self {
        Self::from_entries(vec.len(), vec.into_iter().enumerate())
    }
}

impl<R> SpVec<R> 
where R: Clone + Zero { 
    pub fn generate<F>(dim: usize, init: F) -> Self
    where F: FnOnce(&mut (dyn FnMut(usize, R))) { 
        let mut entries = vec![];
        (init)( &mut |i, a| { 
            entries.push((i, a))
        });
        Self::from_entries(dim, entries)
    }

    pub fn from_entries<T>(dim: usize, entries: T) -> Self
    where T: IntoIterator<Item = (usize, R)> {
        let mut ind = Vec::new();
        let mut val = Vec::new();

        for (i, a) in entries { 
            if a.is_zero() { 
                continue;
            }
            ind.push(i);
            val.push(a);                    
        }
        
        let Ok(cs_vec) = CsVec::new_from_unsorted(dim, ind, val) else { 
            panic!();
        };
        Self::from(cs_vec)
    }

    pub fn reduced(&self) -> Self { 
        Self::from_entries(self.dim(), self.iter().filter_map(|(i, a)| 
            if !a.is_zero() {
                Some((i, a.clone()))
            } else { 
                None
            }
        ))
    }

    pub fn permute<'b>(&self, p: PermView<'b>) -> SpVec<R> { 
        self.view().permute(p).to_owned()
    }

    pub fn subvec(&self, range: Range<usize>) -> SpVec<R> { 
        self.view().subvec(range).to_owned()
    }

    pub fn stack(&self, other: &SpVec<R>) -> SpVec<R> {
        let (n1, n2) = (self.dim(), other.dim());
        Self::from_entries(n1 + n2, Iterator::chain(
            self.iter().map(|(i, a)| (i, a.clone())),
            other.iter().map(|(i, a)| (n1 + i, a.clone()))
        ))
    }

    pub fn to_dense(&self) -> Vec<R> { 
        let mut vec = vec![R::zero(); self.dim()];
        for (i, a) in self.iter() { 
            vec[i] = a.clone();
        }
        vec
    }
}

impl<R> SpVec<R> { 
    pub fn iter(&self) -> impl Iterator<Item = (usize, &R)> { 
        self.cs_vec.iter().map(|(i, a)| {
            (i, a)
        })
    }
}

impl<R> IntoIterator for SpVec<R>
where R: Clone {
    type Item = (usize, R);
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        // MEMO improve this
        self.iter().map(|(i, a)| (i, a.clone())).collect_vec().into_iter()
    }
}

impl<R> Display for SpVec<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let arr = ndarray::Array1::from(self.to_dense());
        arr.fmt(f)
    }
}

impl<R> Default for SpVec<R> {
    fn default() -> Self {
        let cs_vec = CsVec::new(0, vec![], vec![]);
        Self::from(cs_vec)
    }
}

impl<R> SpVec<R> 
where R: Zero { 
    pub fn zero(dim: usize) -> Self { 
        let cs_vec = CsVec::new(dim, vec![], vec![]);
        Self::from(cs_vec)
    }

    pub fn is_zero(&self) -> bool {
        self.cs_vec.data().iter().all(|a| a.is_zero())
    }
}

impl<R> SpVec<R> 
where R: Clone + Zero + One { 
    pub fn unit(n: usize, i: usize) -> Self {
        Self::from_entries(n, [(i, R::one())])
    }
}

impl<R> Neg for SpVec<R>
where R: AddGrp, for<'a> &'a R: AddGrpOps<R> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<R> Neg for &SpVec<R>
where R: AddGrp, for<'a> &'a R: AddGrpOps<R> {
    type Output = SpVec<R>;
    fn neg(self) -> Self::Output {
        let neg = self.cs_vec.map(|a| -a);
        SpVec::from(neg)
    }
}

macro_rules! impl_binop {
    ($trait:ident, $method:ident, $r_trait:ident, $r_op_trait:ident) => {
        #[auto_ops]
        impl<'a, 'b, R> $trait<&'b SpVec<R>> for &'a SpVec<R>
        where R: $r_trait, for<'x> &'x R: $r_op_trait<R> {
            type Output = SpVec<R>;
            fn $method(self, rhs: &'b SpVec<R>) -> Self::Output {
                let res = self.cs_vec().$method(&rhs.cs_vec);
                SpVec::from(res)
            }
        }
    };
}

impl_binop!(Add, add, AddMon, AddMonOps);
impl_binop!(Sub, sub, AddGrp, AddGrpOps);

macro_rules! impl_ops {
    ($trait:ident, $r_trait:ident, $r_op_trait:ident) => {
        impl<R> $trait<SpVec<R>> for SpVec<R>
        where R: $r_trait, for<'x> &'x R: $r_op_trait<R> {}

        impl<R> $trait<SpVec<R>> for &SpVec<R>
        where R: $r_trait, for<'x> &'x R: $r_op_trait<R> {}
    };
}

impl_ops!(AddMonOps, AddMon, AddMonOps);
impl_ops!(AddGrpOps, AddGrp, AddGrpOps);

#[auto_ops(val_val, val_ref, ref_val)]
impl<'a, 'b, R> Mul<&'b SpVec<R>> for &'a SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = SpVec<R>;
    fn mul(self, rhs: &'b SpVec<R>) -> Self::Output {
        let res = self.cs_mat() * rhs.cs_vec();
        SpVec::from(res)
    }
}

impl<R> SpVec<R> {
    pub fn view(&self) -> SpVecView<R> { 
        SpVecView::new(self, self.dim(), |i| Some(i))
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

    pub fn permute(&self, p: PermView<'b>) -> SpVecView<R> { 
        SpVecView::new(self.target, self.dim, move |i| Some(p.at(i)))
    }

    pub fn subvec(&self, range: Range<usize>) -> SpVecView<R> { 
        let (i0, i1) = (range.start, range.end);

        assert!(i0 <= i1 && i1 <= self.dim());

        SpVecView::new(self.target, i1 - i0, move |i| { 
            let Some(i) = (self.trans)(i) else { 
                return None
            };

            if range.contains(&i) {
                Some(i - i0)
            } else { 
                None
            }
        })
    }
}

impl<'a, 'b, R> SpVecView<'a, 'b, R>
where R: Clone + Zero {
    pub fn to_owned(&self) -> SpVec<R> {
        SpVec::from_entries(self.dim(), self.iter().map(|(i, a)|
            (i, a.clone())
        ))
    }

    pub fn to_dense(&self) -> Vec<R> { 
        self.to_owned().to_dense()
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
        assert_eq!(v.cs_vec(), &CsVec::new_from_unsorted(5, vec![0, 2, 3], vec![1, 3, 5]).unwrap());
    }

    #[test]
    fn from_entries() {
        let v = SpVec::from_entries(5, vec![(0, 1), (4, 5), (2, 3)]);
        assert_eq!(v.cs_vec(), &CsVec::new_from_unsorted(5, vec![0, 2, 4], vec![1, 3, 5]).unwrap());
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
        assert_eq!(w.to_owned(), v);
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