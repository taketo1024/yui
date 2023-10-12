use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, Range};
use std::fmt::Display;
use num_traits::Zero;
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

    pub fn cs_vec_into(self) -> CsVec<R> { 
        self.cs_vec
    }

    pub fn dim(&self) -> usize { 
        self.cs_vec.dim()
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, &R)> { 
        self.cs_vec.iter().map(|(i, a)| {
            (i, a)
        })
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
        Self::generate(vec.len(), |set| { 
            vec.into_iter().enumerate().for_each(|(i, a)| { 
                set(i, a)
            })
        })
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

macro_rules! _generate {
    ($n:expr, $f:expr) => {{
        let mut ind = Vec::new();
        let mut val = Vec::new();

        $f(&mut |i, a| { 
            if a.is_zero() { return }
            ind.push(i);
            val.push(a);
        });
        
        let Ok(cs_vec) = CsVec::new_from_unsorted($n, ind, val) else { 
            panic!();
        };
        Self::from(cs_vec)
    }};
}

impl<R> SpVec<R> 
where R: Clone + Zero { 
    pub fn generate<F>(n: usize, f: F) -> Self
    where F: FnOnce(&mut (dyn FnMut(usize, R))) { 
        _generate!(n, f)
    }

    pub fn from_entries<T: IntoIterator<Item = (usize, R)>>(dim: usize, iter: T) -> Self {
        Self::generate(dim, |init| { 
            for (i, a) in iter { 
                init(i, a)
            }
        })
    }

    pub fn to_dense(&self) -> Vec<R> { 
        let mut vec = vec![R::zero(); self.dim()];
        for (i, a) in self.iter() { 
            vec[i] = a.clone();
        }
        vec
    }
}

impl<R> SpVec<R> 
where R: Clone + Zero + Send + Sync { 
    pub fn generate_sync<F>(n: usize, f: F) -> Self
    where F: FnOnce(&mut (dyn FnMut(usize, R) + Send + Sync)) { 
        _generate!(n, f)
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

    pub fn permute<'b>(&self, p: PermView<'b>) -> SpVecView<'_, 'b, R> { 
        SpVecView::new(self, self.dim(), move |i| Some(p.at(i)))
    }

    pub fn subvec(&self, range: Range<usize>) -> SpVecView<R> { 
        let (i0, i1) = (range.start, range.end);

        assert!(i0 <= i1 && i1 <= self.dim());

        SpVecView::new(self, i1 - i0, move |i| { 
            if range.contains(&i) {
                Some(i - i0)
            } else { 
                None
            }
        })
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
        SpVec::generate(self.dim(), |set| { 
            for (i, a) in self.iter() { 
                set(i, a.clone())
            }
        })
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
    fn init() {
        let v = SpVec::generate(5, |set| { 
            set(0, 1);
            set(3, 5);
            set(2, 3);
        });
        assert_eq!(v.cs_vec(), &CsVec::new_from_unsorted(5, vec![0, 2, 3], vec![1, 3, 5]).unwrap());
    }

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
        let w = v.subvec(3..7).to_owned();
        assert_eq!(w, SpVec::from(vec![3,4,5,6]))
    }

    #[test]
    fn subvec2() {
        let v = SpVec::from((0..10).collect_vec());
        let w = v.subvec(1..9);
        let w = w.subvec(1..4).to_owned();
        assert_eq!(w, SpVec::from(vec![2,3,4]))
    }

    #[test]
    fn permute() {
        let p = PermOwned::new(vec![1,3,0,2]);
        let v = SpVec::from(vec![0,1,2,3]);
        let w = v.permute(p.view()).to_owned();
        assert_eq!(w, SpVec::from(vec![2,0,3,1]));
    }
}