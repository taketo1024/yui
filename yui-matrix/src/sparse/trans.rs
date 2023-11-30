use sprs::PermView;
use yui::{RingOps, Ring};
use crate::sparse::{SpMat, MatTrait, SpVec};

#[derive(Clone, Debug)]
pub struct Trans<R> 
where R: Ring, for <'x> &'x R: RingOps<R> {
    src_dim: usize, 
    tgt_dim: usize,
    f_mats: Vec<SpMat<R>>,
    b_mats: Vec<SpMat<R>>,
}

impl<R> Trans<R> 
where R: Ring, for <'x> &'x R: RingOps<R> { 
    pub fn id(n: usize) -> Self { 
        Self { 
            src_dim: n, 
            tgt_dim: n, 
            f_mats: vec![], 
            b_mats: vec![] 
        }
    }

    pub fn zero() -> Self { 
        Self::id(0)
    }

    pub fn new(f: SpMat<R>, b: SpMat<R>) -> Self {
        let mut t = Self::id(f.ncols());
        t.append(f, b);
        t
    }

    pub fn from_perm(p: PermView) -> Self { 
        Self::new(
            SpMat::from_row_perm(p.clone()), 
            SpMat::from_col_perm(p)
        )
    }

    pub fn src_dim(&self) -> usize { 
        self.src_dim
    }

    pub fn tgt_dim(&self) -> usize { 
        self.tgt_dim
    }

    pub fn is_id(&self) -> bool { 
        self.f_mats.is_empty()
    }

    pub fn forward(&self, v: &SpVec<R>) -> SpVec<R> {
        assert_eq!(v.dim(), self.src_dim);
        self.f_mats.iter().fold(v.clone(), |v, f| f * v)
    }

    pub fn backward(&self, v: &SpVec<R>) -> SpVec<R> {
        assert_eq!(v.dim(), self.tgt_dim);
        self.b_mats.iter().rev().fold(v.clone(), |v, f| f * v)
    }

    pub fn permute(&mut self, p: PermView) { 
        assert_eq!(p.dim(), self.tgt_dim);

        let f = self.f_mats.pop().unwrap_or(SpMat::id(self.src_dim));
        let f = f.permute_rows(p.clone());

        let b = self.b_mats.pop().unwrap_or(SpMat::id(self.src_dim));
        let b = b.permute_cols(p);

        self.f_mats.push(f);
        self.b_mats.push(b);
    }

    pub fn append(&mut self, f: SpMat<R>, b: SpMat<R>) { 
        assert_eq!(f.ncols(), b.nrows());
        assert_eq!(f.nrows(), b.ncols());
        assert_eq!(f.ncols(), self.tgt_dim);

        self.tgt_dim = f.nrows();
        self.f_mats.push(f);
        self.b_mats.push(b);
    }

    pub fn merge(&mut self, mut other: Trans<R>) { 
        assert_eq!(self.tgt_dim, other.src_dim);

        self.tgt_dim = other.tgt_dim;
        self.f_mats.append(&mut other.f_mats);
        self.b_mats.append(&mut other.b_mats);
    }

    pub fn modify<F>(&mut self, modify_map: F)
    where F: FnOnce(&mut Vec<SpMat<R>>, &mut Vec<SpMat<R>>) {
        modify_map(&mut self.f_mats, &mut self.b_mats);

        if let Some(f) = self.f_mats.last() { 
            self.tgt_dim = f.nrows();
        } else { 
            self.tgt_dim = self.src_dim;
        }

        if let Some(b) = self.b_mats.last() { 
            assert_eq!(b.ncols(), self.tgt_dim)
        }
    }

    pub fn forward_mat(&self) -> SpMat<R> {
        self.f_mats.iter().fold(
            SpMat::id(self.src_dim), 
            |res, f| f * res
        )
    }

    pub fn backward_mat(&self) -> SpMat<R> {
        self.b_mats.iter().rev().fold(
            SpMat::id(self.tgt_dim), 
            |res, b| b * res
        )
    }
}

impl<R> Default for Trans<R>
where R: Ring, for <'x> &'x R: RingOps<R> {
    fn default() -> Self {
        Self::zero()
    }
}

#[cfg(test)]
mod tests {
    use sprs::PermOwned;

    use super::*;
    use crate::sparse::*;

    #[test]
    fn id() {
        let t = Trans::<i32>::id(5);

        let v = SpVec::from(vec![0,1,2,3,4]);
        let w = t.forward(&v);
        let x = t.backward(&v);

        assert_eq!(w, SpVec::from(vec![0,1,2,3,4]));
        assert_eq!(x, SpVec::from(vec![0,1,2,3,4]));
    }
    
    #[test]
    fn trans() {
        let t = Trans::<i32>::new(
            SpMat::id(5).submat_rows(0..3),
            SpMat::id(5).submat_cols(0..3),
        );

        let v = SpVec::from(vec![0,1,2,3,4]);
        let w = t.forward(&v);
        let x = t.backward(&w);

        assert_eq!(w, SpVec::from(vec![0,1,2]));
        assert_eq!(x, SpVec::from(vec![0,1,2,0,0]));
    }

    #[test]
    fn from_perm() {
        let t = Trans::<i32>::from_perm(
            PermOwned::new(vec![2,3,1,0,4]).view()
        );

        let v = SpVec::from(vec![0,1,2,3,4]);
        let w = t.forward(&v);
        let x = t.backward(&v);

        assert_eq!(w, SpVec::from(vec![3,2,0,1,4]));
        assert_eq!(x, SpVec::from(vec![2,3,1,0,4]));
    }

    #[test]
    fn compose_perm() {
        let mut t = Trans::<i32>::new(
            SpMat::id(5).submat_rows(0..3),
            SpMat::id(5).submat_cols(0..3),
        );
        t.permute(
            PermOwned::new(vec![1,2,0]).view()
        );

        let v = SpVec::from(vec![0,1,2,3,4]);
        let w = t.forward(&v);
        let x = t.backward(&w);

        assert_eq!(w.into_vec(), vec![2,0,1]);
        assert_eq!(x.into_vec(), vec![0,1,2,0,0]);
    }

    #[test]
    fn is_id() { 
        let t = Trans::<i64>::id(10);
        assert!(t.is_id());
    }
}