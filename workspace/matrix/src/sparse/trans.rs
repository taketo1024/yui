use sprs::PermView;
use yui_core::{RingOps, Ring};
use crate::sparse::{SpMat, MatType, SpVec};

#[derive(Clone, Debug)]
pub struct Trans<R> 
where R: Ring, for <'x> &'x R: RingOps<R> {
    f: SpMat<R>,
    b: SpMat<R>,
    f_id: bool,
    b_id: bool
}

impl<R> Trans<R> 
where R: Ring, for <'x> &'x R: RingOps<R> { 
    pub fn new(f: SpMat<R>, b: SpMat<R>) -> Self {
        assert_eq!(f.rows(), b.cols());
        assert_eq!(f.cols(), b.rows());

        let f_id = f.is_id();
        let b_id = b.is_id();
        Self { f, b, f_id, b_id }
    }

    pub fn id(n: usize) -> Self { 
        Self::new(SpMat::id(n), SpMat::id(n))
    }

    pub fn zero() -> Self { 
        Self::id(0)
    }

    pub fn from_perm(p: PermView) -> Self { 
        Self::new(
            SpMat::from_row_perm(p.clone()), 
            SpMat::from_col_perm(p)
        )
    }

    pub fn src_dim(&self) -> usize { 
        self.f.cols()
    }

    pub fn tgt_dim(&self) -> usize { 
        self.f.rows()
    }

    pub fn forward_mat(&self) -> &SpMat<R> {
        &self.f
    }

    pub fn backward_mat(&self) -> &SpMat<R> {
        &self.b
    }

    pub fn forward(&self, v: &SpVec<R>) -> SpVec<R> {
        if self.f_id { 
            v.clone() 
        } else { 
            &self.f * v
        }
    }

    pub fn backward(&self, v: &SpVec<R>) -> SpVec<R> {
        if self.b_id { 
            v.clone()
        } else { 
            &self.b * v
        }
    }

    pub fn is_forward_id(&self) -> bool { 
        self.f_id
    }

    pub fn is_backward_id(&self) -> bool { 
        self.b_id
    }

    pub fn modify<F, B>(&self, f_map: F, b_map: B) -> Self
    where 
        F: FnOnce(&SpMat<R>) -> SpMat<R>,
        B: FnOnce(&SpMat<R>) -> SpMat<R> 
    {
        Self::new(
            f_map(&self.f), 
            b_map(&self.b)
        ) 
    }

    pub fn append(&self, f: &SpMat<R>, b: &SpMat<R>) -> Self {
        assert_eq!(f.cols(), self.tgt_dim());
        assert_eq!(b.rows(), self.tgt_dim());
        assert_eq!(f.rows(), b.cols());

        self.modify(
            |f0|  f * f0, 
            |b0| b0 * b
        )
    }

    pub fn permute(&self, p: PermView) -> Self { 
        assert_eq!(p.dim(), self.f.rows());
        assert_eq!(p.dim(), self.b.cols());

        self.modify(
            |f| f.permute_rows(p.clone()), 
            |b| b.permute_cols(p.clone())
        )
    }

    pub fn compose(&self, other: &Trans<R>) -> Self { 
        assert_eq!(self.tgt_dim(), other.src_dim());
        Self::new(
            other.forward_mat() * self.forward_mat(), 
            self.backward_mat() * other.backward_mat()
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
        let t = Trans::<i32>::new(
            SpMat::id(5).submat_rows(0..3),
            SpMat::id(5).submat_cols(0..3),
        ).permute(
            PermOwned::new(vec![1,2,0]).view()
        );

        let v = SpVec::from(vec![0,1,2,3,4]);
        let w = t.forward(&v);
        let x = t.backward(&w);

        assert_eq!(w, SpVec::from(vec![2,0,1]));
        assert_eq!(x, SpVec::from(vec![0,1,2,0,0]));
    }

    #[test]
    fn is_id() { 
        let t = Trans::<i64>::id(10);
        assert!(t.is_forward_id());
        assert!(t.is_backward_id());
    }
}