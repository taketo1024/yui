use sprs::PermView;
use yui_core::{RingOps, Ring};
use crate::sparse::{SpMat, MatType, SpVec};

#[derive(Clone, Debug)]
pub struct Trans<R> 
where R: Ring, for <'x> &'x R: RingOps<R> {
    f: SpMat<R>,
    b: SpMat<R>
}

impl<R> Trans<R> 
where R: Ring, for <'x> &'x R: RingOps<R> { 
    pub fn new(f: SpMat<R>, b: SpMat<R>) -> Self {
        assert!(f.rows() == b.cols());
        assert!(f.cols() == b.rows());
        Self { f, b }
    }

    pub fn id(n: usize) -> Self { 
        Self::new(SpMat::id(n), SpMat::id(n))
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

    pub fn f_mat(&self) -> &SpMat<R> {
        &self.f
    }

    pub fn b_mat(&self) -> &SpMat<R> {
        &self.b
    }

    pub fn trans_f(&self, v: &SpVec<R>) -> SpVec<R> {
        &self.f * v
    }

    pub fn trans_b(&self, v: &SpVec<R>) -> SpVec<R> {
        &self.b * v
    }

    pub fn map<F, B>(&self, f_map: F, b_map: B) -> Self
    where 
        F: FnOnce(&SpMat<R>) -> SpMat<R>,
        B: FnOnce(&SpMat<R>) -> SpMat<R> 
    {
        Self::new(
            f_map(&self.f), 
            b_map(&self.f)
        ) 
    }

    pub fn compose(&self, f: &SpMat<R>, b: &SpMat<R>) -> Self {
        assert!(f.cols() == self.tgt_dim());
        assert!(b.rows() == self.tgt_dim());
        assert!(f.rows() == b.cols());

        self.map(
            |f0|  f * f0, 
            |b0| b0 * b
        )
    }

    pub fn compose_perm(&self, p: PermView) -> Self { 
        self.map(
            |f| f.permute_rows(p.clone()), 
            |b| b.permute_cols(p.clone())
        )
    }
}
