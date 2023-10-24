use yui_core::{RingOps, Ring};
use yui_matrix::sparse::{SpMat, MatType, SpVec};

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

    pub fn forward(&self, v: &SpVec<R>) -> SpVec<R> {
        &self.f * v
    }

    pub fn backward(&self, v: &SpVec<R>) -> SpVec<R> {
        &self.b * v
    }

    pub fn compose(&self, f: SpMat<R>, b: SpMat<R>) -> Self {
        assert!(f.rows() == b.cols());
        assert!(f.cols() == b.rows());

        Self::new(
            f * &self.f, 
            &self.b * b
        )
    }
}
