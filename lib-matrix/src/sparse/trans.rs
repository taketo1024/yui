use yui_core::abst::{Ring, RingOps};
use yui_core::ext::CloneAnd;
use crate::Perm;
use crate::sparse::{SpMat, MatTrait, SpVec};

/// A composable forward/backward sparse linear map, used to track basis
/// changes through a chain of reductions.
///
/// Internally stores two parallel sequences `(f_0, ..., f_k)` and
/// `(b_0, ..., b_k)` such that the forward map is `f_k * ... * f_0` and the
/// backward map is `b_0 * ... * b_k` (so applying `forward` followed by
/// `backward` recovers a vector in the source space).
///
/// New stages are appended with [`append`](Self::append) /
/// [`append_perm`](Self::append_perm); two `Trans`es can be composed via
/// [`merge`](Self::merge).
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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

    /// Single-stage transform with forward map `f` and backward map `b`.
    pub fn new(f: SpMat<R>, b: SpMat<R>) -> Self {
        let mut t = Self::id(f.n_cols());
        t.append(f, b);
        t
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

    pub fn append(&mut self, f: SpMat<R>, b: SpMat<R>) { 
        assert_eq!(f.n_cols(), b.n_rows());
        assert_eq!(f.n_rows(), b.n_cols());
        assert_eq!(f.n_cols(), self.tgt_dim);

        self.tgt_dim = f.n_rows();
        self.f_mats.push(f);
        self.b_mats.push(b);
    }

    pub fn append_perm(&mut self, p: &Perm) {
        assert_eq!(p.len(), self.tgt_dim);
        let f = SpMat::row_perm_mat(p);
        let b = SpMat::col_perm_mat(p);
        self.append(f, b)
    }

    /// Appends `other`'s stages onto `self`; requires `self.tgt_dim == other.src_dim`.
    pub fn merge(&mut self, mut other: Trans<R>) {
        assert_eq!(self.tgt_dim, other.src_dim);

        self.tgt_dim = other.tgt_dim;
        self.f_mats.append(&mut other.f_mats);
        self.b_mats.append(&mut other.b_mats);
    }

    pub fn merged(&self, other: &Trans<R>) -> Self { 
        self.clone_and(|t| 
            t.merge(other.clone())
        )
    }

    pub fn forward_mat(&self) -> SpMat<R> {
        // f = fn * ... f1 * f0
        if self.f_mats.len() == 1 { 
            self.f_mats[0].clone()
        } else { 
            self.f_mats.iter().rev().fold(
                SpMat::id(self.tgt_dim), 
                |res, f| res * f
            )
        }
    }

    pub fn backward_mat(&self) -> SpMat<R> {
        // b = b0 * b1 * ... * bn
        if self.b_mats.len() == 1 { 
            self.b_mats[0].clone()
        } else { 
            self.b_mats.iter().rev().fold(
                SpMat::id(self.tgt_dim), 
                |res, b| b * res
            )
        }
    }

    /// Collapses the stored stages into a single pair of forward/backward matrices.
    pub fn reduce(&mut self) {
        if self.f_mats.len() > 1 { 
            let f = self.forward_mat();
            self.f_mats = vec![f];
        }

        if self.b_mats.len() > 1 { 
            let b = self.backward_mat();
            self.b_mats = vec![b];
        }
    }

    /// Restricts the target to the given index subset, appending an extra
    /// projection / inclusion stage.
    pub fn sub(&self, indices: &[usize]) -> Self {
        let n = self.tgt_dim();
        let p = indices.len();
        let f = SpMat::from_entries(
            (p, n), 
            indices.iter().enumerate().map(|(i, &j)|
                (i, j, R::one())
            )
        );
        let b = SpMat::from_entries(
            (n, p), 
            indices.iter().enumerate().map(|(i, &j)|
                (j, i, R::one())
            )
        );
        self.clone_and(|sub|
            sub.append(f, b)
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
    fn append_perm() {
        let mut t = Trans::<i32>::new(
            SpMat::id(5).submat_rows(0..3),
            SpMat::id(5).submat_cols(0..3),
        );
        t.append_perm(&Perm::from_indices([1,2,0]));

        let v = SpVec::from(vec![0,1,2,3,4]);
        let w = t.forward(&v);
        let x = t.backward(&w);

        assert_eq!(w.into_dense(), vec![2,0,1]);
        assert_eq!(x.into_dense(), vec![0,1,2,0,0]);
    }

    #[test]
    fn is_id() { 
        let t = Trans::<i64>::id(10);
        assert!(t.is_id());
    }
}