use sprs::{CsMat, PermView, CsVec};
use crate::math::matrix::sparse::{CsMatExt, CsVecExt};
use crate::math::matrix::pivot::{perms_by_pivots, find_pivots_upto};
use crate::math::matrix::schur::schur_partial_upper_triang;
use crate::math::traits::{Ring, RingOps};

//       a0         a1         a2 
// C[0] ----> C[1] ----> C[2] ----> C[3]
//             ::         ::
//             v1         v2

pub struct ChainReducer<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    a0: CsMat<R>, // prev
    a1: CsMat<R>, // target
    a2: CsMat<R>, // next
    v1: Vec<CsVec<R>>,
    v2: Vec<CsVec<R>>,
    step: usize
}

impl<R> ChainReducer<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn reduce(a0: CsMat<R>, a1: CsMat<R>, a2: CsMat<R>)
        -> (CsMat<R>, CsMat<R>, CsMat<R>)
    {
        let (b0, b1, b2, _, _) = Self::reduce_with(a0, a1, a2, vec![], vec![]);
        (b0, b1, b2)
    }

    pub fn reduce_with(a0: CsMat<R>, a1: CsMat<R>, a2: CsMat<R>, v1: Vec<CsVec<R>>, v2: Vec<CsVec<R>>) 
        -> (CsMat<R>, CsMat<R>, CsMat<R>, Vec<CsVec<R>>, Vec<CsVec<R>>) 
    {
        let mut c = Self::new(a0, a1, a2, v1, v2);
        while c.process() {}
        c.result()
    }

    fn new(a0: CsMat<R>, a1: CsMat<R>, a2: CsMat<R>, v1: Vec<CsVec<R>>, v2: Vec<CsVec<R>>) -> Self {
        assert_eq!(a0.rows(), a1.cols());
        assert_eq!(a1.rows(), a2.cols());
        assert!(v1.iter().all(|v| v.dim() == a1.cols()));
        assert!(v2.iter().all(|v| v.dim() == a2.cols()));
        let step = 0;
        Self { a0, a1, a2, v1, v2, step }
    }

    fn result(self) -> (CsMat<R>, CsMat<R>, CsMat<R>, Vec<CsVec<R>>, Vec<CsVec<R>>) {
        (self.a0, self.a1, self.a2, self.v1, self.v2)
    }

    fn set_matrices(&mut self, a0: CsMat<R>, a1: CsMat<R>, a2: CsMat<R>) { 
        (self.a0, self.a1, self.a2) = (a0, a1, a2)
    } 

    fn process(&mut self) -> bool { 
        const MAX_PIVOTS: usize = 300_000;

        let (a0, a1, a2) = (&self.a0, &self.a1, &self.a2);
        let pivs = find_pivots_upto(a1, MAX_PIVOTS);
        let r = pivs.len();
    
        if r == 0 { 
            return false
        }
    
        let (p, q) = perms_by_pivots(a1, &pivs);
        let b1 = a1.permute(p.view(), q.view());
        let sch = schur_partial_upper_triang(b1, r);

        let b1 = sch.complement();
        let b0 = Self::reduce_rows(a0, q.view(), r);
        let b2 = Self::reduce_cols(a2, p.view(), r);

        self.set_matrices(b0, b1, b2);

        // TODO modify vectors

        self.step += 1;
        
        true
    }

    fn reduce_rows(a: &CsMat<R>, p: PermView, r: usize) -> CsMat<R> {
        let (m, n) = a.shape();
        a.permute_rows(p.view()).submatrix(r..m, 0..n)
    }

    fn reduce_cols(a: &CsMat<R>, p: PermView, r: usize) -> CsMat<R> {
        let (m, n) = a.shape();
        a.permute_cols(p.view()).submatrix(0..m, r..n)
    }
}

