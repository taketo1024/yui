use log::*;
use sprs::PermOwned;
use yui_matrix::sparse::*;
use yui_matrix::sparse::pivot::{perms_by_pivots, find_pivots, PivotType};
use yui_matrix::sparse::schur::SchurLT;
use yui_core::{Ring, RingOps};

//          [x]          [a b]
//          [y]          [c d]         [z  w] 
//  C[0] --------> C[1] -------> C[2] -------> C[3]
//    |    [1 a⁻¹b] |             | [a⁻¹    ]    |
//    |    [    1 ] |             | [-ca⁻¹ 1]    |
//    |             V             V              |    
//  C[0] --------> C[1] -------> C[2] -------> C[3]
//    |     [0]     |    [1 0]    |    [0  w]    |
//    |     [y]     |    [0 s]    |              |
//    |             V             V              |
//  C[0] --------> C[1]'-------> C[2]'-------> C[3]
//          [y]           [s]           [w]

pub struct ChainReducer<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    a0: SpMat<R>, // prev
    a1: SpMat<R>, // target
    a2: SpMat<R>, // next
    v1: Vec<SpVec<R>>,
    v2: Vec<SpVec<R>>,
    step: usize
}

impl<R> ChainReducer<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn reduce(a0: SpMat<R>, a1: SpMat<R>, a2: SpMat<R>)
        -> (SpMat<R>, SpMat<R>, SpMat<R>)
    {
        let (b0, b1, b2, _, _) = Self::reduce_with(a0, a1, a2, vec![], vec![]);
        (b0, b1, b2)
    }

    pub fn reduce_with(a0: SpMat<R>, a1: SpMat<R>, a2: SpMat<R>, v1: Vec<SpVec<R>>, v2: Vec<SpVec<R>>) 
        -> (SpMat<R>, SpMat<R>, SpMat<R>, Vec<SpVec<R>>, Vec<SpVec<R>>) 
    {
        let mut c = Self::new(a0, a1, a2, v1, v2);

        info!("reduce: {:?}-{:?}-{:?}", c.a0.shape(), c.a1.shape(), c.a2.shape());
        
        while c.process() {}

        info!("result: {:?}-{:?}-{:?}", c.a0.shape(), c.a1.shape(), c.a2.shape());

        c.result()
    }

    fn new(a0: SpMat<R>, a1: SpMat<R>, a2: SpMat<R>, v1: Vec<SpVec<R>>, v2: Vec<SpVec<R>>) -> Self {
        let (m, n) = a1.shape();

        assert_eq!(a0.rows(), n);
        assert_eq!(a2.cols(), m);
        assert!(v1.iter().all(|v| v.dim() == n));
        assert!(v2.iter().all(|v| v.dim() == m));

        Self { a0, a1, a2, v1, v2, step: 0 }
    }

    fn result(self) -> (SpMat<R>, SpMat<R>, SpMat<R>, Vec<SpVec<R>>, Vec<SpVec<R>>) {
        (self.a0, self.a1, self.a2, self.v1, self.v2)
    }

    fn process(&mut self) -> bool { 
        let (p, q, r) = self.piv_perms();

        if r == 0 { 
            return false
        }
    
        info!("compute schur complement.");

        let s = self.schur(&p, &q, r);

        info!("schur complement: {:?}", s.complement().shape());

        self.reduce_matrices(&p, &q, r);
        self.reduce_vecs(&p, &q, r, &s);
        
        self.a1 = s.complement_into();
        self.step += 1;
        
        true
    }

    fn piv_perms(&self) -> (PermOwned, PermOwned, usize) {
        let a1 = &self.a1;
        let pivs = find_pivots(a1, PivotType::Cols);
        let (p, q) = perms_by_pivots(a1, &pivs);
        let r = pivs.len();

        (p, q, r)
    }

    fn schur(&self, p: &PermOwned, q: &PermOwned, r: usize) -> SchurLT<R> {
        let a1 = &self.a1;
        let b1 = a1.permute(p.view(), q.view());
        SchurLT::from_partial_lower(b1, r)
    }

    fn reduce_matrices(&mut self, p: &PermOwned, q: &PermOwned, r: usize) {
        let (a0, a2) = (&self.a0, &self.a2);

        self.a0 = Self::reduce_rows(a0, q, r);
        self.a2 = Self::reduce_cols(a2, p, r);
    }

    fn reduce_rows(a: &SpMat<R>, p: &PermOwned, r: usize) -> SpMat<R> {
        let m = a.rows();
        a.permute_rows(p.view()).submat_rows(r..m).to_owned()
    }

    fn reduce_cols(a: &SpMat<R>, p: &PermOwned, r: usize) -> SpMat<R> {
        let n = a.cols();
        a.permute_cols(p.view()).submat_cols(r..n).to_owned()
    }

    fn reduce_vecs(&mut self, p: &PermOwned, q: &PermOwned, r: usize, s: &SchurLT<R>) {
        self.reduce_v1(q, r);
        self.reduce_v2(p, s);
    }

    fn reduce_v1(&mut self, q: &PermOwned, r: usize) {
        if self.v1.is_empty() { return }

        let v1 = self.v1.iter().map(|v| {
            let n = v.dim();
            v.permute(q.view()).subvec(r..n).to_owned()
        }).collect();

        self.v1 = v1;
    }

    fn reduce_v2(&mut self, p: &PermOwned, s: &SchurLT<R>) {
        if self.v2.is_empty() { return }

        let v2 = self.v2.iter().map(|v| {
            let v = v.permute(p.view()).to_owned();
            s.trans_vec(v)
        }).collect();
        
        self.v2 = v2;
    }
}

#[cfg(test)]
mod tests {
    use yui_matrix::sparse::SpVec;

    use crate::math::homology::complex::{tests::TestChainComplex, ChainComplex};
    use yui_matrix::sparse::*;

    use super::ChainReducer;

    #[test]
    fn s2_2nd() {
        let c = TestChainComplex::<i32>::s2();
        let d3 = c.d_matrix(3); // zero
        let d2 = c.d_matrix(2);
        let d1 = c.d_matrix(1);

        let v = SpVec::from(vec![-1, 1, -1, 1]);

        let (d3, d2, d1, v_red, _) = ChainReducer::reduce_with(d3, d2, d1, vec![v], vec![]);

        assert_eq!((&d2 * &d3).is_zero(), true);
        assert_eq!((&d1 * &d2).is_zero(), true);

        assert_eq!(v_red.len(), 1);
        assert_eq!(v_red[0].is_zero(), false);
        assert_eq!((&d2 * &v_red[0]).is_zero(), true);
    }

    #[test]
    fn t2_2nd() {
        let c = TestChainComplex::<i32>::t2();
        let d3 = c.d_matrix(3); // zero
        let d2 = c.d_matrix(2);
        let d1 = c.d_matrix(1);

        let v = SpVec::from(vec![1,1,-1,1,-1,-1,1,-1,1,-1,-1,1,1,-1,1,-1,-1,1]);

        let (d3, d2, d1, v_red, _) = ChainReducer::reduce_with(d3, d2, d1, vec![v], vec![]);

        assert_eq!((&d2 * &d3).is_zero(), true);
        assert_eq!((&d1 * &d2).is_zero(), true);
        
        assert_eq!(v_red.len(), 1);
        assert_eq!(v_red[0].is_zero(), false);
        assert_eq!((&d2 * &v_red[0]).is_zero(), true);
    }

    #[test]
    fn t2_1st_right() {
        let c = TestChainComplex::<i32>::t2();
        let d3 = c.d_matrix(3);
        let d2 = c.d_matrix(2);
        let d1 = c.d_matrix(1);

        let v = SpVec::from(vec![0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,-1,0,0,1,-1,1,0,0]);
        let w = SpVec::from(vec![0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,1,0,-1,0,0,0,-1,1]);

        let (d3, d2, d1, _, v_red) = ChainReducer::reduce_with(d3, d2, d1, vec![], vec![v, w]);

        assert_eq!((&d2 * &d3).is_zero(), true);
        assert_eq!((&d1 * &d2).is_zero(), true);
        
        assert_eq!(v_red.len(), 2);
        for i in 0..2 { 
            assert_eq!(v_red[i].is_zero(), false);
            assert_eq!((&d1 * &v_red[i]).is_zero(), true);
        }
    }

    #[test]
    fn t2_1st_left() {
        let c = TestChainComplex::<i32>::t2();
        let d2 = c.d_matrix(2);
        let d1 = c.d_matrix(1);
        let d0 = c.d_matrix(0);

        let v = SpVec::from(vec![0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,-1,0,0,1,-1,1,0,0]);
        let w = SpVec::from(vec![0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,1,0,-1,0,0,0,-1,1]);

        let (d2, d1, d0, v_red, _) = ChainReducer::reduce_with(d2, d1, d0, vec![v, w], vec![]);

        assert_eq!((&d1 * &d2).is_zero(), true);
        assert_eq!((&d0 * &d1).is_zero(), true);
        
        assert_eq!(v_red.len(), 2);
        for i in 0..2 { 
            assert_eq!(v_red[i].is_zero(), false);
            assert_eq!((&d1 * &v_red[i]).is_zero(), true);
        }

        dbg!(v_red);
    }

    #[test]
    fn t2_1st_both() {
        let c = TestChainComplex::<i32>::t2();
        let d3 = c.d_matrix(3); // zero
        let d2 = c.d_matrix(2);
        let d1 = c.d_matrix(1);
        let d0 = c.d_matrix(0);

        let v = SpVec::from(vec![0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,-1,0,0,1,-1,1,0,0]);
        let w = SpVec::from(vec![0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,1,0,-1,0,0,0,-1,1]);

        let (d3, d2, d1, _, v_red) = ChainReducer::reduce_with(d3, d2, d1, vec![], vec![v, w]);
        let (d2, d1, d0, v_red, _) = ChainReducer::reduce_with(d2, d1, d0, v_red, vec![]);

        assert_eq!((&d2 * &d3).is_zero(), true);
        assert_eq!((&d1 * &d2).is_zero(), true);
        assert_eq!((&d0 * &d1).is_zero(), true);
        
        assert_eq!(v_red.len(), 2);

        for i in 0..2 { 
            assert_eq!(v_red[i].is_zero(), false);
            assert_eq!((&d1 * &v_red[i]).is_zero(), true);
        }
    }
}