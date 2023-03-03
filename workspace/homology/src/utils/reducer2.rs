use std::collections::HashMap;
use log::*;
use sprs::PermOwned;
use yui_matrix::sparse::*;
use yui_matrix::sparse::pivot::{perms_by_pivots, find_pivots, PivotType};
use yui_matrix::sparse::schur::SchurLT;
use yui_core::{Ring, RingOps};

use crate::{ChainComplex, RModStr};

pub struct ChainReducer<'a, R, C>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    C: ChainComplex<R = R>,
    C::Output: RModStr<R = R>
 { 
    complex: &'a C,
    mats: HashMap<C::Idx, SpMat<R>>,
    vecs: HashMap<C::Idx, Vec<SpVec<R>>>
}

impl<'a, R, C> ChainReducer<'a, R, C>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    C: ChainComplex<R = R>,
    C::Output: RModStr<R = R>
 { 
    pub fn new(complex: &'a C) -> Self { 
        let mats = HashMap::new();
        let vecs = HashMap::new();
        Self { complex, mats, vecs }
    }

    pub fn set_vec(&mut self, i: C::Idx, v: SpVec<R>) {
        if !self.vecs.contains_key(&i) { 
            self.vecs.insert(i, vec![]);
        }
        self.vecs.get_mut(&i).unwrap().push(v);
    }

    pub fn process(&mut self) { 
        for i in self.complex.indices() { 
            self.process_at(i)
        }
    }

    fn process_at(&mut self, i: C::Idx) { 
        let c = &self.complex;
        let deg = c.d_degree();

        let (i0, i1, i2) = (i - deg, i, i + deg);

        let a0 = self.take_matrix(i0); // prev
        let a1 = self.take_matrix(i1); // target
        let a2 = self.take_matrix(i2); // next

        let v1 = self.take_vecs(i1);
        let v2 = self.take_vecs(i2);

        let (a0, a1, a2, v1, v2) = process_triple(a0, a1, a2, v1, v2);

        self.mats.insert(i0, a0);
        self.mats.insert(i1, a1);
        self.mats.insert(i2, a2);

        self.vecs.insert(i1, v1);
        self.vecs.insert(i2, v2);
    }

    pub fn matrix(&self, i: C::Idx) -> Option<&SpMat<R>> {
        self.mats.get(&i)
    }

    pub fn take_matrix(&mut self, i: C::Idx) -> SpMat<R> {
        self.mats.remove(&i).unwrap_or(self.complex.d_matrix(i))
    }

    pub fn vecs(&self, i: C::Idx) -> Option<&Vec<SpVec<R>>> {
        self.vecs.get(&i)
    }

    pub fn take_vecs(&mut self, i: C::Idx) -> Vec<SpVec<R>> {
        self.vecs.remove(&i).unwrap_or(vec![])
    }
}

//       a0 = [x]      a1 = [a b]      a2 = [z w]
//            [y]           [c d]     
//  C[0] --------> C[1] ---------> C[2] -------> C[3]
//    |    [1 a⁻¹b] |             | [a⁻¹    ]    |
//    |    [    1 ] |             | [-ca⁻¹ 1]    |
//    |             V             V              |    
//  C[0] --------> C[1] ---------> C[2] -------> C[3]
//    |     [0]     |    [1 0]    |    [0  w]    |
//    |     [y]     |    [0 s]    |              |
//    |             V             V              |
//  C[0] --------> C[1]'---------> C[2]'-------> C[3]
//          [y]           [s]           [w]

fn process_triple<R>(a0: SpMat<R>, a1: SpMat<R>, a2: SpMat<R>, v1: Vec<SpVec<R>>, v2: Vec<SpVec<R>>) -> (SpMat<R>, SpMat<R>, SpMat<R>, Vec<SpVec<R>>, Vec<SpVec<R>>)
where R: Ring, for<'x> &'x R: RingOps<R> { 
    let (p, q, r) = pivots(&a1);
    if r == 0 { 
        return (a0, a1, a2, v1, v2)
    }

    info!("compute schur complement.");

    let s = schur(&a1, &p, &q, r);

    info!("schur complement: {:?}", s.complement().shape());

    let v1 = reduce_src_vecs(v1, &q, r);
    let v2 = reduce_tgt_vecs(v2, &p, &s);
    
    let a0 = reduce_mat_rows(a0, &q, r);
    let a1 = s.complement_into();
    let a2 = reduce_mat_cols(a2, &p, r);

    process_triple(a0, a1, a2, v1, v2)
}

fn pivots<R>(a: &SpMat<R>) -> (PermOwned, PermOwned, usize)
where R: Ring, for<'x> &'x R: RingOps<R> { 
    let pivs = find_pivots(a, PivotType::Cols);
    let (p, q) = perms_by_pivots(a, &pivs);
    let r = pivs.len();
    (p, q, r)
}

fn schur<R>(a: &SpMat<R>, p: &PermOwned, q: &PermOwned, r: usize) -> SchurLT<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    let b1 = a.permute(p.view(), q.view());
    SchurLT::from_partial_lower(b1, r)
}

fn reduce_mat_rows<R>(a: SpMat<R>, p: &PermOwned, r: usize) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    let m = a.rows();
    a.permute_rows(p.view()).submat_rows(r..m).to_owned()
}

fn reduce_mat_cols<R>(a: SpMat<R>, p: &PermOwned, r: usize) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    let n = a.cols();
    a.permute_cols(p.view()).submat_cols(r..n).to_owned()
}

fn reduce_src_vecs<R>(vs: Vec<SpVec<R>>, p: &PermOwned, r: usize) -> Vec<SpVec<R>>
where R: Ring, for<'x> &'x R: RingOps<R> {
    vs.into_iter().map(|v| {
        let n = v.dim();
        v.permute(p.view()).subvec(r..n).to_owned()
    }).collect()
}

fn reduce_tgt_vecs<R>(vs: Vec<SpVec<R>>, p: &PermOwned, s: &SchurLT<R>) -> Vec<SpVec<R>>
where R: Ring, for<'x> &'x R: RingOps<R> {
    vs.iter().map(|v| {
        let v = v.permute(p.view()).to_owned();
        s.trans_vec(v)
    }).collect()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::test::TestChainComplex;

    #[test]
    fn s2() {
        let c = TestChainComplex::<i32>::s2();

        let mut red = ChainReducer::new(&c);
        red.process();

        let d2 = red.matrix(2).unwrap();
        let d1 = red.matrix(1).unwrap();

        assert!( (d1 * d2).is_zero() );
    }

    #[test]
    fn t2() {
        let c = TestChainComplex::<i32>::t2();

        let mut red = ChainReducer::new(&c);
        red.process();

        let d2 = red.matrix(2).unwrap();
        let d1 = red.matrix(1).unwrap();

        assert!( (d1 * d2).is_zero() );
    }

    #[test]
    fn rp2() {
        let c = TestChainComplex::<i32>::rp2();

        let mut red = ChainReducer::new(&c);
        red.process();

        let d2 = red.matrix(2).unwrap();
        let d1 = red.matrix(1).unwrap();

        assert!( (d1 * d2).is_zero() );
    }

    #[test]
    fn s2_top() {
        let c = TestChainComplex::<i32>::s2();
        let v = SpVec::from(vec![-1, 1, -1, 1]);

        let mut red = ChainReducer::new(&c);
        red.set_vec(2, v);
        red.process();

        let d2 = red.matrix(2).unwrap();
        let vs = red.vecs(2).unwrap();

        assert_eq!(vs.len(), 1);

        let v = &vs[0];

        assert_eq!(v.is_zero(), false);
        assert!( (d2 * v).is_zero() );
    }

    #[test]
    fn t2_top() {
        let c = TestChainComplex::<i32>::t2();
        let v = SpVec::from(vec![1,1,-1,1,-1,-1,1,-1,1,-1,-1,1,1,-1,1,-1,-1,1]);

        let mut red = ChainReducer::new(&c);
        red.set_vec(2, v);
        red.process();

        let d2 = red.matrix(2).unwrap();
        let vs = red.vecs(2).unwrap();

        assert_eq!(vs.len(), 1);

        let v = &vs[0];

        assert_eq!(v.is_zero(), false);
        assert!( (d2 * v).is_zero() );
    }

    #[test]
    fn t2_1st() {
        let c = TestChainComplex::<i32>::t2();
        let v0 = SpVec::from(vec![0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,-1,0,0,1,-1,1,0,0]);
        let v1 = SpVec::from(vec![0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,1,0,-1,0,0,0,-1,1]);

        let mut red = ChainReducer::new(&c);
        red.set_vec(1, v0);
        red.set_vec(1, v1);
        red.process();

        let d1 = red.matrix(1).unwrap();
        let vs = red.vecs(1).unwrap();

        assert_eq!(vs.len(), 2);

        let v0 = &vs[0];
        let v1 = &vs[1];

        assert_eq!(v0.is_zero(), false);
        assert_eq!(v1.is_zero(), false);
        assert!( (d1 * v0).is_zero() );
        assert!( (d1 * v1).is_zero() );
    }
}