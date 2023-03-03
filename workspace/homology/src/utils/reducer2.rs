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
    mats: HashMap<C::Idx, SpMat<R>>
}

impl<'a, R, C> ChainReducer<'a, R, C>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    C: ChainComplex<R = R>,
    C::Output: RModStr<R = R>
 { 
    pub fn new(complex: &'a C) -> Self { 
        let mats = HashMap::new();
        Self { complex, mats }
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

        let (b0, b1, b2) = process_triple(a0, a1, a2);

        self.set_matrix(i0, b0);
        self.set_matrix(i1, b1);
        self.set_matrix(i2, b2);
    }

    pub fn take_matrix(&mut self, i: C::Idx) -> SpMat<R> {
        self.mats.remove(&i).unwrap_or(self.complex.d_matrix(i))
    }

    fn set_matrix(&mut self, i: C::Idx, d: SpMat<R>) {
        self.mats.insert(i, d);
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

fn process_triple<R>(a0: SpMat<R>, a1: SpMat<R>, a2: SpMat<R>) -> (SpMat<R>, SpMat<R>, SpMat<R>)
where R: Ring, for<'x> &'x R: RingOps<R> { 
    let (p, q, r) = pivots(&a1);
    if r == 0 { 
        return (a0, a1, a2)
    }

    info!("compute schur complement.");

    let s = schur(&a1, &p, &q, r);

    info!("schur complement: {:?}", s.complement().shape());

    let a0 = reduce_rows(a0, &q, r);
    let a1 = s.complement_into();
    let a2 = reduce_cols(a2, &p, r);
    
    process_triple(a0, a1, a2)
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

fn reduce_rows<R>(a: SpMat<R>, p: &PermOwned, r: usize) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    let m = a.rows();
    a.permute_rows(p.view()).submat_rows(r..m).to_owned()
}

fn reduce_cols<R>(a: SpMat<R>, q: &PermOwned, r: usize) -> SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    let n = a.cols();
    a.permute_cols(q.view()).submat_cols(r..n).to_owned()
}

#[cfg(_test)]
mod tests {
    use super::*;
    use crate::ChainComplex;
    use crate::test::TestChainComplex;
    use yui_matrix::sparse::SpVec;

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