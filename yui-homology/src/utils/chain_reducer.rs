use std::collections::HashMap;
use itertools::Itertools;
use log::*;
use sprs::PermOwned;

use yui_matrix::sparse::*;
use yui_matrix::sparse::pivot::{perms_by_pivots, find_pivots, PivotType};
use yui_matrix::sparse::schur::Schur;
use yui::{Ring, RingOps};

use crate::{GridDeg, ChainComplexTrait, ChainComplexBase, Grid, SimpleRModStr, GridIter, GridTrait};

//       a0 = [x]      a1 = [a b]      a2 = [z w]
//            [y]           [c d]     
//  C[0] --------> C[1] ---------> C[2] -------> C[3]
//    |    [1 a⁻¹b] |             | [ 1     ]    |
//    |    [    1 ] |             | [-ca⁻¹ 1]    |
//    |             V             V              |    
//  C[0] --------> C[1] ---------> C[2] -------> C[3]
//    |     [0]     |    [a 0]    |    [0  w]    |
//    |     [y]     |    [0 s]    |              |
//    |             V             V              |
//  C[0] --------> C[1]'---------> C[2]'-------> C[3]
//          [y]           [s]           [w]

pub struct ChainReducer<I, R>
where 
    I: GridDeg,
    R: Ring, for<'x> &'x R: RingOps<R>,
{ 
    support: Vec<I>,
    d_deg: I,
    with_trans: bool,
    mats: HashMap<I, SpMat<R>>,
    trans: HashMap<I, Trans<R>>
}

impl<I, R> ChainReducer<I, R>
where 
    I: GridDeg,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    pub fn reduce<C>(complex: &C, with_trans: bool) -> ChainComplexBase<I, R> 
    where C: GridTrait<I> + ChainComplexTrait<I, R = R> {
        let r = Self::from(complex, with_trans);
        r.into_complex()
    }

    pub fn from<C>(complex: &C, with_trans: bool) -> Self 
    where C: GridTrait<I> + ChainComplexTrait<I, R = R> {
        let support = complex.support();
        let d_deg = complex.d_deg();

        let mut reducer = Self::new(support, d_deg, with_trans);

        for i in reducer.support.clone() { 
            for j in [i, i + d_deg] { 
                if !reducer.is_set(j) {
                    let d = complex.d_matrix(j);
                    reducer.set_matrix(j, d);
                }
            }
            reducer.reduce_at(i);
        }

        reducer
    }

    pub fn new<Itr>(support: Itr, d_deg: I, with_trans: bool) -> Self
    where Itr: Iterator<Item = I> {
        let support = support.collect_vec();
        let mats = HashMap::new();
        let trans = HashMap::new();
        Self { support, d_deg, with_trans, mats, trans }
    }

    pub fn support(&self) -> GridIter<I> { 
        self.support.clone().into_iter()
    }

    pub fn matrix(&self, i: I) -> Option<&SpMat<R>> {
        self.mats.get(&i)
    }

    pub fn trans(&self, i: I) -> Option<&Trans<R>> {
        self.trans.get(&i)
    }

    pub fn trans_mut(&mut self, i: I) -> Option<&mut Trans<R>> {
        self.trans.get_mut(&i)
    }

    pub fn rank(&self, i: I) -> Option<usize> { 
        self.matrix(i).map(|d| d.ncols())
    }

    pub fn take_matrix(&mut self, i: I) -> Option<SpMat<R>> {
        self.mats.remove(&i)
    }

    pub fn take_trans(&mut self, i: I) -> Option<Trans<R>> {
        self.trans.remove(&i)
    }

    pub fn is_set(&self, i: I) -> bool { 
        self.mats.contains_key(&i)
    }

    pub fn set_matrix(&mut self, i: I, d: SpMat<R>) {
        if self.with_trans { 
            let (m, n) = d.shape();

            if let Some(t) = self.trans(i) {
                assert_eq!(t.tgt_dim(), n)
            } else { 
                self.trans.insert(i, Trans::id(n));
            }

            let i1 = i + self.d_deg;
            if let Some(t) = self.trans(i1) {
                assert_eq!(t.tgt_dim(), m)
            } else { 
                self.trans.insert(i1, Trans::id(m));
            }
        }

        self.mats.insert(i, d);
    }

    pub fn reduce_at(&mut self, i: I) { 
        assert!(self.is_set(i), "not initialized at {i}");
        self.reduce_at_itr(i, 0)
    }

    fn reduce_at_itr(&mut self, i: I, itr: usize) { 
        let a = self.matrix(i).unwrap();
        if a.is_zero() { 
            return;
        }

        info!("red C[{i}]: {:?} (itr: {itr}) ..", a.shape());

        let (p, q, r) = pivots(a);

        if r == 0 { 
            info!("no more pivots.");
            return 
        }

        let s = schur(a, &p, &q, r, self.with_trans);

        info!("red C[{i}]: {:?} -> {:?}.", a.shape(), s.complement().shape());

        if self.with_trans { 
            self.update_trans(i, &p, &q, r, &s);
        }
        self.update_mats(i, &p, &q, r, s);

        // to next iteration
        self.reduce_at_itr(i, itr + 1)
    }

    fn update_trans(&mut self, i: I, p: &PermOwned, q: &PermOwned, r: usize, s: &Schur<R>) {
        assert!(self.with_trans);

        let (m, n) = s.orig_shape(); // (m, n)
        let (_, i1, i2) = self.deg_trip(i);
        
        let t1 = self.trans_mut(i1).unwrap();
        t1.permute(q.view());
        t1.modify(|fs, bs| { 
            let f = fs.pop().unwrap().submat_rows(r..n);
            fs.push(f);

            let b = s.trans_in().unwrap().clone();
            bs.push(b);
        });

        let t2 = self.trans_mut(i2).unwrap();
        t2.permute(p.view());
        t2.modify(|fs, bs| { 
            let f = s.trans_out().unwrap().clone();
            fs.push(f);

            let b = bs.pop().unwrap().submat_cols(r..m);
            bs.push(b)
        });
    }

    fn update_mats(&mut self, i: I, p: &PermOwned, q: &PermOwned, r: usize, s: Schur<R>) {
        let (i0, i1, i2) = self.deg_trip(i);

        if let Some(a0) = self.matrix(i0) {
            let a0 = reduce_mat_rows(a0, q, r);
            self.mats.insert(i0, a0);
        }

        let a1 = s.complement_into();
        self.mats.insert(i1, a1);

        if let Some(a2) = self.matrix(i2) { 
            let a2 = reduce_mat_cols(a2, p, r);
            self.mats.insert(i2, a2);
        }
    }

    fn deg_trip(&self, i: I) -> (I, I, I) { 
        let deg = self.d_deg;
        (i - deg, i, i + deg)
    }

    pub fn into_complex(mut self) -> ChainComplexBase<I, R> { 
        let summands = Grid::generate(
            self.support(), 
            |i| { 
                let rank = self.rank(i).unwrap();
                let trans = self.take_trans(i); // optional
                SimpleRModStr::new( rank, vec![], trans )
            }
        );

        let d_deg = self.d_deg;
        let d_matrices = Grid::generate(
            self.support(),
            |i| self.take_matrix(i).unwrap()
        );

        ChainComplexBase::new(
            summands,
            d_deg,
            d_matrices
        )
    }
}

fn pivots<R>(a: &SpMat<R>) -> (PermOwned, PermOwned, usize) 
where R: Ring, for<'x> &'x R: RingOps<R> {
    let pivs = find_pivots(a, PivotType::Cols);
    let (p, q) = perms_by_pivots(a, &pivs);
    let r = pivs.len();
    (p, q, r)
}

fn schur<R>(a: &SpMat<R>, p: &PermOwned, q: &PermOwned, r: usize, with_trans: bool) -> Schur<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    let b = a.permute(p.view(), q.view());
    Schur::from_partial_lower(&b, r, with_trans)
}

fn reduce_mat_rows<R>(a: &SpMat<R>, p: &PermOwned, r: usize) -> SpMat<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    let m = a.nrows();
    a.view().permute_rows(p.view()).submat_rows(r..m).collect()
}

fn reduce_mat_cols<R>(a: &SpMat<R>, p: &PermOwned, r: usize) -> SpMat<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    let n = a.ncols();
    a.view().permute_cols(p.view()).submat_cols(r..n).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ChainComplex, RModStr, ChainComplexCommon};

    #[test]
    fn zero() { 
        let c = ChainComplex::<i32>::zero();
        let r = ChainReducer::reduce(&c, true);

        r.check_d_all();

        assert!(r[0].is_zero());
    }

    #[test]
    fn acyclic() { 
        let c = ChainComplex::<i32>::one_one(1);
        let r = ChainReducer::reduce(&c, true);

        r.check_d_all();

        assert!(r[0].is_zero());
        assert!(r[1].is_zero());
    }

    #[test]
    fn tor() { 
        let c = ChainComplex::<i32>::one_one(2);
        let r = ChainReducer::reduce(&c, true);

        r.check_d_all();

        assert_eq!(r[0].rank(), 1);
        assert_eq!(r[1].rank(), 1);
    }
    
    #[test]
    fn d3() {
        let c = ChainComplex::<i32>::d3();
        let r = ChainReducer::reduce(&c, false);

        r.check_d_all();

        assert_eq!(r[0].rank(), 1);
        assert_eq!(r[1].rank(), 0);
        assert_eq!(r[2].rank(), 0);
        assert_eq!(r[3].rank(), 0);

        assert!(r.d_matrix(0).is_zero());
        assert!(r.d_matrix(1).is_zero());
        assert!(r.d_matrix(2).is_zero());
        assert!(r.d_matrix(3).is_zero());
    }

    #[test]
    fn s2() {
        let c = ChainComplex::<i32>::s2();
        let r = ChainReducer::reduce(&c, false);

        r.check_d_all();

        assert_eq!(r[0].rank(), 1);
        assert_eq!(r[1].rank(), 0);
        assert_eq!(r[2].rank(), 1);

        assert!(r.d_matrix(0).is_zero());
        assert!(r.d_matrix(1).is_zero());
        assert!(r.d_matrix(2).is_zero());
    }

    #[test]
    fn t2() {
        let c = ChainComplex::<i32>::t2();
        let r = ChainReducer::reduce(&c, false);

        r.check_d_all();

        assert_eq!(r[0].rank(), 1);
        assert_eq!(r[1].rank(), 2);
        assert_eq!(r[2].rank(), 1);

        assert!(r.d_matrix(0).is_zero());
        assert!(r.d_matrix(1).is_zero());
        assert!(r.d_matrix(2).is_zero());
    }

    #[test]
    fn rp2() {
        let c = ChainComplex::<i32>::rp2();
        let r = ChainReducer::reduce(&c, false);

        r.check_d_all();

        assert_eq!(r[0].rank(), 1);
        assert_eq!(r[1].rank(), 1);
        assert_eq!(r[2].rank(), 1);

        assert!( r.d_matrix(0).is_zero());
        assert!( r.d_matrix(1).is_zero());
        assert!(!r.d_matrix(2).is_zero());

        let a = r.d_matrix(2).into_dense()[(0, 0)];
        assert!(a == 2 || a == -2);
    }

    #[test]
    fn s2_trans() {
        let c = ChainComplex::<i32>::s2();
        let r = ChainReducer::reduce(&c, true);

        let t0 = r[0].trans().unwrap();
        let t1 = r[1].trans().unwrap();
        let t2 = r[2].trans().unwrap();

        assert_eq!(t0.src_dim(), 4);
        assert_eq!(t1.src_dim(), 6);
        assert_eq!(t2.src_dim(), 4);

        assert_eq!(t0.tgt_dim(), 1);
        assert_eq!(t1.tgt_dim(), 0);
        assert_eq!(t2.tgt_dim(), 1);

        assert_eq!(t0.forward_mat() * t0.backward_mat(), SpMat::id(1));
        assert_eq!(t1.forward_mat() * t1.backward_mat(), SpMat::id(0));
        assert_eq!(t2.forward_mat() * t2.backward_mat(), SpMat::id(1));

        let v = SpVec::unit(1, 0);
        let w = t2.backward(&v);
        
        assert!(c.d(2, &w).is_zero());
    }

    #[test]
    fn t2_trans() {
        let c = ChainComplex::<i32>::t2();
        let r = ChainReducer::reduce(&c, true);

        let t0 = r[0].trans().unwrap();
        let t1 = r[1].trans().unwrap();
        let t2 = r[2].trans().unwrap();

        assert_eq!(t0.src_dim(), 9);
        assert_eq!(t1.src_dim(), 27);
        assert_eq!(t2.src_dim(), 18);

        assert_eq!(t0.tgt_dim(), 1);
        assert_eq!(t1.tgt_dim(), 2);
        assert_eq!(t2.tgt_dim(), 1);

        assert_eq!(t0.forward_mat() * t0.backward_mat(), SpMat::id(1));
        assert_eq!(t1.forward_mat() * t1.backward_mat(), SpMat::id(2));
        assert_eq!(t2.forward_mat() * t2.backward_mat(), SpMat::id(1));

        let v = SpVec::unit(1, 0);
        let w = t2.backward(&v);
        
        assert!(c.d(2, &w).is_zero());

        let a = SpVec::unit(2, 0);
        let b = SpVec::unit(2, 1);
        let a = t1.backward(&a);
        let b = t1.backward(&b);
        
        assert!(c.d(1, &a).is_zero());
        assert!(c.d(1, &b).is_zero());
    }

    #[test]
    fn rp2_trans() {
        let c = ChainComplex::<i32>::rp2();
        let r = ChainReducer::reduce(&c, true);

        let t0 = r[0].trans().unwrap();
        let t1 = r[1].trans().unwrap();
        let t2 = r[2].trans().unwrap();

        assert_eq!(t0.src_dim(), 6);
        assert_eq!(t1.src_dim(), 15);
        assert_eq!(t2.src_dim(), 10);

        assert_eq!(t0.tgt_dim(), 1);
        assert_eq!(t1.tgt_dim(), 1);
        assert_eq!(t2.tgt_dim(), 1);

        assert_eq!(t0.forward_mat() * t0.backward_mat(), SpMat::id(1));
        assert_eq!(t1.forward_mat() * t1.backward_mat(), SpMat::id(1));
        assert_eq!(t2.forward_mat() * t2.backward_mat(), SpMat::id(1));

        let v = SpVec::unit(1, 0);
        let w = t2.backward(&v);
        let dw = c.d(2, &w);
        let dv = t1.forward(&dw);

        assert_eq!(dv.to_dense()[0].abs(), 2);

        let v = SpVec::unit(1, 0);
        let w = t1.backward(&v);
        
        assert!(c.d(1, &w).is_zero());
    }
}