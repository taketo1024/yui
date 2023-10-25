use std::collections::HashMap;
use log::*;
use sprs::PermOwned;

use yui_matrix::sparse::*;
use yui_matrix::sparse::pivot::{perms_by_pivots, find_pivots, PivotType};
use yui_matrix::sparse::schur::Schur;
use yui_core::{Ring, RingOps, Deg};

use super::complex::{ChainComplexTrait, ChainComplexBase};
use super::graded::Graded;

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

pub struct ChainReducer<'a, I, R>
where 
    I: Deg,
    R: Ring, for<'x> &'x R: RingOps<R>,
{ 
    complex: &'a ChainComplexBase<I, R>,
    mats: HashMap<I, SpMat<R>>,
    with_trans: bool,
    trans: Option<HashMap<I, Trans<R>>>
}

impl<'a, I, R> ChainReducer<'a, I, R>
where 
    I: Deg,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    pub fn reduce(complex: &'a ChainComplexBase<I, R>) -> ChainComplexBase<I, R> {
        let mut r = Self::new(complex, false);
        r.process_all();
        r.as_complex()
    }

    pub fn new(complex: &'a ChainComplexBase<I, R>, with_trans: bool) -> Self { 
        let mats = HashMap::new();
        let trans = if with_trans { 
            Some(HashMap::new())
        } else { 
            None
        };
        Self { complex, mats, with_trans, trans }
    }

    pub fn as_complex(mut self) -> ChainComplexBase<I, R> { 
        ChainComplexBase::new(
            self.complex.support(), 
            self.complex.d_deg(), 
            move |i| self.mats.remove(&i).unwrap()
        )
    }

    pub fn matrix(&self, i: I) -> &SpMat<R> {
        self.mats.get(&i).unwrap_or(self.complex.d_matrix(i))
    }

    pub fn trans(&self, i: I) -> Option<&Trans<R>> {
        self.trans.as_ref().and_then(|t| t.get(&i))
    }

    pub fn process_all(&mut self) { 
        for i in self.complex.support() { 
            self.process_at(i)
        }
    }

    pub fn process_at(&mut self, i: I) { 
        self.process_at_itr(i, 0)
    }

    fn process_at_itr(&mut self, i: I, itr: usize) { 
        let a = self.matrix(i);

        info!("reduce at C[{i}] (itr: {itr}), size: {:?}.", a.shape());

        let (p, q, r) = pivots(a);

        if r == 0 { 
            info!("no pivots found.");
            return 
        }

        info!("compute schur complement.");

        let s = schur(&a, &p, &q, r, self.with_trans);

        info!("reduced {:?} -> {:?}.", a.shape(), s.complement().shape());

        if self.with_trans { 
            self.update_trans(i, &p, &q, r, &s);
        }
        self.update_mats(i, &p, &q, r, s);

        // to next iteration
        self.process_at_itr(i, itr + 1)
    }

    fn update_trans(&mut self, i: I, p: &PermOwned, q: &PermOwned, r: usize, s: &Schur<R>) {
        assert!(self.with_trans);

        let (m, n) = self.matrix(i).shape(); // (m, n)
        let (_, i1, i2) = self.deg_trip(i);
        
        if self.trans(i1).is_none() { 
            self.trans.as_mut().unwrap().insert(i1, Trans::id(n));
        }
        if self.trans(i2).is_none() { 
            self.trans.as_mut().unwrap().insert(i2, Trans::id(m));
        }

        let t1 = self.trans(i1).unwrap().permute(q.view()).modify(
            |a| a.submat_rows(r..n),
            |b| b * s.trans_in().unwrap()
        );

        let t2 = self.trans(i2).unwrap().permute(p.view()).modify(
            |a| s.trans_out().unwrap() * a,
            |b| b.submat_cols(r..m)
        );

        let ts = self.trans.as_mut().unwrap();
        ts.insert(i1, t1);
        ts.insert(i2, t2);
    }

    fn update_mats(&mut self, i: I, p: &PermOwned, q: &PermOwned, r: usize, s: Schur<R>) {
        let (i0, i1, i2) = self.deg_trip(i);
        let a0 = self.matrix(i0);
        let a2 = self.matrix(i2);

        let a0 = reduce_mat_rows(a0, &q, r);
        let a1 = s.complement_into();
        let a2 = reduce_mat_cols(a2, &p, r);

        self.mats.insert(i0, a0);
        self.mats.insert(i1, a1);
        self.mats.insert(i2, a2);
    }

    fn deg_trip(&self, i: I) -> (I, I, I) { 
        let deg = self.complex.d_deg();
        (i - deg, i, i + deg)
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
    let m = a.rows();
    a.view().permute_rows(p.view()).submat_rows(r..m).to_owned()
}

fn reduce_mat_cols<R>(a: &SpMat<R>, p: &PermOwned, r: usize) -> SpMat<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    let n = a.cols();
    a.view().permute_cols(p.view()).submat_cols(r..n).to_owned()
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::complex::tests::*;

    #[test]
    fn s2() {
        let c = Samples::<i32>::s2();

        let mut red = ChainReducer::new(&c, false);
        red.process_all();

        let d2 = red.matrix(2);
        let d1 = red.matrix(1);

        assert!( (d1 * d2).is_zero() );
    }

    #[test]
    fn t2() {
        let c = Samples::<i32>::t2();

        let mut red = ChainReducer::new(&c, false);
        red.process_all();

        let d2 = red.matrix(2);
        let d1 = red.matrix(1);

        assert!( (d1 * d2).is_zero() );
    }

    #[test]
    fn rp2() {
        let c = Samples::<i32>::rp2();

        let mut red = ChainReducer::new(&c, false);
        red.process_all();

        let d2 = red.matrix(2);
        let d1 = red.matrix(1);

        assert!( (d1 * d2).is_zero() );
    }

    #[test]
    fn s2_trans() {
        let c = Samples::<i32>::s2();

        let mut red = ChainReducer::new(&c, true);
        red.process_all();

        let t0 = red.trans(0).unwrap();
        let t1 = red.trans(1).unwrap();
        let t2 = red.trans(2).unwrap();

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
        
        assert!(c.is_cycle(2, &w));
    }

    #[test]
    fn t2_trans() {
        let c = Samples::<i32>::t2();

        let mut red = ChainReducer::new(&c, true);
        red.process_all();

        let t0 = red.trans(0).unwrap();
        let t1 = red.trans(1).unwrap();
        let t2 = red.trans(2).unwrap();

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
        
        assert!(c.is_cycle(2, &w));

        let a = SpVec::unit(2, 0);
        let b = SpVec::unit(2, 1);
        let a = t1.backward(&a);
        let b = t1.backward(&b);
        
        assert!(c.is_cycle(1, &a));
        assert!(c.is_cycle(1, &b));
    }

    #[test]
    fn rp2_trans() {
        let c = Samples::<i32>::rp2();

        let mut red = ChainReducer::new(&c, true);
        red.process_all();

        let t0 = red.trans(0).unwrap();
        let t1 = red.trans(1).unwrap();
        let t2 = red.trans(2).unwrap();

        assert_eq!(t0.src_dim(), 6);
        assert_eq!(t1.src_dim(), 15);
        assert_eq!(t2.src_dim(), 10);

        assert_eq!(t0.tgt_dim(), 1);
        assert_eq!(t1.tgt_dim(), 1);
        assert_eq!(t2.tgt_dim(), 1);

        assert_eq!(t0.forward_mat() * t0.backward_mat(), SpMat::id(1));
        assert_eq!(t1.forward_mat() * t1.backward_mat(), SpMat::id(1));
        assert_eq!(t2.forward_mat() * t2.backward_mat(), SpMat::id(1));

        let v = SpVec::unit(1, 0); // generates 2
        let w = t2.backward(&v);
        let dw = c.d(2, &w);
        let dv = t1.forward(&dw);

        assert!(
            dv == SpVec::from(vec![ 2]) || 
            dv == SpVec::from(vec![-2])
        );

        let v = SpVec::unit(1, 0);
        let w = t1.backward(&v);
        
        assert!(c.is_cycle(1, &w));
    }
}