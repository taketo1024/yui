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
    trans: Option<HashMap<I, SpMat<R>>>,
    trans_rev: Option<HashMap<I, SpMat<R>>>,
}

impl<'a, I, R> ChainReducer<'a, I, R>
where 
    I: Deg,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    pub fn reduce(complex: &'a ChainComplexBase<I, R>) -> ChainComplexBase<I, R> {
        let mut r = Self::new(complex, false);
        r.process();
        r.as_complex()
    }

    fn new(complex: &'a ChainComplexBase<I, R>, with_trans: bool) -> Self { 
        let mats = HashMap::new();
        let (trans, trans_rev) = if with_trans { 
            (Some(HashMap::new()), Some(HashMap::new()))
        } else { 
            (None, None)
        };
        Self { complex, mats, with_trans, trans, trans_rev }
    }

    fn as_complex(mut self) -> ChainComplexBase<I, R> { 
        ChainComplexBase::new(
            self.complex.support(), 
            self.complex.d_deg(), 
            move |i| self.take_matrix(i)
        )
    }

    #[cfg(test)]
    fn matrix(&self, i: I) -> Option<&SpMat<R>> {
        self.mats.get(&i)
    }

    fn take_matrix(&mut self, i: I) -> SpMat<R> {
        self.mats.remove(&i).unwrap_or(self.complex.d_matrix(i).clone())
    }

    fn process(&mut self) { 
        for i in self.complex.support() { 
            self.process_at(i)
        }
    }

    fn process_at(&mut self, i: I) { 
        let c = &self.complex;
        let deg = c.d_deg();

        let (i0, i1, i2) = (i - deg, i, i + deg);

        let a0 = self.take_matrix(i0); // prev
        let a1 = self.take_matrix(i1); // target
        let a2 = self.take_matrix(i2); // next

        info!("i = {i}, reduce: {:?}-{:?}-{:?}", a0.shape(), a1.shape(), a2.shape());

        let (a0, a1, a2) = Self::process_triple(a0, a1, a2);

        info!("result: {:?}-{:?}-{:?}", a0.shape(), a1.shape(), a2.shape());

        self.mats.insert(i0, a0);
        self.mats.insert(i1, a1);
        self.mats.insert(i2, a2);
    }

    fn process_triple(a0: SpMat<R>, a1: SpMat<R>, a2: SpMat<R>) -> (SpMat<R>, SpMat<R>, SpMat<R>) {
        let (p, q, r) = Self::pivots(&a1);
        if r == 0 { 
            return (a0, a1, a2)
        }

        info!("compute schur complement.");

        let s = Self::schur(&a1, &p, &q, r);

        info!("schur complement: {:?}", s.complement().shape());

        let a0 = Self::reduce_mat_rows(a0, &q, r);
        let a1 = s.complement_into();
        let a2 = Self::reduce_mat_cols(a2, &p, r);

        // to next iteration.
        Self::process_triple(a0, a1, a2)
    }

    fn pivots(a: &SpMat<R>) -> (PermOwned, PermOwned, usize) {
        let pivs = find_pivots(a, PivotType::Cols);
        let (p, q) = perms_by_pivots(a, &pivs);
        let r = pivs.len();
        (p, q, r)
    }

    fn schur(a: &SpMat<R>, p: &PermOwned, q: &PermOwned, r: usize) -> Schur<R> {
        let b = a.permute(p.view(), q.view()).to_owned();
        Schur::from_partial_lower(&b, r, false)
    }

    fn reduce_mat_rows(a: SpMat<R>, p: &PermOwned, r: usize) -> SpMat<R> {
        let m = a.rows();
        a.permute_rows(p.view()).submat_rows(r..m).to_owned()
    }

    fn reduce_mat_cols(a: SpMat<R>, p: &PermOwned, r: usize) -> SpMat<R> {
        let n = a.cols();
        a.permute_cols(p.view()).submat_cols(r..n).to_owned()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::complex::tests::*;

    #[test]
    fn s2() {
        let c = Samples::<i32>::s2();

        let mut red = ChainReducer::new(&c, false);
        red.process();

        let d2 = red.matrix(2).unwrap();
        let d1 = red.matrix(1).unwrap();

        assert!( (d1 * d2).is_zero() );
    }

    #[test]
    fn t2() {
        let c = Samples::<i32>::t2();

        let mut red = ChainReducer::new(&c, false);
        red.process();

        let d2 = red.matrix(2).unwrap();
        let d1 = red.matrix(1).unwrap();

        assert!( (d1 * d2).is_zero() );
    }

    #[test]
    fn rp2() {
        let c = Samples::<i32>::rp2();

        let mut red = ChainReducer::new(&c, false);
        red.process();

        let d2 = red.matrix(2).unwrap();
        let d1 = red.matrix(1).unwrap();

        assert!( (d1 * d2).is_zero() );
    }
}