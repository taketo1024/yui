use std::marker::PhantomData;
use log::*;

use yui_core::{EucRing, EucRingOps};
use yui_matrix::dense::{*, snf::*};
use yui_matrix::sparse::*;

use crate::HomologySummand;

pub struct HomologyCalc<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    _r: PhantomData<R>
}

impl<R> HomologyCalc<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    //            d1             d2
    //    C1 ----------> C2 -----------> C3
    //     |              |               |
    //     |           p1 |               |
    //     V      d1'     V               |
    //    C11 ---------> C21              |
    //     ⊕              ⊕      d2'      |
    //    C11'           C21'----------> C3
    //                    |               |
    //                q2⁻¹|               |
    //                    V      d2''     V
    //                   C22 ----------> C31
    //                    ⊕               ⊕
    //                   C22'            C31'
    // 
    //  H2 = Ker(d2) / Im(d1)
    //     ≅ C22' (free) ⊕ (C21 / Im(d1')) (tor)

    pub fn calculate(d1: &SpMat<R>, d2: &SpMat<R>, with_trans: bool) -> HomologySummand<R> {
        info!("calculate homology: {:?}-{:?}", d1.shape(), d2.shape());

        let (s1, s2) = Self::process_snf(d1, d2, with_trans);
        let (rank, tors) = Self::result(&s1, &s2);

        let trans = if with_trans { 
            Some( Self::trans(&s1, &s2) )
        } else {
            None
        };

        HomologySummand::new(rank, tors, trans)
    }

    fn process_snf(d1: &SpMat<R>, d2: &SpMat<R>, with_trans: bool) -> (SnfResult<R>, SnfResult<R>) {
        let n = d1.rows();

        let d1_dns = d1.to_dense();
        let s1 = snf_in_place(d1_dns, [with_trans, true, false, false]);
        let r1 = s1.rank();

        let d2_dns = if r1 > 0 { 
            let p1_inv = s1.pinv().unwrap().to_sparse();
            let t2 = p1_inv.submat_cols(r1..n);
            let d2 = d2 * &t2; // d2': C21' -> C3
            d2.to_dense()
        } else {
            d2.to_dense()
        };

        let s2 = snf_in_place(d2_dns, [false, false, with_trans, with_trans]);

        (s1, s2)
    }

    fn result(s1: &SnfResult<R>, s2: &SnfResult<R>) -> (usize, Vec<R>) {
        let n = s1.result().rows();
        let (r1, r2) = (s1.rank(), s2.rank());
        let rank = n - r1 - r2;

        let tors = s1.factors().into_iter().filter_map(|a| {
            if !a.is_unit() {
                Some(a.clone())
            } else {
                None
            }
        }).collect();

        (rank, tors)
    }

    fn trans(s1: &SnfResult<R>, s2: &SnfResult<R>) -> Trans<R> {
        let n = s1.result().rows();
        let (r1, r2) = (s1.rank(), s2.rank());
        let r = n - r1 - r2;
        let t = s1.factors().iter().filter(|a| !a.is_unit()).count();

        let p1 = s1.p().unwrap().to_sparse();     // size = (n, n)
        let p11 = p1.submat_rows(r1..n);          // size = (n - r1, n)
                
        let p2 = s2.qinv().unwrap().to_sparse();  // size = (n - r1, n - r1)
        let p22 = p2.submat_rows(r2..n-r1);       // size = (n - (r1 + r2), n - r1)

        let p_free = p22 * p11;                   // size = (n - (r1 + r2), n)
        let p_tor = p1.submat_rows(r1-t..r1);     // size = (t, n)
        let p = p_free.stack(&p_tor);             // size = (r + t, n)

        assert_eq!(p.shape(), (r + t, n));

        let q1 = s1.pinv().unwrap().to_sparse();  // size = (n, n)
        let q12 = q1.submat_cols(r1..n);          // size = (n, n - r1)

        let q2 = s2.q().unwrap().to_sparse();     // size = (n - r1, n - r1)
        let q22 = q2.submat_cols(r2..n-r1);       // size = (n - r1, n - (r1 + r2))

        let q_free = q12 * q22;                   // size = (n, n - (r1 + r2))
        let q_tor = q1.submat_cols(r1-t..r1);     // size = (n, t)
        let q = q_free.concat(&q_tor);            // size = (n, r + t)

        assert_eq!(q.shape(), (n, r + t));

        Trans::new(p, q)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ChainComplexTrait;
    use crate::complex::tests::*;
 
    #[test]
    fn s2_0th() {
        let c = Samples::<i32>::s2();
        let d1 = c.d_matrix(1);
        let d0 = c.d_matrix(0); // zero

        let h = HomologyCalc::calculate(&d1, &d0, true);

        assert_eq!(h.rank(), 1);
        assert_eq!(h.tors().len(), 0);

        let v = h.gen(0);

        assert_eq!(v.is_zero(), false);
        assert_eq!(c.is_cycle(0, &v), true);
        assert_eq!(h.trans_forward(&v), SpVec::unit(1, 0));
    }

    #[test]
    fn s2_1st() {
        let c = Samples::<i32>::s2();
        let d1 = c.d_matrix(2);
        let d0 = c.d_matrix(1);

        let h = HomologyCalc::calculate(&d1, &d0, true);

        assert_eq!(h.is_zero(), true);
    }

    #[test]
    fn s2_2nd() {
        let c = Samples::<i32>::s2();
        let d3 = c.d_matrix(3); // zero
        let d2 = c.d_matrix(2);

        let h = HomologyCalc::calculate(&d3, &d2, true);

        assert_eq!(h.rank(), 1);
        assert_eq!(h.tors().len(), 0);

        let v = h.gen(0);

        assert_eq!(v.is_zero(), false);
        assert_eq!(c.is_cycle(2, &v), true);
        assert_eq!(h.trans_forward(&v), SpVec::unit(1, 0));
    }

    #[test]
    fn t2_0th() {
        let c = Samples::<i32>::t2();
        let d1 = c.d_matrix(1);
        let d0 = c.d_matrix(0); // zero

        let h = HomologyCalc::calculate(&d1, &d0, true);

        assert_eq!(h.rank(), 1);
        assert_eq!(h.tors().len(), 0);

        let v = h.gen(0);

        assert_eq!(v.is_zero(), false);
        assert_eq!(c.is_cycle(0, &v), true);
        assert_eq!(h.trans_forward(&v), SpVec::unit(1, 0));
    }

    #[test]
    fn t2_1st() {
        let c = Samples::<i32>::t2();
        let d2 = c.d_matrix(2);
        let d1 = c.d_matrix(1);

        let h = HomologyCalc::calculate(&d2, &d1, true);

        assert_eq!(h.rank(), 2);
        assert_eq!(h.tors().len(), 0);

        for i in 0..2 { 
            let v = h.gen(i);
            assert_eq!(v.is_zero(), false);
            assert_eq!(c.is_cycle(1, &v), true);
            assert_eq!(h.trans_forward(&v), SpVec::unit(2, i));
        }
    }

    #[test]
    fn t2_2nd() {
        let c = Samples::<i32>::t2();
        let d3 = c.d_matrix(3); // zero
        let d2 = c.d_matrix(2);

        let h = HomologyCalc::calculate(&d3, &d2, true);

        assert_eq!(h.rank(), 1);
        assert_eq!(h.tors().len(), 0);

        let v = h.gen(0);

        assert_eq!(v.is_zero(), false);
        assert_eq!(c.is_cycle(2, &v), true);
        assert_eq!(h.trans_forward(&v), SpVec::unit(1, 0));
    }

    #[test]
    fn rp2_0th() {
        let c = Samples::<i32>::rp2();
        let d1 = c.d_matrix(1);
        let d0 = c.d_matrix(0); // zero

        let h = HomologyCalc::calculate(&d1, &d0, true);

        assert_eq!(h.rank(), 1);
        assert_eq!(h.tors().len(), 0);

        let v = h.gen(0);

        assert_eq!(v.is_zero(), false);
        assert_eq!(c.is_cycle(0, &v), true);
        assert_eq!(h.trans_forward(&v), SpVec::unit(1, 0));
    }

    #[test]
    fn rp2_1st() {
        let c = Samples::<i32>::rp2();
        let d2 = c.d_matrix(2);
        let d1 = c.d_matrix(1);

        let h = HomologyCalc::calculate(&d2, &d1, true);

        assert_eq!(h.rank(), 0);
        assert_eq!(h.tors(), &vec![2]);

        let v = h.gen(0);

        assert_eq!(v.is_zero(), false);
        assert_eq!(c.is_cycle(1, &v), true);
        assert_eq!(h.trans_forward(&v), SpVec::unit(1, 0));
    }
}