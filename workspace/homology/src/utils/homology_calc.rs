use itertools::Itertools;
use log::*;
use yui_core::{EucRing, EucRingOps};
use yui_matrix::dense::{*, snf::*};
use yui_matrix::sparse::*;
use crate::GenericRModStr;

pub type HomologyCalcResult<R> = (usize, Vec<R>, Option<SpMat<R>>, Option<SpMat<R>>);

pub struct HomologyCalc<'a, 'b, R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    d1: &'a SpMat<R>,
    d2: &'b SpMat<R>,
    with_trans: bool
}

impl<'a, 'b, R> HomologyCalc<'a, 'b, R>
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

    pub fn calculate(d1: &'a SpMat<R>, d2: &'b SpMat<R>) -> GenericRModStr<R> {
        let (rank, tor, _) = Self::_calculate(d1, d2, false);
        GenericRModStr::new(rank, tor)
    }

    pub fn calculate_with_trans(d1: &'a SpMat<R>, d2: &'b SpMat<R>) -> (GenericRModStr<R>, SpMat<R>, SpMat<R>) {
        let (rank, tor, trans) = Self::_calculate(d1, d2, true);
        let h = GenericRModStr::new(rank, tor);
        let (p, q) = trans.unwrap();
        (h, p, q)
    }

    fn _calculate(d1: &'a SpMat<R>, d2: &'b SpMat<R>, with_trans: bool) -> (usize, Vec<R>, Option<(SpMat<R>, SpMat<R>)>) {
        info!("calculate homology: {:?}-{:?}", d1.shape(), d2.shape());

        let mut c = Self::new(&d1, &d2, with_trans);
        let (s1, s2) = c.process();
        let (rank, tors) = c.result(&s1, &s2);

        let trans = if with_trans {
            Some(c.trans(&s1, &s2))
        } else {
            None
        };

        (rank, tors, trans)
    }

    fn new(d1: &'a SpMat<R>, d2: &'b SpMat<R>, with_trans: bool) -> Self { 
        assert_eq!(d1.rows(), d2.cols());
        Self { d1, d2, with_trans }
    }

    fn process(&mut self) -> (SnfResult<R>, SnfResult<R>) {
        let with_trans = self.with_trans;
        let (d1, d2) =  (self.d1, self.d2);
        let n = d2.cols();

        let d1_dns = d1.to_dense();
        let s1 = snf_in_place(d1_dns, [with_trans, true, false, false]);
        let r1 = s1.rank();

        let d2_dns = if r1 > 0 { 
            let p1_inv = s1.pinv().unwrap().to_sparse();
            let t2 = p1_inv.submat_cols(r1..n).to_owned();
            let d2 = d2 * &t2; // d2': C21' -> C3
            d2.to_dense()
        } else {
            d2.to_dense()
        };

        let s2 = snf_in_place(d2_dns, [false, false, with_trans, with_trans]);

        (s1, s2)
    }

    fn result(&self, s1: &SnfResult<R>, s2: &SnfResult<R>) -> (usize, Vec<R>) {
        let n = self.d2.cols();
        let (r1, r2) = (s1.rank(), s2.rank());

        let rank = n - r1 - r2;

        let tors = s1.factors().into_iter().filter_map(|a| {
            if !a.is_unit() {
                Some(a.clone())
            } else {
                None
            }
        }).collect_vec();

        (rank, tors)
    }

    fn trans(&self, s1: &SnfResult<R>, s2: &SnfResult<R>) -> (SpMat<R>, SpMat<R>) {
        let n = self.d2.cols();
        let (r1, r2) = (s1.rank(), s2.rank());
        let r = n - r1 - r2;
        let t = s1.factors().iter().filter(|a| !a.is_unit()).count();

        let p1 = s1.p().unwrap().to_sparse();     // size = (n, n)
        let p11 = p1.submat_rows(r1..n)           // size = (n - r1, n)
                    .to_owned();
                
        let p2 = s2.qinv().unwrap().to_sparse();  // size = (n - r1, n - r1)
        let p22 = p2.submat_rows(r2..n-r1)        // size = (n - (r1 + r2), n - r1)
                    .to_owned();

        let p_free = p22 * p11;                   // size = (n - (r1 + r2), n)
        let p_tor = p1.submat_rows(r1-t..r1);     // size = (t, n)
        let p = p_free.stack(&p_tor.to_owned());  // size = (r + t, n)

        assert_eq!(p.shape(), (r + t, n));

        let q1 = s1.pinv().unwrap().to_sparse();  // size = (n, n)
        let q12 = q1.submat_cols(r1..n)           // size = (n, n - r1)
                    .to_owned();

        let q2 = s2.q().unwrap().to_sparse();     // size = (n - r1, n - r1)
        let q22 = q2.submat_cols(r2..n-r1)        // size = (n - r1, n - (r1 + r2))
                    .to_owned();

        let q_free = q12 * q22;                   // size = (n, n - (r1 + r2))
        let q_tor = q1.submat_cols(r1-t..r1);     // size = (n, t)
        let q = q_free.concat(&q_tor.to_owned()); // size = (n, r + t)

        assert_eq!(q.shape(), (n, r + t));

        (p, q)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{RModStr, ChainComplex};
    use crate::test::TestChainComplex;
 
    #[test]
    fn trans_s2_0th() {
        let c = TestChainComplex::<i32>::s2();
        let d1 = c.d_matrix(1);
        let d0 = c.d_matrix(0); // zero

        let (h, p, q) = HomologyCalc::calculate_with_trans(&d1, &d0);

        assert_eq!(h.rank(), 1);
        assert_eq!(h.tors().len(), 0);

        assert!((&p * &q).is_id());

        let v = q.col_vec(0);
        assert_eq!(v.is_zero(), false);
    }

    #[test]
    fn trans_s2_1st() {
        let c = TestChainComplex::<i32>::s2();
        let d1 = c.d_matrix(2);
        let d0 = c.d_matrix(1);

        let (h, p, q) = HomologyCalc::calculate_with_trans(&d1, &d0);

        assert_eq!(h.rank(), 0);
        assert_eq!(h.tors().len(), 0);

        assert!((&p * &q).is_id());
    }

    #[test]
    fn trans_s2_2nd() {
        let c = TestChainComplex::<i32>::s2();
        let d3 = c.d_matrix(3); // zero
        let d2 = c.d_matrix(2);

        let (h, p, q) = HomologyCalc::calculate_with_trans(&d3, &d2);

        assert_eq!(h.rank(), 1);
        assert_eq!(h.tors().len(), 0);

        assert!((&p * &q).is_id());

        let v = q.col_vec(0);
        assert_eq!(v.is_zero(), false);
        assert_eq!((&d2 * &v).is_zero(), true);
    }

    #[test]
    fn trans_t2_0th() {
        let c = TestChainComplex::<i32>::t2();
        let d1 = c.d_matrix(1);
        let d0 = c.d_matrix(0); // zero

        let (h, p, q) = HomologyCalc::calculate_with_trans(&d1, &d0);

        assert_eq!(h.rank(), 1);
        assert_eq!(h.tors().len(), 0);

        assert!((&p * &q).is_id());

        let v = q.col_vec(0);
        assert_eq!(v.is_zero(), false);
    }

    #[test]
    fn trans_t2_1st() {
        let c = TestChainComplex::<i32>::t2();
        let d2 = c.d_matrix(2);
        let d1 = c.d_matrix(1);

        let (h, p, q) = HomologyCalc::calculate_with_trans(&d2, &d1);

        assert_eq!(h.rank(), 2);
        assert_eq!(h.tors().len(), 0);

        assert!((&p * &q).is_id());

        for i in 0..2 { 
            let v = q.col_vec(i);
            assert_eq!(v.is_zero(), false);
            assert_eq!((&d1 * &v).is_zero(), true);
        }
    }

    #[test]
    fn trans_t2_2nd() {
        let c = TestChainComplex::<i32>::t2();
        let d3 = c.d_matrix(3); // zero
        let d2 = c.d_matrix(2);

        let (h, p, q) = HomologyCalc::calculate_with_trans(&d3, &d2);

        assert_eq!(h.rank(), 1);
        assert_eq!(h.tors().len(), 0);

        assert!((&p * &q).is_id());

        let v = q.col_vec(0);
        assert_eq!(v.is_zero(), false);
        assert_eq!((&d2 * &v).is_zero(), true);
    }


    #[test]
    fn trans_rp2_0th() {
        let c = TestChainComplex::<i32>::rp2();
        let d1 = c.d_matrix(1);
        let d0 = c.d_matrix(0); // zero

        let (h, p, q) = HomologyCalc::calculate_with_trans(&d1, &d0);

        assert_eq!(h.rank(), 1);
        assert_eq!(h.tors().len(), 0);

        assert!((&p * &q).is_id());

        let v = q.col_vec(0);
        assert_eq!(v.is_zero(), false);
    }

    #[test]
    fn trans_rp2_1st() {
        let c = TestChainComplex::<i32>::rp2();
        let d2 = c.d_matrix(2);
        let d1 = c.d_matrix(1);

        let (h, p, q) = HomologyCalc::calculate_with_trans(&d2, &d1);

        assert_eq!(h.rank(), 0);
        assert_eq!(h.tors(), &vec![2]);

        assert!((&p * &q).is_id());

        let v = q.col_vec(0);
        assert_eq!(v.is_zero(), false);
        assert_eq!((&d1 * &v).is_zero(), true);
    }
}