use std::marker::PhantomData;
use log::*;

use yui_core::{EucRing, EucRingOps};
use yui_matrix::dense::{*, snf::*};
use yui_matrix::sparse::*;

use super::homology::HomologySummand;

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

        let n = d2.cols();
        let (s1, s2) = Self::process_snf(d1, d2, with_trans);
        let (rank, tors) = Self::result(n, &s1, &s2);

        if with_trans { 
            let (p, q) = Self::trans(n, &s1, &s2);
            let gens = Self::gens(&q);
            HomologySummand::new(rank, tors, gens, p)
        } else { 
            HomologySummand::new_no_gens(rank, tors)
        }
    }

    fn process_snf(d1: &SpMat<R>, d2: &SpMat<R>, with_trans: bool) -> (SnfResult<R>, SnfResult<R>) {
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

    fn result(n: usize, s1: &SnfResult<R>, s2: &SnfResult<R>) -> (usize, Vec<R>) {
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

    fn trans(n: usize, s1: &SnfResult<R>, s2: &SnfResult<R>) -> (SpMat<R>, SpMat<R>) {
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

    fn gens(q: &SpMat<R>) -> Vec<SpVec<R>> {
        let l = q.shape().1;
        (0..l).map(|j| 
            q.col_vec(j)
        ).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::complex::*;
    use super::super::complex::tests::*;
 
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
        assert_eq!(c.is_cycle(0, v), true);
        assert_eq!(h.vectorize(v), SpVec::unit(1, 0));
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
        assert_eq!(c.is_cycle(2, v), true);
        assert_eq!(h.vectorize(v), SpVec::unit(1, 0));
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
        assert_eq!(c.is_cycle(0, v), true);
        assert_eq!(h.vectorize(v), SpVec::unit(1, 0));
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
            assert_eq!(c.is_cycle(1, v), true);
            assert_eq!(h.vectorize(v), SpVec::unit(2, i));
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
        assert_eq!(c.is_cycle(2, v), true);
        assert_eq!(h.vectorize(v), SpVec::unit(1, 0));
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
        assert_eq!(c.is_cycle(0, v), true);
        assert_eq!(h.vectorize(v), SpVec::unit(1, 0));
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
        assert_eq!(c.is_cycle(1, v), true);
        assert_eq!(h.vectorize(v), SpVec::unit(1, 0));
    }
}