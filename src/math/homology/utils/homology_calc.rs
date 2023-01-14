use itertools::Itertools;
use sprs::CsMat;

use crate::math::{traits::{EucRing, EucRingOps}, matrix::{DnsMat, snf_in_place, sparse::CsMatExt, snf::SnfResult}};

pub type HomologyCalcResult<R> = (usize, Vec<R>, Option<CsMat<R>>, Option<CsMat<R>>);

pub struct HomologyCalc<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    d1: CsMat<R>,
    d2: CsMat<R>,
    with_trans: bool
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

    pub fn calculate(d1: CsMat<R>, d2: CsMat<R>, with_trans: bool) -> HomologyCalcResult<R> {
        let c = Self::new(d1, d2, with_trans);

        let (s1, s2) = c.process();
        let res = c.result(&s1, &s2);

        res
    }

    fn new(d1: CsMat<R>, d2: CsMat<R>, with_trans: bool) -> Self { 
        assert_eq!(d1.rows(), d2.cols());
        Self { d1, d2, with_trans }
    }

    fn process(&self) -> (SnfResult<R>, SnfResult<R>) {
        let with_trans = self.with_trans;
        let (d1, d2) = (&self.d1, &self.d2);
        let n = d2.cols();

        let d1_dns = DnsMat::from(d1);
        let s1 = snf_in_place(d1_dns, [with_trans, true, false, false]);
        let r1 = s1.rank();

        let d2_dns = if r1 > 0 { 
            let p1_inv = s1.pinv().unwrap().to_sparse();
            let t2 = p1_inv.slice_outer(r1..n);
            let d2 = d2 * &t2; // d2': C21' -> C3
            DnsMat::from(&d2)
        } else {
            DnsMat::from(d2)
        };

        let s2 = snf_in_place(d2_dns, [false, false, with_trans, with_trans]);

        (s1, s2)
    }

    fn result(&self, s1: &SnfResult<R>, s2: &SnfResult<R>) -> HomologyCalcResult<R> {
        let n = self.d2.cols();
        let (r1, r2) = (s1.rank(), s2.rank());
        let r = n - r1 - r2;

        let tors = s1.factors().into_iter().filter_map(|a| {
            if !a.is_unit() {
                Some(a.clone())
            } else {
                None
            }
        }).collect_vec();

        let (p, q) = if self.with_trans { 
            let t = tors.len();
            let (p, q) = self.trans(s1, s2, r, t);
            (Some(p), Some(q))
        } else {
            (None, None)
        };

        (r, tors, p, q)
    }

    fn trans(&self, s1: &SnfResult<R>, s2: &SnfResult<R>, r: usize, t: usize) -> (CsMat<R>, CsMat<R>) {
        let n = self.d2.cols();
        let (r1, r2) = (s1.rank(), s2.rank());

        let p1 = s1.p().unwrap().to_sparse();       // size = (n, n)
        let p11 = p1.submatrix(r1..n, 0..n);        // size = (n - r1, n)
        let p2 = s2.qinv().unwrap().to_sparse();    // size = (n - r1, n - r1)
        let p22 = p2.submatrix(r2..n-r1, 0..n-r1);  // size = (n - (r1 + r2), n - r1)

        let p_free = &p22 * &p11;                   // size = (n - (r1 + r2), n)
        let p_tor = p1.submatrix(r1-t..r1, 0..n);   // size = (t, n)
        let p = CsMat::stack(&p_free, &p_tor);      // size = (r + t, n)

        assert_eq!(p.shape(), (r + t, n));

        let q1 = s1.pinv().unwrap().to_sparse();    // size = (n, n)
        let q12 = q1.submatrix(0..n, r1..n);        // size = (n, n - r1)
        let q2 = s2.q().unwrap().to_sparse();       // size = (n - r1, n - r1)
        let q22 = q2.submatrix(0..n-r1, r2..n-r1);  // size = (n - r1, n - (r1 + r2))

        let q_free = &q12 * &q22;                   // size = (n, n - (r1 + r2))
        let q_tor = q1.submatrix(0..n, r1-t..r1);   // size = (n, t)
        let q = CsMat::concat(&q_free, &q_tor);     // size = (n, r + t)

        assert_eq!(q.shape(), (n, r + t));

        (p, q)
    }
}

#[cfg(test)]
mod tests {
    use crate::math::homology::complex::{tests::TestChainComplex, ChainComplex};

    use super::HomologyCalc;
 
    #[test]
    fn trans_s2() {
        let c = TestChainComplex::<i32>::s2();
        let d3 = c.d_matrix(1); // zero
        let d2 = c.d_matrix(0);

        let (r, tors, p, q) = HomologyCalc::calculate(d3, d2, true);
        let p = p.unwrap();
        let q = q.unwrap();

        assert_eq!(r, 1);
        assert_eq!(tors.len(), 0);

        dbg!((&p * &q).to_dense());
    }
}