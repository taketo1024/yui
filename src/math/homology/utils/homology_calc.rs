use itertools::Itertools;
use sprs::CsMat;

use crate::math::{traits::{EucRing, EucRingOps}, matrix::{DnsMat, snf_in_place, sparse::CsMatExt}};

pub struct HomologyCalc{}

impl HomologyCalc { 
    //         d1                d2
    //  C1 ------------> C2 ------------> C3
    //  |                 |               |
    //  |              p1 |               |
    //  V      d1'        V               |
    //  C1 ------------> C21              |
    //                    ⊕      d2'      |
    //                   C21'----------> C31
    //                    |               |
    //                q2⁻¹|               |
    //                    V      d2''     |
    //                   C22 ----------> C31 
    //                    ⊕
    //                   C23
    // 
    //  H2 = Ker(d2) / Im(d1)
    //     ≅ C23 (free) ⊕ (C21 / Im(d1')) (tor)

    pub fn calculate<R>(d1: CsMat<R>, d2: CsMat<R>, with_trans: bool) -> HomologyCalcResult<R>
    where R: EucRing, for<'x> &'x R: EucRingOps<R> {
        let (n2, _n1) = d1.shape();

        assert_eq!(d2.cols(), n2);

        if n2 == 0 { 
            return HomologyCalcResult::new(0, vec![], None, None) // TODO
        }

        let d1_dns = DnsMat::from(&d1);
        let s1 = snf_in_place(d1_dns, [with_trans, true, false, false]);
        let r1 = s1.rank();

        let d2_dns = if r1 > 0 { 
            let p1_inv = s1.pinv().unwrap().to_sparse();
            let t2 = p1_inv.slice_outer(r1..n2);
            let b2 = &d2 * &t2; // d2': C21' ---> C3
            DnsMat::from(&b2)
        } else {
            DnsMat::from(&d2)
        };

        let s2 = snf_in_place(d2_dns, [false, false, with_trans, with_trans]);
        let r2 = s2.rank();
        let r = n2 - r1 - r2;

        let tors = s1.factors().into_iter().filter_map(|a| {
            if !a.is_unit() {
                Some(a.clone())
            } else {
                None
            }
        }).collect_vec();

        let (p, q) = if with_trans { 
            let t = tors.len();

            let p1 = s1.p().unwrap().to_sparse();         // size = (n2, n2)
            let p11 = p1.submatrix(r1..n2, 0..n2);        // size = (n2 - r1, n2)
            let p2 = s2.qinv().unwrap().to_sparse();      // size = (n2 - r1, n2 - r1)
            let p22 = p2.submatrix(r2..n2-r1, 0..n2-r1);  // size = (n2 - (r1 + r2), n2 - r1)

            let p_free = &p22 * &p11;                     // size = (n2 - (r1 + r2), n2)
            let p_tor = p1.submatrix(r1-t..r1, 0..n2);    // size = (t, n2)
            let p = CsMat::stack(&p_free, &p_tor);        // size = (r + t, n2)

            assert_eq!(p.shape(), (r + t, n2));

            let q1 = s1.pinv().unwrap().to_sparse();      // size = (n2, n2)
            let q12 = q1.submatrix(0..n2, r1..n2);        // size = (n2, n2 - r1)
            let q2 = s2.q().unwrap().to_sparse();         // size = (n2 - r1, n2 - r1)
            let q22 = q2.submatrix(0..n2-r1, r2..n2-r1);  // size = (n2 - r1, n2 - (r1 + r2))

            let q_free = &q12 * &q22;                     // size = (n2, n2 - (r1 + r2))
            let q_tor = q1.submatrix(0..n2, r1-t..r1);    // size = (n2, t)
            let q = CsMat::concat(&q_free, &q_tor);       // size = (n2, r + t)

            assert_eq!(q.shape(), (n2, r + t));

            (Some(p), Some(q))
        } else {
            (None, None)
        };

        HomologyCalcResult::new(r, tors, p, q)
    }
}

pub struct HomologyCalcResult<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    rank: usize,
    tors: Vec<R>,
    p: Option<CsMat<R>>,
    q: Option<CsMat<R>>,
}

impl<R> HomologyCalcResult<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    pub fn new(rank: usize, tors: Vec<R>, p: Option<CsMat<R>>, q: Option<CsMat<R>>) -> Self { 
        Self { rank, tors, p, q }
    }

    pub fn into(self) -> (usize, Vec<R>) { 
        (self.rank, self.tors)
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

        let h = HomologyCalc::calculate(d3, d2, true);
        let p = h.p.unwrap();
        let q = h.q.unwrap();

        dbg!((&p * &q).to_dense());
    }
}