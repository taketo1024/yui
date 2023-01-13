use sprs::CsMat;

use crate::math::{traits::{EucRing, EucRingOps}, matrix::{DnsMat, snf_in_place}};

pub struct HomologyCalc{}

impl HomologyCalc { 
    //        d₁            d₂
    //  Rⁿ¹---------> Rᵐ²--------> Rᵖ
    //  |             | P₁         | 
    //  |             V            V
    //  |             *  --------> * ⊕ ..
    //  |             ⊕ 
    //  |             Rᶠ \ 
    //  V      s₁     ⊕   | Z
    //  *  ---------> Rᵇ /
    //  ⊕
    //  :
    //
    //  H = Ker(d₂) / Im(d₁)
    //    ≅ Rᶠ ⊕ (Rᵇ / Im(s₁))

    pub fn calculate<R>(d1: CsMat<R>, d2: CsMat<R>, with_trans: bool) -> HomologyCalcResult<R>
    where R: EucRing, for<'x> &'x R: EucRingOps<R> {
        let (n2, _n1) = d1.shape();

        assert_eq!(d2.cols(), n2);

        if n2 == 0 { 
            return HomologyCalcResult::new(0, vec![])
        }

        let d1_dns = DnsMat::from(&d1);
        let s1 = snf_in_place(d1_dns, [with_trans, true, false, false]);
        let r1 = s1.rank();
        let p1_inv = s1.pinv().unwrap().to_sparse();

        let d2_dns = if r1 > 0 { 
            let t2 = p1_inv.slice_outer(r1..n2);
            DnsMat::from(&(&d2 * &t2))
        } else {
            DnsMat::from(&d2)
        };

        let s2 = snf_in_place(d2_dns, [false, false, false, with_trans]);
        let r2 = s2.rank();

        let rank = n2 - r1 - r2;
        let tors = s1.factors().into_iter().filter_map(|a| {
            if !a.is_unit() {
                Some(a.clone())
            } else {
                None
            }
        }).collect();

        HomologyCalcResult::new(rank, tors)
    }
}

pub struct HomologyCalcResult<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    rank: usize,
    tors: Vec<R>,
}

impl<R> HomologyCalcResult<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    pub fn new(rank: usize, tors: Vec<R>) -> Self { 
        Self { rank, tors }
    }

    pub fn into(self) -> (usize, Vec<R>) { 
        (self.rank, self.tors)
    }
}