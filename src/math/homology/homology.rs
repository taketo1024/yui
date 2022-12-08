use std::ops::{Index, RangeInclusive};
use sprs::CsMat;

use crate::math::traits::{Ring, RingOps, EucRing, EucRingOps};
use crate::math::matrix::{snf_in_place, DnsMat};
use crate::math::matrix::sparse::*;
use super::complex::{ChainComplex, ChainComplexSparseD};

pub trait HomologyComputable: ChainComplex
where 
    Self::R: Ring, 
    for<'x> &'x Self::R: RingOps<Self::R>  
{
    fn homology_at(&self, i: isize) -> HomologySummand<Self::R>;
}

impl<C> HomologyComputable for C
where 
    C: ChainComplexSparseD,
    C::R: EucRing + CsMatElem, 
    for<'x> &'x C::R: EucRingOps<C::R>  
{ 
    fn homology_at(&self, k: isize) -> HomologySummand<Self::R> {
        let d1 = self.d_matrix(k - self.d_degree());
        let d2 = self.d_matrix(k);
        compute_homology(&d1, &d2, false)
    }
}

#[derive(Debug, Clone)]
pub struct HomologySummand<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    rank: usize,
    tors: Vec<R>
}

impl<R> HomologySummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(rank: usize, tors: Vec<R>) -> Self { 
        Self { rank, tors }
    }

    pub fn empty() -> Self { 
        Self { rank: 0, tors: vec![] }
    }

    pub fn rank(&self) -> usize { 
        self.rank
    }

    pub fn tors(&self) -> &Vec<R> {
        &self.tors
    }

    pub fn is_free(&self) -> bool { 
        self.tors.is_empty()
    }
}

pub trait Homology: Index<isize>
where Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> {
    type R;
    fn range(&self) -> RangeInclusive<isize>;
}

//       d₁            d₂
// Rⁿ¹---------> Rᵐ²--------> Rᵖ
// |             | P₁         | 
// |             V            V
// |             *  --------> * ⊕ ..
// |             ⊕ 
// |             Rᶠ \ 
// V      s₁     ⊕   | Z
// *  ---------> Rᵇ /
// ⊕
// :
// H = Ker(d₂) / Im(d₁)
// ≅ Rᶠ ⊕ (Rᵇ / Im(s₁))

pub fn compute_homology<R>(d1: &CsMat<R>, d2: &CsMat<R>, with_trans: bool) 
    -> HomologySummand<R>
where 
    R: EucRing + CsMatElem,
    for<'x> &'x R: EucRingOps<R>
{
    assert_eq!(d2.cols(), d1.rows());
    debug_assert!((d2 * d1).is_zero());

    let d1_dns = DnsMat::from(d1);
    let s1 = snf_in_place(d1_dns, [with_trans, true, false, false]);
    let r1 = s1.rank();
    let p1_inv = s1.pinv().unwrap().to_sparse();

    let n2 = d2.cols();
    let d2_dns = if r1 > 0 { 
        let t2 = p1_inv.slice_outer(r1..n2);
        DnsMat::from(&(d2 * &t2))
    } else {
        DnsMat::from(d2)
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

    HomologySummand{ rank, tors }
}