use std::collections::HashMap;
use std::fmt::Display;
use std::ops::Index;
use sprs::CsMat;

use crate::math::traits::{Ring, RingOps, EucRing, EucRingOps};
use crate::math::matrix::{snf_in_place, DnsMat};
use crate::math::matrix::sparse::*;
use super::complex::{ChainComplex, Graded};

pub trait HomologyComputable: ChainComplex
where 
    Self::R: Ring, 
    for<'x> &'x Self::R: RingOps<Self::R>  
{
    fn homology_at(&self, i: Self::Index) -> HomologySummand<Self::R>;
    fn homology<H>(self) -> H where H: Homology<R = Self::R>, H: From<Self>, Self: Sized {
        H::from(self)
    }
}

impl<C> HomologyComputable for C
where 
    C: ChainComplex,
    C::R: EucRing + CsMatElem, 
    for<'x> &'x C::R: EucRingOps<C::R>  
{ 
    fn homology_at(&self, k: C::Index) -> HomologySummand<Self::R> {
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

    pub fn zero() -> Self { 
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

    pub fn is_zero(&self) -> bool { 
        self.rank == 0 && self.tors.is_empty()
    }
}

impl<R> Display for HomologySummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use grouping_by::GroupingBy;
        use crate::utils::format::superscript;

        if self.is_zero() { return write!(f, "0") }

        let mut res = vec![];
        let symbol = R::symbol();
        let r = self.rank as isize;

        if r > 0 {
            let str = if r > 1 { 
                format!("{}{}", symbol, superscript(r))
            }  else {
                format!("{}", symbol)
            };
            res.push(str);
        }

        for (t, r) in self.tors.iter().counter(|&t| t) { 
            let str = if r > 1 { 
                format!("({}/{}){}", symbol, t, superscript(r as isize))
            } else { 
                format!("({}/{})", symbol, t)
            };
            res.push(str);
        }

        write!(f, "{}", res.join(" ⊕ "))
    }
}

pub trait Homology: 
    Graded + 
    Index<Self::Index, Output = HomologySummand<Self::R>>
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> 
{
    type R;
    fn fmt_default(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in self.range() { 
            write!(f, "H[{}]: {}\n", i, self[i])?
        }
        Ok(())
    }
}

pub struct GenericHomology<C>
where 
    C: ChainComplex,
    C::R: Ring + CsMatElem, 
    for<'x> &'x C::R: RingOps<C::R>  
{
    complex: C,
    summands: HashMap<C::Index, HomologySummand<C::R>>,
    zero: HomologySummand<C::R>
}

impl<C> From<C> for GenericHomology<C>
where
    C: HomologyComputable,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>  
{
    fn from(c: C) -> Self {
        let summands = c.range().map(|k| 
            (k, c.homology_at(k))
        ).collect();
        let zero = HomologySummand::zero();
        Self { complex: c, summands, zero }
    }
}

impl<C> Graded for GenericHomology<C>
where
    C: ChainComplex,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>  
{
    type Index = C::Index;
    type IndexRange = C::IndexRange;

    fn in_range(&self, k: Self::Index) -> bool {
        self.complex.in_range(k)
    }

    fn range(&self) -> Self::IndexRange {
        self.complex.range()
    }
}

impl<C> Index<C::Index> for GenericHomology<C>
where
    C: ChainComplex,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>  
{
    type Output = HomologySummand<C::R>;

    fn index(&self, k: C::Index) -> &Self::Output {
        if self.in_range(k) {
            &self.summands[&k]
        } else {
            &self.zero
        }
    }
}

impl<C> Homology for GenericHomology<C>
where
    C: ChainComplex,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>  
{
    type R = C::R;
}

impl<C> Display for GenericHomology<C>
where
    C: ChainComplex,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>  
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt_default(f)
    }
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

#[cfg(test)]
mod tests { 
    use sprs::CsMat;
    use crate::math::matrix::sparse::CsMatExt;
    use super::*;
    use super::super::complex::tests::*;

    #[test]
    fn cancel_pair() { 
        let c = TestChainComplex::<i32>::new(
            -1,
            vec![ CsMat::csc_from_vec((1, 1), vec![1]) ],
        );

        let h = GenericHomology::from(c);
        
        assert_eq!(h[0].rank(), 0);
        assert!(h[0].is_free());
    }

    #[test]
    fn torsion() { 
        let c = TestChainComplex::<i32>::new( 
            -1,
            vec![ CsMat::csc_from_vec((1, 1), vec![2]) ],
        );

        let h = GenericHomology::from(c);
        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[0].tors(), &vec![2]);
    }

    #[test]
    fn homology_d3() {
        let c = TestChainComplex::<i32>::d3();
        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[2].is_free(), true);

        assert_eq!(h[3].rank(), 0);
        assert_eq!(h[3].is_free(), true);
    }

    #[test]
    fn homology_s2() {
        let c = TestChainComplex::<i32>::s2();
        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn homology_t2() {
        let c = TestChainComplex::<i32>::t2();
        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn homology_rp2() {
        let c = TestChainComplex::<i32>::rp2();
        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].tors(), &vec![2]);
        assert_eq!(h[1].is_free(), false);

        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[2].is_free(), true);
    }
}