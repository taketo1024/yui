use std::fmt::Display;
use std::ops::Index;
use sprs::CsMat;

use crate::math::traits::{Ring, RingOps, EucRing, EucRingOps};
use crate::math::matrix::{snf_in_place, DnsMat};
use crate::math::matrix::sparse::*;
use super::base::{GradedRModStr, RModStr, GenericRModStr, RModGrid, AdditiveIndex, AdditiveIndexRange};
use super::complex::ChainComplex;

pub trait HomologyComputable<S>: ChainComplex
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>,
    S: RModStr<R = Self::R>
{
    fn homology_at(&self, i: Self::Index) -> S;
}

impl<R, C> HomologyComputable<GenericRModStr<R>> for C
where 
    R: EucRing + CsMatElem, for<'x> &'x R: EucRingOps<R>,
    C: ChainComplex<R = R>,
    C::Output: RModStr<R = R>
{ 
    fn homology_at(&self, k: C::Index) -> GenericRModStr<Self::R> {
        let d1 = self.d_matrix(k - self.d_degree());
        let d2 = self.d_matrix(k);
        let (rank, tors) = compute_homology(&d1, &d2, false);
        GenericRModStr::new(rank, tors)
    }
}

pub trait Homology: GradedRModStr
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>
{
    fn is_zero(&self) -> bool {
        self.range().all(|i| self[i].is_zero())
    }

    fn is_free(&self) -> bool {
        self.range().all(|i| self[i].is_free())
    }

    fn fmt_default(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in self.range() { 
            write!(f, "H[{}]: {}\n", i, self[i])?
        }
        Ok(())
    }
}

pub struct GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    grid: RModGrid<GenericRModStr<R>, I>
}

impl<R, C> From<C> for GenericHomology<R, C::IndexRange>
where
    R: Ring, for<'x> &'x R: RingOps<R>,
    C: HomologyComputable<GenericRModStr<R>, R = R>,
    C::IndexRange: Clone,
    C::Output: RModStr<R = C::R>
{
    fn from(c: C) -> Self {
        let range = c.range();
        let grid = RModGrid::new(range, |i| {
            let h_i = c.homology_at(i);
            if !h_i.is_zero() { Some(h_i) } else { None }
        });
        Self { grid }
    }
}

impl<R, I> GradedRModStr for GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    type R = R;
    type Index = I::Item;
    type IndexRange = I;

    fn in_range(&self, k: Self::Index) -> bool {
        self.grid.in_range(k)
    }

    fn range(&self) -> Self::IndexRange {
        self.grid.range()
    }
}

impl<R, I> Index<I::Item> for GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
    type Output = GenericRModStr<R>;

    fn index(&self, k: I::Item) -> &Self::Output {
        &self.grid[k]
    }
}

impl<R, I> Homology for GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
{
}

impl<R, I> Display for GenericHomology<R, I>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: AdditiveIndexRange,
    I::Item: AdditiveIndex
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

pub fn compute_homology<R>(d1: &CsMat<R>, d2: &CsMat<R>, with_trans: bool) -> (usize, Vec<R>)
where R: EucRing + CsMatElem, for<'x> &'x R: EucRingOps<R> {
    assert_eq!(d2.cols(), d1.rows());
    debug_assert!((d2 * d1).is_zero());

    if d1.rows() == 0 { 
        return (0, vec![])
    }

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

    (rank, tors)
}

#[cfg(test)]
mod tests { 
    use sprs::CsMat;
    use crate::math::matrix::sparse::CsMatExt;
    use super::*;
    use super::super::complex::tests::*;

    #[test]
    fn cancel_pair() { 
        let c = TestChainComplex::<i32>::descending(
            vec![ CsMat::csc_from_vec((1, 1), vec![1]) ],
        );

        let h = GenericHomology::from(c);
        
        assert_eq!(h[0].rank(), 0);
        assert!(h[0].is_free());
    }

    #[test]
    fn torsion() { 
        let c = TestChainComplex::<i32>::descending( 
            vec![ CsMat::csc_from_vec((1, 1), vec![2]) ],
        );

        let h = GenericHomology::from(c);
        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[0].tors(), &vec![2]);
    }

    #[test]
    fn d3() {
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
    fn s2() {
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
    fn t2() {
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
    fn rp2() {
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

    #[test]
    fn d3_dual() {
        let c = TestChainComplex::<i32>::d3().dual();
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
    fn s2_dual() {
        let c = TestChainComplex::<i32>::s2().dual();
        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn t2_dual() {
        let c = TestChainComplex::<i32>::t2().dual();
        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn rp2_dual() {
        let c = TestChainComplex::<i32>::rp2().dual();
        let h = GenericHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[2].tors(), &vec![2]);
        assert_eq!(h[2].is_free(), false);
    }
}