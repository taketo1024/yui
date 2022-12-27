use std::collections::HashMap;
use sprs::{CsMat, PermView};

use crate::math::matrix::pivot::{perms_by_pivots, find_pivots_upto};
use crate::math::matrix::schur::schur_partial_upper_triang;
use crate::math::matrix::sparse::CsMatExt;
use crate::math::traits::{Ring, RingOps};
use crate::math::matrix::CsMatElem;
use super::base::Graded;
use super::complex::ChainComplex;

pub struct Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring + CsMatElem, 
    for<'x> &'x C::R: RingOps<C::R>  
{ 
    d_matrices: HashMap<C::Index, CsMat<C::R>>,
    original: C,
}

impl<C> Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring + CsMatElem, 
    for<'x> &'x C::R: RingOps<C::R>  
{ 
    fn reduce(a0: Option<CsMat<C::R>>, a1: CsMat<C::R>, a2: CsMat<C::R>, k: C::Index, step: usize) -> (Option<CsMat<C::R>>, CsMat<C::R>, CsMat<C::R>) {
        const MAX_PIVOTS: usize = 300_000;

        if let Some(a0) = a0.as_ref() {
            assert_eq!(a0.rows(), a1.cols());
        }
        assert_eq!(a1.rows(), a2.cols());

        let pivs = find_pivots_upto(&a1, MAX_PIVOTS);
        let r = pivs.len();

        if r == 0 { 
            // println!("k = {k} ({step}), no reduce: {}", a1.cols());
            return (a0, a1, a2)
        }
    
        // println!("k = {k} ({step}), reduce: {} -> {}", a1.cols(), a1.cols() - r);
        
        let (p, q) = perms_by_pivots(&a1, &pivs);
        let b1 = a1.permute(p.view(), q.view());
        let b1 = schur_partial_upper_triang(b1, r);
        
        let b0 = a0.map(|a0| {
            Self::reduce_rows(a0, q.view(), r)
        });
        let b2 = Self::reduce_cols(a2, p.view(), r);

        Self::reduce(b0, b1, b2, k, step + 1)
    }

    fn reduce_rows(a: CsMat<C::R>, p: PermView, r: usize) -> CsMat<C::R> {
        let (m, n) = a.shape();
        a.permute_rows(p.view()).submatrix(r..m, 0..n)
    }

    fn reduce_cols(a: CsMat<C::R>, p: PermView, r: usize) -> CsMat<C::R> {
        let (m, n) = a.shape();
        a.permute_cols(p.view()).submatrix(0..m, r..n)
    }
} 

impl<C> From<C> for Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring + CsMatElem, 
    for<'x> &'x C::R: RingOps<C::R>  
{
    fn from(c: C) -> Self {
        let mut d_matrices = HashMap::new();
        for k in c.range() {
            let deg = c.d_degree();
            let (k0, k1, k2) = (k - deg, k, k + deg);

            let a0 = d_matrices.remove(&k0); // prev(optional)
            let a1 = d_matrices.remove(&k1).unwrap_or(c.d_matrix(k1)); // target
            let a2 = c.d_matrix(k2); // next

            let (b0, b1, b2) = Self::reduce(a0, a1, a2, k, 1);
            
            if let Some(b0) = b0 {
                d_matrices.insert(k0, b0);
            }
            d_matrices.insert(k1, b1);
            d_matrices.insert(k2, b2);
        }

        Self { d_matrices, original: c }
    }
}

impl<C> Graded for Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring + CsMatElem, 
    for<'x> &'x C::R: RingOps<C::R>  
{
    type Index = C::Index;
    type IndexRange = C::IndexRange;

    fn in_range(&self, k: Self::Index) -> bool {
        self.original.in_range(k)
    }

    fn range(&self) -> Self::IndexRange {
        self.original.range()
    }
}

impl<C> ChainComplex for Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring + CsMatElem, 
    for<'x> &'x C::R: RingOps<C::R>  
{
    type R = C::R;

    fn rank(&self, k: Self::Index) -> usize {
        if self.in_range(k) {
            self.d_matrices[&k].cols()
        } else {
            0
        }
    }

    fn d_degree(&self) -> Self::Index {
        self.original.d_degree()
    }

    fn d_matrix(&self, k: Self::Index) -> sprs::CsMat<Self::R> {
        if self.in_range(k) { 
            self.d_matrices[&k].clone()
        } else {
            let m = self.rank(k + self.d_degree());
            let n = self.rank(k);
            CsMat::zero((m, n))
        }
    }
}

#[cfg(test)]
mod tests { 
    use crate::math::homology::base::RModStr;
    use crate::math::homology::homology::GenericHomology;

    use super::*;
    use super::super::complex::tests::*;

    #[test]
    fn cancel_pair() { 
        let c = TestChainComplex::<i32>::new(
            -1,
            vec![ CsMat::csc_from_vec((1, 1), vec![1]) ],
        );
        let c = Reduced::from(c);

        assert_eq!(c.rank(0), 0);
        assert_eq!(c.rank(1), 0);
    }

    #[test]
    fn torsion() { 
        let c = TestChainComplex::<i32>::new( 
            -1,
            vec![ CsMat::csc_from_vec((1, 1), vec![2]) ],
        );
        let c = Reduced::from(c);

        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 1);
    }

    #[test]
    fn homology_d3() {
        let c = TestChainComplex::<i32>::d3();
        let c = Reduced::from(c);

        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 0);
        assert_eq!(c.rank(2), 0);
        assert_eq!(c.rank(3), 0);

        let h = GenericHomology::from(&c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[3].rank(), 0);
    }

    #[test]
    fn homology_s2() {
        let c = TestChainComplex::<i32>::s2();
        let c = Reduced::from(c);

        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 0);
        assert_eq!(c.rank(2), 1);

        let h = GenericHomology::from(&c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 1);
    }

    #[test]
    fn homology_t2() {
        let c = TestChainComplex::<i32>::t2();
        let c = Reduced::from(c);

        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);

        let h = GenericHomology::from(&c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[2].rank(), 1);
    }

    #[test]
    fn homology_rp2() {
        let c = TestChainComplex::<i32>::rp2();
        let c = Reduced::from(c);

        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 1);
        assert_eq!(c.rank(2), 1);

        let h = GenericHomology::from(&c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].tors(), &vec![2]);
        assert_eq!(h[2].rank(), 0);
    }
}