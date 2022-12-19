use std::collections::HashMap;

use sprs::{CsMat, PermView};

use crate::math::matrix::pivot::{find_pivots, perms_by_pivots};
use crate::math::matrix::schur::schur_partial_upper_triang;
use crate::math::matrix::sparse::CsMatExt;
use crate::math::traits::{Ring, RingOps};
use crate::math::matrix::CsMatElem;
use super::complex::ChainComplex;

pub struct Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring + CsMatElem, 
    for<'x> &'x C::R: RingOps<C::R>  
{ 
    d_matrices: HashMap<isize, CsMat<C::R>>,
    original: C,
}

impl<C> Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring + CsMatElem, 
    for<'x> &'x C::R: RingOps<C::R>  
{ 
    pub fn new(c: C) -> Self {
        let mut d_matrices = HashMap::new();
        let range: Box<dyn Iterator<Item = _>> = if c.d_degree() > 0 { 
            Box::new(c.range())
        } else {
            Box::new(c.range().rev())
        };

        for k in range {
            let deg = c.d_degree();
            let (k0, k1, k2) = (k - deg, k, k + deg);

            let a0 = d_matrices.remove(&k0); // prev(optional)
            let a1 = d_matrices.remove(&k1).unwrap_or(c.d_matrix(k1)); // target
            let a2 = c.d_matrix(k2); // next

            let (b0, b1, b2) = Self::reduce(k, a0, a1, a2);
            
            if let Some(b0) = b0 {
                d_matrices.insert(k0, b0);
            }
            d_matrices.insert(k1, b1);
            d_matrices.insert(k2, b2);
        }

        Self { d_matrices, original: c }
    }

    fn reduce(k: isize, a0: Option<CsMat<C::R>>, a1: CsMat<C::R>, a2: CsMat<C::R>) -> (Option<CsMat<C::R>>, CsMat<C::R>, CsMat<C::R>) {
        if let Some(a0) = a0.as_ref() {
            assert_eq!(a0.rows(), a1.cols());
        }
        assert_eq!(a1.rows(), a2.cols());

        let pivs = find_pivots(&a1);
        let r = pivs.len();

        println!("k = {k}, reduce: {r} / {}", a1.cols());
    
        if r == 0 { 
            return (a0, a1, a2)
        }
    
        let (p, q) = perms_by_pivots(&a1, &pivs);
        let b1 = a1.permute(p.view(), q.view());
        let b1 = schur_partial_upper_triang(b1, r);
        
        let b0 = a0.map(|a0| {
            Self::reduce_rows(a0, q.view(), r)
        });
        let b2 = Self::reduce_cols(a2, p.view(), r);

        (b0, b1, b2)
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

impl<C> ChainComplex for Reduced<C>
where 
    C: ChainComplex,
    C::R: Ring + CsMatElem, 
    for<'x> &'x C::R: RingOps<C::R>  
{
    type R = C::R;
    type Generator = C::Generator;

    fn range(&self) -> std::ops::RangeInclusive<isize> {
        self.original.range()
    }

    fn rank(&self, k: isize) -> usize {
        if self.range().contains(&k) {
            self.d_matrices[&k].cols()
        } else {
            0
        }
    }

    fn d_degree(&self) -> isize {
        self.original.d_degree()
    }

    fn d_matrix(&self, k: isize) -> sprs::CsMat<Self::R> {
        if self.range().contains(&k) { 
            self.d_matrices[&k].clone()
        } else {
            let m = self.rank(k + self.d_degree());
            let n = self.rank(k);
            CsMat::zero((m, n))
        }
    }

    fn generators(&self, _k: isize) -> Vec<&Self::Generator> {
        todo!()
    }

    fn differentiate(&self, _k: isize, _x:&Self::Generator) -> Vec<(Self::Generator, Self::R)> {
        todo!()
    }
}