use std::ops::{RangeInclusive, Index};
use itertools::Itertools;

use crate::math::traits::{Ring, RingOps};
use super::homology::{HomologySummand, HomologyComputable, Homology};

pub struct SimpleHomology<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    summands: Vec<HomologySummand<R>>,
    zero: HomologySummand<R>,
    shift: isize
}

impl<R> SimpleHomology<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(summands: Vec<HomologySummand<R>>, shift: isize) -> Self { 
        let zero = HomologySummand::zero();
        SimpleHomology { summands, zero, shift }
    }
}

impl<C> From<C> for SimpleHomology<C::R>
where
    C: HomologyComputable,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>  
{
    fn from(c: C) -> Self {
        let summands = c.range().map(|k| c.homology_at(k)).collect_vec();
        let shift = *c.range().start();
        Self::new(summands, shift)
    }
}

impl<R> Index<isize> for SimpleHomology<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = HomologySummand<R>;

    fn index(&self, k: isize) -> &Self::Output {
        if self.range().contains(&k) {
            let index = (k - self.shift) as usize;
            &self.summands[index]
        } else {
            &self.zero
        }
    }
}

impl<R> Homology for SimpleHomology<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;
    fn range(&self) -> RangeInclusive<isize> {
        let s = self.shift;
        let n = self.summands.len() as isize;
        s ..= s + n - 1
    }
}

#[cfg(test)]
mod tests { 
    use sprs::CsMat;
    use crate::math::matrix::sparse::CsMatExt;

    use super::*;
    use super::super::simple_complex::*;
    use super::super::simple_complex::tests::*;

    #[test]
    fn cancel_pair() { 
        let c = SimpleChainComplex::new(
            vec![ CsMat::csc_from_vec((1, 1), vec![1]) ],
            -1
        );

        let h = SimpleHomology::from(c);
        
        assert_eq!(h[0].rank(), 0);
        assert!(h[0].is_free());
    }

    #[test]
    fn torsion() { 
        let c = SimpleChainComplex::new( 
            vec![ CsMat::csc_from_vec((1, 1), vec![2]) ],
            -1
        );

        let h = SimpleHomology::from(c);
        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[0].tors(), &vec![2]);
    }

    #[test]
    fn homology_d3() {
        let c = sample_d3();
        let h = SimpleHomology::from(c);

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
        let c = sample_s2();
        let h = SimpleHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn homology_t2() {
        let c = sample_t2();
        let h = SimpleHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn homology_rp2() {
        let c = sample_rp2();
        let h = SimpleHomology::from(c);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].tors(), &vec![2]);
        assert_eq!(h[1].is_free(), false);

        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[2].is_free(), true);
    }
}