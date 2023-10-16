use std::cell::OnceCell;
use std::collections::HashMap;

use itertools::Itertools;
use yui_core::{Ring, RingOps, EucRing, EucRingOps};
use yui_matrix::sparse::SpMat;

use super::deg::{Deg, isize2, isize3};
use super::graded::Graded;
use super::homology::HomologyBase;

pub type ChainComplex<R>  = ChainComplexBase<isize,  R>;
pub type ChainComplex2<R> = ChainComplexBase<isize2, R>;
pub type ChainComplex3<R> = ChainComplexBase<isize3, R>;

pub struct ChainComplexBase<D, R>
where D: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    support: Vec<D>,
    d_deg: D,
    f_rank: Box<dyn Fn(D) -> usize>,
    f_d:    Box<dyn Fn(D) -> SpMat<R>>,
    cache_rank: HashMap<D, OnceCell<usize>>,
    cache_d   : HashMap<D, OnceCell<SpMat<R>>>,
    zero_d: SpMat<R>
}

impl<I, R> ChainComplexBase<I, R> 
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new<It, F1, F2>(support: It, d_deg: I, f_rank: F1, f_d: F2) -> Self
    where 
        It: Iterator<Item = I>, 
        F1: Fn(I) -> usize + 'static,
        F2: Fn(I) -> SpMat<R> + 'static 
    {
        let support = support.collect_vec();
        let f_rank = Box::new( f_rank );
        let f_d = Box::new( f_d );
        let cache_rank = HashMap::from_iter( support.iter().map(|&i| (i, OnceCell::new())) );
        let cache_d    = HashMap::from_iter( support.iter().map(|&i| (i, OnceCell::new())) );
        let zero_d = SpMat::zero((0, 0));

        Self { 
            support, d_deg, f_rank, f_d, cache_rank, cache_d, zero_d
        }
    }

    pub fn rank(&self, i: I) -> usize { 
        if self.is_supported(i) { 
            self.cache_rank[&i].get_or_init(|| { 
                (self.f_rank)(i)
            }).clone()
        } else { 
            0
        }
    }

    pub fn d_deg(&self) -> I { 
        self.d_deg
    }

    pub fn d_matrix(&self, i: I) -> &SpMat<R> { 
        if self.is_supported(i) { 
            self.cache_d[&i].get_or_init(|| { 
                (self.f_d)(i)
            })
        } else { 
            &self.zero_d
        }
    }

    pub fn is_supported(&self, i: I) -> bool { 
        self.support.contains(&i)
    }

    pub fn check_d_at(&self, i: I) { 
        let i1 = i + self.d_deg;
        if !(self.is_supported(i) && self.is_supported(i1)) {
            return 
        }

        let d0 = self.d_matrix(i);
        let d1 = self.d_matrix(i1);
        let res = d1 * d0;

        assert!( res.is_zero(), "d² is non-zero at {i}." );
    }

    pub fn check_d_all(&self) {
        for &i in self.support.iter() { 
            self.check_d_at(i);
        }
    }
}

impl<I, R> ChainComplexBase<I, R> 
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology(self) -> HomologyBase<I, R> {
        HomologyBase::new(self)
    }
}

impl<D, R> Graded<D> for ChainComplexBase<D, R>
where D: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type Itr = std::vec::IntoIter<D>;
    fn support(&self) -> Self::Itr {
        self.support.clone().into_iter()
    }

    fn display(&self, i: D) -> String {
        use yui_utils::superscript;

        let symbol = R::set_symbol();
        let rank = self.rank(i);
        if rank > 1 {
            format!("{}{}", symbol, superscript(rank as isize))
        } else if rank == 1 { 
            symbol
        } else { 
            String::from("0")
        }
    }
}

#[cfg(test)]
pub(crate) mod tests { 
    use super::*;

    impl<R> ChainComplexBase<isize, R>
    where R: Ring, for<'x> &'x R: RingOps<R> {
        pub fn from_mats(d_deg: isize, mats: Vec<SpMat<R>>) -> Self { 
            use yui_matrix::sparse::MatType;
    
            let n = mats.len() as isize;
            let ranks = mats.iter().map(|d| d.shape().1).collect_vec();
            Self::new(
                0..n, d_deg, 
                move |i| if 0 <= i && i < n { ranks[i as usize] } else { 0 }, 
                move |i| if 0 <= i && i < n { mats[i as usize].clone() } else { SpMat::zero((0, 0)) }
            )
        }
    }
    
    pub(crate) struct Samples<R>
    where R: Ring, for<'x> &'x R: RingOps<R> {
        _r: R
    }

    impl<R> Samples<R> 
    where R: Ring, for<'x> &'x R: RingOps<R> {
        fn mat(shape: (usize, usize), entries: Vec<i32>) -> SpMat<R> { 
            SpMat::from_vec(shape, entries.into_iter().map(|x| R::from(x)).collect())
        }

        pub fn d3() -> ChainComplex<R> {
            ChainComplex::from_mats(-1,
                vec![
                    Self::mat((0, 4), vec![]),
                    Self::mat((4, 6), vec![-1, -1, 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 1] ),
                    Self::mat((6, 4), vec![1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1] ),
                    Self::mat((4, 1), vec![-1, 1, -1, 1]),
                    Self::mat((1, 0), vec![]),
                ]
            )
        }

        pub fn s2() -> ChainComplex<R> {
            ChainComplex::from_mats(-1,
                vec![
                    Self::mat((0, 4), vec![]),
                    Self::mat((4, 6), vec![-1, -1, 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 1]),
                    Self::mat((6, 4), vec![1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1] ),
                    Self::mat((4, 0), vec![])
                ]
            )
        }
    
        pub fn t2() -> ChainComplex<R> {
            ChainComplex::from_mats(-1,
                vec![
                    Self::mat((0, 9), vec![]),
                    Self::mat((9, 27), vec![-1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, -1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1]),
                    Self::mat((27, 18), vec![1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1]),
                    Self::mat((18, 0), vec![])
                ]
            )
        }
    
        pub fn rp2() -> ChainComplex<R> { 
            ChainComplex::from_mats(-1,
                vec![
                    Self::mat((0, 6), vec![]),
                    Self::mat((6, 15), vec![-1, -1, 0, 0, 0, 0, 0, -1, -1, 0, -1, 0, 0, 0, 0, 1, 0, -1, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 1, 1, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1] ),
                    Self::mat((15, 10), vec![1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1] ),
                    Self::mat((10, 0), vec![]),
                ]
            )
        }
    }

    #[test]
    fn d3() { 
        let c = Samples::<i64>::d3();

        assert_eq!(c.rank(0), 4);
        assert_eq!(c.rank(1), 6);
        assert_eq!(c.rank(2), 4);
        assert_eq!(c.rank(3), 1);

        c.check_d_all();
    }

    #[test]
    fn s2() { 
        let c = Samples::<i64>::s2();

        assert_eq!(c.rank(0), 4);
        assert_eq!(c.rank(1), 6);
        assert_eq!(c.rank(2), 4);
        assert_eq!(c.rank(3), 0);

        c.check_d_all();
    }

    #[test]
    fn t2() { 
        let c = Samples::<i64>::t2();

        assert_eq!(c.rank(0), 9);
        assert_eq!(c.rank(1), 27);
        assert_eq!(c.rank(2), 18);
        assert_eq!(c.rank(3), 0);

        c.check_d_all();
    }

    #[test]
    fn rp2() { 
        let c = Samples::<i64>::rp2();
        
        assert_eq!(c.rank(0), 6);
        assert_eq!(c.rank(1), 15);
        assert_eq!(c.rank(2), 10);
        assert_eq!(c.rank(3), 0);

        c.check_d_all();
    }
}