use std::collections::HashMap;
use std::mem;

use itertools::Itertools;
use yui_core::{Ring, RingOps, EucRing, EucRingOps};
use yui_matrix::sparse::{SpMat, SpVec, MatType};

use super::deg::{Deg, isize2, isize3};
use super::graded::Graded;
use super::homology::HomologyBase;

pub type ChainComplex<R>  = ChainComplexBase<isize,  R>;
pub type ChainComplex2<R> = ChainComplexBase<isize2, R>;
pub type ChainComplex3<R> = ChainComplexBase<isize3, R>;

pub trait ChainComplexTrait<I>: Graded<I>
where I: Deg, Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> { 
    type R;

    fn d_deg(&self) -> I;
    fn d_matrix(&self, i: I) -> &SpMat<Self::R>;

    fn rank(&self, i: I) -> usize { 
        if self.is_supported(i) { 
            self.d_matrix(i).shape().1 // column
        } else { 
            0
        }
    }

    fn check_d_at(&self, i: I) { 
        let i1 = i + self.d_deg();
        if !(self.is_supported(i) && self.is_supported(i1)) {
            return 
        }

        let d0 = self.d_matrix(i);
        let d1 = self.d_matrix(i1);
        let res = d1 * d0;

        assert!( res.is_zero(), "d² is non-zero at {i}." );
    }

    fn check_d_all(&self) {
        for i in self.support() { 
            self.check_d_at(i);
        }
    }
    
    fn display_d_at(&self, i: I) -> String {
        let c = |i| self.display(i);
        let c0 = c(i);
        let c1 = c(i + self.d_deg());
        let d = self.d_matrix(i).to_dense();
        format!("C[{i}]: {c0} -> {c1}\n{d}") 
    }

    fn display_d(&self) -> String { 
        self.support().filter_map(|i| 
            if self.rank(i) > 0 && self.rank(i + self.d_deg()) > 0 {
                Some(self.display_d_at(i))
            } else { 
                None
            }
        ).join("\n\n")
    }

    fn print_d(&self) {
        println!("{}", self.display_d());
    }
}

pub struct ChainComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    support: Vec<I>,
    d_deg: I,
    d_matrix: HashMap<I, SpMat<R>>,
    zero_d: SpMat<R>
}

impl<I, R> ChainComplexBase<I, R> 
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new<It, F>(support: It, d_deg: I, mut d_matrix: F) -> Self
    where 
        It: Iterator<Item = I>, 
        F: FnMut(I) -> SpMat<R>
    {
        let support = support.collect_vec();
        let mut d_matrix = HashMap::from_iter( 
            support.iter().map(|&i| (i, d_matrix(i))) 
        );

        // add zero maps
        for &i in support.iter() {
            if !d_matrix.contains_key(&(i - d_deg)) { 
                let r = d_matrix[&i].shape().1;
                let d = SpMat::zero((r, 0));
                d_matrix.insert(i - d_deg, d);
            }
        }
        let zero_d = SpMat::zero((0, 0));

        Self { 
            support, d_deg, d_matrix, zero_d
        }
    }

    pub fn differentiate(&self, i: I, v: &SpVec<R>) -> SpVec<R> {
        assert_eq!(self.rank(i), v.dim());
        let d = self.d_matrix(i);
        d * v
    }

    pub fn is_cycle(&self, i: I, v: &SpVec<R>) -> bool { 
        self.differentiate(i, v).is_zero()
    }
}

impl<I, R> ChainComplexBase<I, R> 
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology(self) -> HomologyBase<I, R> {
        HomologyBase::new(self, false)
    }

    pub fn homology_with_gens(self) -> HomologyBase<I, R> {
        HomologyBase::new(self, true)
    }
}

impl<I, R> Graded<I> for ChainComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type Itr = std::vec::IntoIter<I>;

    fn support(&self) -> Self::Itr {
        self.support.clone().into_iter()
    }

    fn display(&self, i: I) -> String {
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

impl<I, R> ChainComplexTrait<I> for ChainComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn d_deg(&self) -> I { 
        self.d_deg
    }

    fn d_matrix(&self, i: I) -> &SpMat<R> { 
        if let Some(d) = self.d_matrix.get(&i) {
            d
        } else { 
            &self.zero_d
        }
    }
}


impl<R> ChainComplexBase<isize, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_mats(d_deg: isize, offset: isize, mut mats: Vec<SpMat<R>>) -> Self { 
        let n = mats.len() as isize;
        let range = offset .. offset + n;
        Self::new(
            range, d_deg, 
            move |i| {
                let i = (i - offset) as usize;
                mem::take(&mut mats[i])
            }
        )
    }
}

#[cfg(test)]
pub(crate) mod tests { 
    use super::*;

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
            ChainComplex::from_mats(-1, 0,
                vec![
                    Self::mat((0, 4), vec![]),
                    Self::mat((4, 6), vec![-1, -1, 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 1] ),
                    Self::mat((6, 4), vec![1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1] ),
                    Self::mat((4, 1), vec![-1, 1, -1, 1]),
                ]
            )
        }

        pub fn s2() -> ChainComplex<R> {
            ChainComplex::from_mats(-1, 0,
                vec![
                    Self::mat((0, 4), vec![]),
                    Self::mat((4, 6), vec![-1, -1, 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 1]),
                    Self::mat((6, 4), vec![1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1] ),
                ]
            )
        }
    
        pub fn t2() -> ChainComplex<R> {
            ChainComplex::from_mats(-1, 0,
                vec![
                    Self::mat((0, 9), vec![]),
                    Self::mat((9, 27), vec![-1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, -1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1]),
                    Self::mat((27, 18), vec![1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1]),
                ]
            )
        }
    
        pub fn rp2() -> ChainComplex<R> { 
            ChainComplex::from_mats(-1, 0,
                vec![
                    Self::mat((0, 6), vec![]),
                    Self::mat((6, 15), vec![-1, -1, 0, 0, 0, 0, 0, -1, -1, 0, -1, 0, 0, 0, 0, 1, 0, -1, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 1, 1, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1] ),
                    Self::mat((15, 10), vec![1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1] ),
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