use std::mem;
use std::ops::Index;

use itertools::{Itertools, Either};
use delegate::delegate;
use yui_core::{Ring, RingOps, Deg, isize2, isize3};
use yui_matrix::sparse::{SpMat, SpVec, MatType, Trans};

use crate::{GridBase, GridIter, RModStr, SimpleRModStr};

use super::grid::GridTrait;
use super::utils::ChainReducer;

pub type ChainComplexSummand<R> = SimpleRModStr<R>;

pub type ChainComplex<R>  = ChainComplexBase<isize,  R>;
pub type ChainComplex2<R> = ChainComplexBase<isize2, R>;
pub type ChainComplex3<R> = ChainComplexBase<isize3, R>;

pub trait ChainComplexTrait<I>: GridTrait<I> + Sized
where 
    I: Deg, 
    Self::E: RModStr<R = Self::R>,
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> 
{ 
    type R;

    fn d_deg(&self) -> I;
    fn d_matrix(&self, i: I) -> SpMat<Self::R>;

    fn check_d_at(&self, i0: I) { 
        let i1 = i0 + self.d_deg();
        if !(self.is_supported(i0) && self.is_supported(i1)) {
            return 
        }

        let d0 = self.d_matrix(i0);
        let d1 = self.d_matrix(i1);
        let res = d1 * d0;

        assert!( res.is_zero(), "d² is non-zero at {i0}." );
    }

    fn check_d_all(&self) {
        for i in self.support() { 
            self.check_d_at(i);
        }
    }

    fn as_generic(&self) -> ChainComplexBase<I, Self::R> {
        ChainComplexBase::new(self.support(), self.d_deg(), |i| self.d_matrix(i))
    }
}

pub trait ChainComplexDisplay<I>: ChainComplexTrait<I>
where 
    I: Deg, 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::E: RModStr<R = Self::R>,
{
    fn display_d_at(&self, i: I) -> String {
        let c = |i| self.get(i).math_symbol();
        let c0 = c(i);
        let c1 = c(i + self.d_deg());
        let d = self.d_matrix(i).to_dense();

        if d.is_zero() { 
            format!("C[{i}] {c0} -> {c1}; zero.")
        } else { 
            format!("C[{i}] {c0} -> {c1}\n{d}.")
        }
    }

    fn display_d(&self) -> String { 
        self.support().filter_map(|i| 
            if self.get(i).rank() > 0 && self.get(i + self.d_deg()).rank() > 0 {
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
    summands: GridBase<I, ChainComplexSummand<R>>,
    d_deg: I,
    d_matrices: GridBase<I, SpMat<R>>
}

impl<I, R> ChainComplexBase<I, R> 
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new<It, F>(support: It, d_deg: I, mut d_matrix: F) -> Self
    where 
        It: Iterator<Item = I>, 
        F: FnMut(I) -> SpMat<R>
    {
        Self::new_with_trans(support, d_deg, move |i| {
            let d = d_matrix(i);
            let r = d.cols();
            let t = Trans::id(r);
            (d, Some(t))
        })
    }

    pub fn new_with_trans<It, F>(support: It, d_deg: I, mut d_matrix: F) -> Self
    where 
        It: Iterator<Item = I>, 
        F: FnMut(I) -> (SpMat<R>, Option<Trans<R>>)
    {
        use std::mem::take;

        let support = support.collect_vec().into_iter();
        let mut data = GridBase::new(support.clone(), |i| d_matrix(i));

        let mut d_matrices = GridBase::new(support.clone(), |i| 
            take(&mut data.get_mut(i).unwrap().0)
        );

        let summands = GridBase::new(support.clone(), |i| {
            let r = d_matrices[i].cols();
            let t = take(&mut data.get_mut(i).unwrap().1);
            ChainComplexSummand::new(r, vec![], t) 
        });

        for i in support.clone() { 
            let r = summands[i].rank();
            if r > 0 && !d_matrices.is_supported(i - d_deg) { 
                let d = SpMat::zero((r, 0));
                d_matrices.insert(i - d_deg, d);
            }
        }

        Self { summands, d_deg, d_matrices }
    }

    pub fn d_matrix_ref(&self, i: I) -> &SpMat<R> {
        &self.d_matrices[i]
    }

    pub fn d(&self, i: I, v: &SpVec<R>) -> SpVec<R> {
        assert_eq!(self.get(i).rank(), v.dim());
        let d = self.d_matrix_ref(i);
        d * v
    }
    
    pub fn reduced(&self, with_trans: bool) -> ChainComplexBase<I, R> { 
        ChainReducer::reduce(self, with_trans)
    }
}

impl<I, R> GridTrait<I> for ChainComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type Itr = GridIter<I>;
    type E = ChainComplexSummand<R>;

    delegate! { 
        to self.summands { 
            fn support(&self) -> Self::Itr;        
            fn is_supported(&self, i: I) -> bool;
            fn get(&self, i: I) -> &Self::E;
        }
    }
}

impl<I, R> ChainComplexTrait<I> for ChainComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn d_deg(&self) -> I { 
        self.d_deg
    }

    fn d_matrix(&self, i: I) -> SpMat<Self::R> {
        self.d_matrices[i].clone()
    }
}

impl<I, R> ChainComplexDisplay<I> for ChainComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {}

impl<R> ChainComplexBase<isize, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_mats(d_deg: isize, offset: isize, mut mats: Vec<SpMat<R>>) -> Self { 
        let n = mats.len() as isize;
        let range = offset .. offset + n;
        let range = if d_deg.is_positive() { 
            Either::Left(range) 
        } else { 
            Either::Right(range.rev())
        };
        Self::new(
            range, d_deg, 
            move |i| {
                let i = (i - offset) as usize;
                mem::take(&mut mats[i])
            }
        )
    }
}

impl<I, R> Index<I> for ChainComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = ChainComplexSummand<R>;
    fn index(&self, i: I) -> &Self::Output {
        self.get(i)
    }
}

impl<R> Index<(isize, isize)> for ChainComplex2<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = ChainComplexSummand<R>;
    fn index(&self, i: (isize, isize)) -> &Self::Output {
        self.get(i.into())
    }
}

impl<R> Index<(isize, isize, isize)> for ChainComplex3<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = ChainComplexSummand<R>;
    fn index(&self, i: (isize, isize, isize)) -> &Self::Output {
        self.get(i.into())
    }
}

#[cfg(test)]
pub(crate) mod tests { 
    use super::*;

    #[test]
    fn d3() { 
        let c = ChainComplex::<i64>::d3();

        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 4);
        assert_eq!(c[3].rank(), 1);

        c.check_d_all();
    }

    #[test]
    fn s2() { 
        let c = ChainComplex::<i64>::s2();

        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 4);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();
    }

    #[test]
    fn t2() { 
        let c = ChainComplex::<i64>::t2();

        assert_eq!(c[0].rank(), 9);
        assert_eq!(c[1].rank(), 27);
        assert_eq!(c[2].rank(), 18);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();
    }

    #[test]
    fn rp2() { 
        let c = ChainComplex::<i64>::rp2();
        
        assert_eq!(c[0].rank(), 6);
        assert_eq!(c[1].rank(), 15);
        assert_eq!(c[2].rank(), 10);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();
    }
}