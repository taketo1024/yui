use std::mem;
use std::ops::Index;

use itertools::{Itertools, Either};
use delegate::delegate;
use yui::{Ring, RingOps};
use yui_matrix::sparse::{SpMat, SpVec, MatType, Trans};

use crate::{GridTrait, GridDeg, Grid, GridIter, isize2, isize3};
use crate::utils::ChainReducer;
use super::{RModStr, SimpleRModStr, rmod_str_symbol};

pub type ChainComplexSummand<R> = SimpleRModStr<R>;

pub type ChainComplex<R>  = ChainComplexBase<isize,  R>;
pub type ChainComplex2<R> = ChainComplexBase<isize2, R>;
pub type ChainComplex3<R> = ChainComplexBase<isize3, R>;

pub trait ChainComplexTrait<I>: Sized
where 
    I: GridDeg, 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> 
{ 
    type R;
    type Element;

    fn rank(&self, i: I) -> usize;
    fn d_deg(&self) -> I;
    fn d(&self, i: I, z: &Self::Element) -> Self::Element;
    fn d_matrix(&self, i: I) -> SpMat<Self::R>;
}

pub trait ChainComplexCommon<I>: GridTrait<I> + ChainComplexTrait<I>
where 
    I: GridDeg, 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>
{
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

    fn display_d_at(&self, i: I) -> String {
        let c = |i| rmod_str_symbol(self.rank(i), &[], "0");
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

    fn as_generic(&self) -> ChainComplexBase<I, Self::R> {
        ChainComplexBase::generate(self.support(), self.d_deg(), |i| self.d_matrix(i))
    }
}

impl<I, C> ChainComplexCommon<I> for C 
where 
    I: GridDeg,
    C: GridTrait<I> + ChainComplexTrait<I>,
    C::R: Ring, for<'x> &'x C::R: RingOps<C::R>
{}

pub struct ChainComplexBase<I, R>
where I: GridDeg, R: Ring, for<'x> &'x R: RingOps<R> {
    summands: Grid<I, ChainComplexSummand<R>>,
    d_deg: I,
    d_matrices: Grid<I, SpMat<R>>
}

impl<I, R> ChainComplexBase<I, R> 
where I: GridDeg, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(summands: Grid<I, ChainComplexSummand<R>>, d_deg: I, mut d_matrices: Grid<I, SpMat<R>>) -> Self { 
        assert!(summands.iter().all(|(_, s)| s.is_free()));
        
        for i in summands.support() { 
            let r = summands[i].rank();
            if r > 0 && !d_matrices.is_supported(i - d_deg) { 
                let d = SpMat::zero((r, 0));
                d_matrices.insert(i - d_deg, d);
            }
        }

        Self { summands, d_deg, d_matrices }
    }

    pub fn generate<It, F>(support: It, d_deg: I, d_matrix_map: F) -> Self
    where 
        It: IntoIterator<Item = I>, 
        F: FnMut(I) -> SpMat<R>
    {
        let d_matrices = Grid::generate(support, d_matrix_map);

        let summands = Grid::generate(d_matrices.support(), |i| {
            let r = d_matrices[i].cols();
            let t = Trans::id(r);
            ChainComplexSummand::new(r, vec![], Some(t)) 
        });

        Self::new(summands, d_deg, d_matrices)
    }

    pub fn summands(&self) -> &Grid<I, ChainComplexSummand<R>> { 
        &self.summands
    }

    pub fn d_matrix_ref(&self, i: I) -> &SpMat<R> {
        &self.d_matrices[i]
    }

    pub fn reduced(&self, with_trans: bool) -> ChainComplexBase<I, R> { 
        ChainReducer::reduce(self, with_trans)
    }
}

impl<I, R> GridTrait<I> for ChainComplexBase<I, R>
where I: GridDeg, R: Ring, for<'x> &'x R: RingOps<R> {
    type Itr = GridIter<I>;
    type Output = ChainComplexSummand<R>;

    delegate! { 
        to self.summands { 
            fn support(&self) -> Self::Itr;        
            fn is_supported(&self, i: I) -> bool;
            fn get(&self, i: I) -> &Self::Output;
        }
    }
}

impl<I, R> ChainComplexTrait<I> for ChainComplexBase<I, R>
where I: GridDeg, R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;
    type Element = SpVec<R>;

    fn rank(&self, i: I) -> usize {
        self[i].rank()
    }

    fn d_deg(&self) -> I { 
        self.d_deg
    }

    fn d_matrix(&self, i: I) -> SpMat<Self::R> {
        self.d_matrices[i].clone()
    }

    fn d(&self, i: I, v: &Self::Element) -> Self::Element {
        assert_eq!(self.rank(i), v.dim());
        let d = self.d_matrix_ref(i);
        d * v
    }
}

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
        
        Self::generate(
            range, d_deg, 
            move |i| {
                let i = (i - offset) as usize;
                mem::take(&mut mats[i])
            }
        )
    }
}

impl<I, R> Index<I> for ChainComplexBase<I, R>
where I: GridDeg, R: Ring, for<'x> &'x R: RingOps<R> {
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