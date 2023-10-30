use std::collections::HashMap;
use std::mem;
use std::ops::Index;

use itertools::{Itertools, Either};
use delegate::delegate;
use yui_core::{Ring, RingOps, EucRing, EucRingOps, Deg, isize2, isize3};
use yui_matrix::sparse::{SpMat, SpVec, MatType, Trans};

use crate::utils::r_mod_str;
use crate::{ReducedComplexBase, GridBase, GridIter, DisplayForGrid};

use super::grid::GridTrait;
use super::utils::ChainReducer;
use super::homology::HomologyBase;

pub type ChainComplex<R>  = ChainComplexBase<isize,  R>;
pub type ChainComplex2<R> = ChainComplexBase<isize2, R>;
pub type ChainComplex3<R> = ChainComplexBase<isize3, R>;

pub trait ChainComplexSummandTrait: DisplayForGrid
where Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> { 
    type R;
    fn rank(&self) -> usize;
    fn d_matrix(&self) -> &SpMat<Self::R>;

    fn d(&self, v: &SpVec<Self::R>) -> SpVec<Self::R> {
        assert_eq!(self.rank(), v.dim());
        let d = self.d_matrix();
        d * v
    }

    fn is_cycle(&self, v: &SpVec<Self::R>) -> bool { 
        self.d(v).is_zero()
    }

    fn module_str(&self) -> String { 
        r_mod_str(self.rank(), vec![].iter())
    }
}

pub trait ChainComplexTrait<I>: GridTrait<I> + Sized
where 
    I: Deg, 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> 
{ 
    type R;

    fn rank(&self, i: I) -> usize;
    fn d_deg(&self) -> I;
    fn d_matrix(&self, i: I) -> &SpMat<Self::R>;

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
    
    fn reduced(&self, with_trans: bool) -> ReducedComplexBase<I, Self::R> { 
        ChainReducer::reduce(self, with_trans)
    }

    fn reduced_by<F>(&self, trans: F) -> ReducedComplexBase<I, Self::R>
    where F: FnMut(I) -> Trans<Self::R> {
        ReducedComplexBase::reduced_by(self, trans)
    }
}

pub trait ChainComplexDisplay<I>: ChainComplexTrait<I>
where 
    I: Deg, 
    Self::E: DisplayForGrid,
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> 
{
    fn display_d_at(&self, i: I) -> String {
        let c = |i| self.get(i).display_for_grid();
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
}

#[derive(Default)]
pub struct ChainComplexSummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    d_matrix: SpMat<R>
}

impl<R> ChainComplexSummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(d_matrix: SpMat<R>) -> Self { 
        Self { d_matrix }
    }
}

impl<R> DisplayForGrid for ChainComplexSummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn display_for_grid(&self) -> String {
        self.module_str()
    }
}

impl<R> ChainComplexSummandTrait for ChainComplexSummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn rank(&self) -> usize {
        self.d_matrix.cols()
    }

    fn d_matrix(&self) -> &SpMat<Self::R> {
        &self.d_matrix
    }
}

pub struct ChainComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    d_deg: I,
    summands: GridBase<I, ChainComplexSummand<R>>
}

impl<I, R> ChainComplexBase<I, R> 
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new<It, F>(support: It, d_deg: I, mut d_matrix: F) -> Self
    where 
        It: Iterator<Item = I>, 
        F: FnMut(I) -> SpMat<R>
    {
        let support = support.collect_vec();
        let mut data = HashMap::<I, SpMat<R>>::from_iter( 
            support.iter().map(|&i| (i, d_matrix(i))) 
        );

        // add zero maps
        for &i in support.iter() {
            let r = data[&i].shape().1;
            if r > 0 && !data.contains_key(&(i - d_deg)) { 
                let d = SpMat::zero((r, 0));
                data.insert(i - d_deg, d);
            }
        }

        let summands = GridBase::new_raw(
            support, 
            data.into_iter().map(|(i, d)| 
                (i, ChainComplexSummand::new(d))
            ).collect()
        );

        Self { d_deg, summands }
    }
}

impl<I, R> Index<I> for ChainComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = ChainComplexSummand<R>;
    fn index(&self, i: I) -> &Self::Output {
        self.get(i)
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

    fn rank(&self, i: I) -> usize {
        self[i].rank()
    }

    fn d_matrix(&self, i: I) -> &SpMat<Self::R> {
        self[i].d_matrix()
    }
}

impl<I, R> ChainComplexDisplay<I> for ChainComplexBase<I, R>
where I: Deg, R: Ring, for<'x> &'x R: RingOps<R> {}

impl<I, R> ChainComplexBase<I, R> 
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology(&self, with_trans: bool) -> HomologyBase<I, R> {
        HomologyBase::compute_from(self, with_trans)
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