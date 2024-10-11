use std::ops::{Index, RangeInclusive};
use std::sync::Arc;

use delegate::delegate;
use num_traits::Zero;
use yui::{Ring, RingOps};
use yui::lc::{Gen, Lc};
use yui_matrix::sparse::{SpMat, SpVec};

use crate::utils::ChainReducer;
use crate::{isize2, isize3, ChainComplexTrait, Grid, Grid1, GridDeg, GridIter, GridTrait, SummandTrait};
use super::Summand;

#[cfg(feature = "multithread")]
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

pub type ChainComplex <X, R> = ChainComplexBase<isize,  X, R>;
pub type ChainComplex2<X, R> = ChainComplexBase<isize2, X, R>;
pub type ChainComplex3<X, R> = ChainComplexBase<isize3, X, R>;

pub type XChainComplexSummand<X, R> = Summand<X, R>;

#[derive(Clone)]
pub struct ChainComplexBase<I, X, R>
where 
    I: GridDeg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    summands: Grid<I, XChainComplexSummand<X, R>>,
    d_deg: I,
    d_map: Arc<dyn Fn(I, &Lc<X, R>) -> Lc<X, R> + Send + Sync>,
}

impl<I, X, R> ChainComplexBase<I, X, R>
where 
    I: GridDeg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    pub fn new<F>(summands: Grid<I, XChainComplexSummand<X, R>>, d_deg: I, d_map: F) -> Self
    where F: Fn(I, &Lc<X, R>) -> Lc<X, R> + Send + Sync + 'static {
        assert!(summands.iter().all(|(_, s)| s.is_free()));

        let d_map = Arc::new(d_map);
        Self { summands, d_deg, d_map }
    }

    pub fn zero() -> Self { 
        Self::new(Grid::default(), I::zero(), |_, _| Lc::zero())
    }

    pub fn summands(&self) -> &Grid<I, XChainComplexSummand<X, R>> { 
        &self.summands
    }

    fn d_matrix(&self, i: I) -> SpMat<R> { 
        let m = self[i + self.d_deg].rank();
        let n = self[i].rank();

        cfg_if::cfg_if! { 
            if #[cfg(feature = "multithread")] {
                let cols = (0..n).into_par_iter().map(|j|
                    self.d_matrix_col(i, j)
                ).collect::<Vec<_>>();
                SpMat::from_col_vecs(m, cols)
            } else { 
                let cols = (0..n).map(|j|
                    self.d_matrix_col(i, j)
                );
                SpMat::from_col_vecs(m, cols)
            }
        }
    }

    #[inline(never)] // for profilability
    fn d_matrix_col(&self, i: I, j: usize) -> SpVec<R> { 
        let z = self[i].gen(j);
        let w = self.d(i, &z);
        self[i + self.d_deg].vectorize(&w)
    }

    pub fn reduced(&self) -> ChainComplexBase<I, X, R> { 
        let r = ChainReducer::reduce(self, true);

        let summands = Grid::generate(
            self.summands.support(),
            |i| {
                let c = &self[i];
                Summand::new(
                    c.raw_gens().clone(), 
                    r.rank(i).unwrap(), 
                    vec![], 
                    c.trans().merged(r.trans(i).unwrap())
                )
            } 
        );

        let d_deg = self.d_deg;
        let d_map = self.d_map.clone();
        Self { summands, d_deg, d_map }
    }
}

impl<X, R> ChainComplex<X, R>
where 
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    pub fn truncated(&self, range: RangeInclusive<isize>) -> Self { 
        let d_deg = self.d_deg;
        let d_map = self.d_map.clone();
        let summands = self.summands.truncated(range.clone());

        Self::new(summands, d_deg, move |i, z| 
            if range.contains(&(i + d_deg)) { 
                d_map(i, z)
            } else { 
                Lc::zero()
            }
        )
    }
}

impl<I, X, R> GridTrait<I> for ChainComplexBase<I, X, R>
where 
    I: GridDeg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    type Itr = GridIter<I>;
    type Output = XChainComplexSummand<X, R>;
    
    delegate! { 
        to self.summands { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: I) -> bool;
            fn get(&self, i: I) -> &Self::Output;
        }
    }
}

impl<I, X, R> ChainComplexTrait<I> for ChainComplexBase<I, X, R>
where 
    I: GridDeg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    type R = R;
    type Element = Lc<X, R>;

    fn rank(&self, i: I) -> usize {
        self[i].rank()
    }
    
    fn d_deg(&self) -> I {
        self.d_deg
    }

    fn d(&self, i: I, z: &Lc<X, R>) -> Lc<X, R> { 
        (self.d_map)(i, z)
    }

    fn d_matrix(&self, i: I) -> SpMat<Self::R> { 
        self.d_matrix(i)
    }
}

impl<I, X, R> Index<I> for ChainComplexBase<I, X, R>
where I: GridDeg, X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = XChainComplexSummand<X, R>;
    fn index(&self, i: I) -> &Self::Output {
        self.get(i)
    }
}

impl<X, R> Index<(isize, isize)> for ChainComplex2<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = XChainComplexSummand<X, R>;
    fn index(&self, i: (isize, isize)) -> &Self::Output {
        self.get(i.into())
    }
}

impl<X, R> Index<(isize, isize, isize)> for ChainComplex3<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = XChainComplexSummand<X, R>;
    fn index(&self, i: (isize, isize, isize)) -> &Self::Output {
        self.get(i.into())
    }
}