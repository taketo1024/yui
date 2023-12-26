use std::ops::Index;
use std::sync::Arc;

use delegate::delegate;
use yui::{Ring, RingOps};
use yui::lc::{Gen, Lc};
use yui_matrix::sparse::pivot::PivotCondition;
use yui_matrix::sparse::{SpMat, SpVec};

use crate::utils::ChainReducer;
use crate::{GridTrait, GridDeg, Grid, GridIter, ChainComplexTrait, RModStr, isize2, isize3};
use super::XModStr;

#[cfg(feature = "multithread")]
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

pub type XChainComplex <X, R> = XChainComplexBase<isize,  X, R>;
pub type XChainComplex2<X, R> = XChainComplexBase<isize2, X, R>;
pub type XChainComplex3<X, R> = XChainComplexBase<isize3, X, R>;

pub type XChainComplexSummand<X, R> = XModStr<X, R>;

pub struct XChainComplexBase<I, X, R>
where 
    I: GridDeg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    summands: Grid<I, XChainComplexSummand<X, R>>,
    d_deg: I,
    d_map: Arc<dyn Fn(I, &Lc<X, R>) -> Lc<X, R> + Send + Sync>,
}

impl<I, X, R> XChainComplexBase<I, X, R>
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
        let z = self[i].gen_chain(j);
        let w = self.d(i, &z);
        self[i + self.d_deg].vectorize(&w)
    }

    pub fn reduced(&self) -> XChainComplexBase<I, X, R> { 
        let reduced = ChainReducer::reduce(self, PivotCondition::AnyUnit, true);
        let summands = Grid::generate(
            self.summands.support(),
            |i| { 
                let mut s = self.summands[i].clone();
                s.merge(reduced[i].clone(), false);
                s
            }
        );

        let d_deg = self.d_deg;
        let d_map = self.d_map.clone();
        Self { summands, d_deg, d_map }
    }
}

impl<I, X, R> GridTrait<I> for XChainComplexBase<I, X, R>
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

impl<I, X, R> ChainComplexTrait<I> for XChainComplexBase<I, X, R>
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

impl<I, X, R> Index<I> for XChainComplexBase<I, X, R>
where I: GridDeg, X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = XChainComplexSummand<X, R>;
    fn index(&self, i: I) -> &Self::Output {
        self.get(i)
    }
}

impl<X, R> Index<(isize, isize)> for XChainComplex2<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = XChainComplexSummand<X, R>;
    fn index(&self, i: (isize, isize)) -> &Self::Output {
        self.get(i.into())
    }
}

impl<X, R> Index<(isize, isize, isize)> for XChainComplex3<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = XChainComplexSummand<X, R>;
    fn index(&self, i: (isize, isize, isize)) -> &Self::Output {
        self.get(i.into())
    }
}

#[cfg(test)]
pub(crate) mod tests { 
    use num_traits::Zero;
    use yui::lc::Free;
    use yui_matrix::sparse::SpVec;

    use crate::{RModStr, ChainComplex, ChainComplexCommon};

    use super::*;

    type X = Free<i64>;
    fn e(i: usize) -> X { 
        X::from(i as i64)
    }

    impl From<ChainComplex<i64>> for XChainComplex<X, i64> {
        fn from(c: ChainComplex<i64>) -> Self {
            Self::new(
                c.summands().map(|s| {
                    let n = s.rank();
                    XModStr::free((0..n).map(e))
                }), 
                c.d_deg(), 
                move |i, z| {
                    let n = c[i].rank();
                    z.apply(|x| {
                        let v = SpVec::unit(n, x.0 as usize);
                        let dv = c.d(i, &v);
                        dv.iter().map(|(i, a)| (e(i), *a)).collect()
                    })
                }
            )
        }
    }

    #[test]
    fn d3() { 
        let c = XChainComplex::from(ChainComplex::d3());

        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 4);
        assert_eq!(c[3].rank(), 1);

        c.check_d_all();
    }

    #[test]
    fn s2() { 
        let c = XChainComplex::from(ChainComplex::s2());

        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 4);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();
    }

    #[test]
    fn t2() { 
        let c = XChainComplex::from(ChainComplex::t2());

        assert_eq!(c[0].rank(), 9);
        assert_eq!(c[1].rank(), 27);
        assert_eq!(c[2].rank(), 18);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();
    }

    #[test]
    fn rp2() { 
        let c = XChainComplex::from(ChainComplex::rp2());
        
        assert_eq!(c[0].rank(), 6);
        assert_eq!(c[1].rank(), 15);
        assert_eq!(c[2].rank(), 10);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();
    }

    #[test]
    fn s2_reduced() { 
        let c = XChainComplex::from(ChainComplex::s2());
        let r = c.reduced();

        assert_eq!(r[0].rank(), 1);
        assert_eq!(r[1].rank(), 0);
        assert_eq!(r[2].rank(), 1);

        r.check_d_all();

        let z = r[0].gen_chain(0);
        let dz = r.d(0, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());

        let z = r[2].gen_chain(0);
        let dz = r.d(2, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());
    }

    #[test]
    fn t2_reduced() { 
        let c = XChainComplex::from(ChainComplex::t2());
        let r = c.reduced();

        assert_eq!(r[0].rank(), 1);
        assert_eq!(r[1].rank(), 2);
        assert_eq!(r[2].rank(), 1);

        r.check_d_all();

        let z = r[0].gen_chain(0);
        let dz = r.d(0, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());

        let z1 = r[1].gen_chain(0);
        let z2 = r[1].gen_chain(1);
        let dz1 = c.d(1, &z1);
        let dz2 = c.d(1, &z2);

        assert!(!z1.is_zero());
        assert!(!z2.is_zero());
        assert!(dz1.is_zero());
        assert!(dz2.is_zero());

        let z = r[2].gen_chain(0);
        let dz = c.d(2, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());
    }

    #[test]
    fn rp2_reduced() { 
        let c = XChainComplex::from(ChainComplex::rp2());
        let r = c.reduced();
        
        assert_eq!(r[0].rank(), 1);
        assert_eq!(r[1].rank(), 1);
        assert_eq!(r[2].rank(), 1);

        r.check_d_all();

        let z = r[0].gen_chain(0);
        let dz = c.d(0, &z);
        assert!(!z.is_zero());
        assert!(dz.is_zero());

        let z = r[1].gen_chain(0);
        let dz = c.d(1, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());

        let z = r[2].gen_chain(0);
        let dz = c.d(2, &z);

        assert!(!z.is_zero());
        assert!(!dz.is_zero());
        assert_eq!(r[1].vectorize(&dz).to_dense()[0].abs(), 2);
    }
}