use std::ops::Index;
use std::sync::Arc;

use delegate::delegate;

use itertools::Itertools;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use yui_core::{Ring, RingOps, Deg, isize2, isize3};
use yui_lin_comb::{Gen, LinComb};
use yui_matrix::sparse::{SpMat, Trans, MatType};

use crate::utils::ChainReducer;
use crate::{Grid, GridIter, ChainComplexDisplay, XModStr, RModStr, SimpleRModStr};

use super::grid::GridTrait;
use super::complex::ChainComplexTrait;

pub type XChainComplex <X, R> = XChainComplexBase<isize,  X, R>;
pub type XChainComplex2<X, R> = XChainComplexBase<isize2, X, R>;
pub type XChainComplex3<X, R> = XChainComplexBase<isize3, X, R>;

pub type XChainComplexSummand<X, R> = XModStr<X, R>;

pub struct XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    summands: Grid<I, XChainComplexSummand<X, R>>,
    d_deg: I,
    d_map: Arc<dyn Fn(I, &LinComb<X, R>) -> LinComb<X, R> + Send + Sync>,
}

impl<I, X, R> XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    pub fn new<F>(summands: Grid<I, XChainComplexSummand<X, R>>, d_deg: I, d_map: F) -> Self
    where F: Fn(I, &LinComb<X, R>) -> LinComb<X, R> + Send + Sync + 'static {
        let d_map = Arc::new(d_map);
        Self { summands, d_deg, d_map }
    }

    pub fn d(&self, i: I, z: &LinComb<X, R>) -> LinComb<X, R> { 
        (self.d_map)(i, z)
    }

    fn d_matrix_for(&self, i: I, q: &SpMat<R>) -> SpMat<R> { 
        assert_eq!(q.rows(), self[i].rank());

        let i1 = i + self.d_deg;
        let (m, n) = (self[i1].rank(), q.cols());

        let entries = (0..n).into_par_iter().flat_map(|j| {
            let v = q.col_vec(j);
            let z = v.iter().map(|(k, a)|
                self[i].gen_chain(k) * a
            ).sum::<LinComb<_, _>>();

            let dz = self.d(i, &z);
            let w = self[i1].vectorize(&dz);

            w.iter().filter_map(|(i, a)| { 
                if !a.is_zero() {
                    Some((i, j, a.clone()))
                } else { 
                    None
                }
            }).collect_vec()
        }).collect::<Vec<_>>();

        SpMat::from_entries((m, n), entries)
    }

    pub fn reduced(&self) -> XChainComplexBase<I, X, R> { 
        let mut reducer = ChainReducer::new(self.support(), self.d_deg, true);

        for i in self.support() {
            let d = if let Some(t) = reducer.trans(i) {
                self.d_matrix_for(i, t.backward_mat())
            } else { 
                self.d_matrix(i)
            };
            reducer.set_matrix(i, d);
            reducer.reduce_at(i);
        }

        self.reduced_by(|i| 
            reducer.take_trans(i).unwrap()
        )
    }

    pub fn reduced_by<F>(&self, mut trans_map: F) -> XChainComplexBase<I, X, R>
    where F: FnMut(I) -> Trans<R> { 
        let d_deg = self.d_deg;
        let d_map = self.d_map.clone();

        let summands = self.summands.map(|i, summand| { 
            let t = trans_map(i);

            assert_eq!(t.src_dim(), self[i].rank());

            let r = t.tgt_dim();
            let s = SimpleRModStr::new(r, vec![], Some(t));

            summand.compose(&s)
        });

        Self { summands, d_deg, d_map }
    }
}

impl<I, X, R> GridTrait<I> for XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    type Itr = GridIter<I>;
    type E = XChainComplexSummand<X, R>;
    
    delegate! { 
        to self.summands { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: I) -> bool;
            fn get(&self, i: I) -> &Self::E;
        }
    }
}

impl<I, X, R> ChainComplexTrait<I> for XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    type R = R;

    fn rank(&self, i: I) -> usize {
        self[i].rank()
    }
    
    fn d_deg(&self) -> I {
        self.d_deg
    }

    fn d_matrix(&self, i: I) -> SpMat<Self::R> { 
        let n = self[i].rank();
        let q = SpMat::id(n);
        self.d_matrix_for(i, &q)
    }
}

impl<I, X, R> ChainComplexDisplay<I> for XChainComplexBase<I, X, R>
where I: Deg, X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {}

impl<I, X, R> Index<I> for XChainComplexBase<I, X, R>
where I: Deg, X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
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
    use yui_lin_comb::Free;
    use yui_matrix::sparse::SpVec;

    use crate::{RModStr, ChainComplex, Grid1};

    use super::*;

    type X = Free<i64>;
    fn e(i: usize) -> X { 
        X::from(i as i64)
    }

    impl From<ChainComplex<i64>> for XChainComplex<X, i64> {
        fn from(c: ChainComplex<i64>) -> Self {
            let summands = Grid1::generate(c.support(), |i| { 
                let n = c[i].rank();
                let gens = (0..n).map(|j| e(j));
                XModStr::free(gens)
            });

            Self::new(summands, c.d_deg(), move |i, z| {
                let n = c[i].rank();
                z.apply(|x| {
                    let v = SpVec::unit(n, x.0 as usize);
                    let dv = c.d(i, &v);
                    dv.iter().map(|(i, a)| (e(i), a.clone())).collect()
                })
            })
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