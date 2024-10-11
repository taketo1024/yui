use derive_more::derive::Display;
use itertools::Either;
use yui::lc::Gen;
use yui::Elem;
use yui::{Ring, RingOps};
use yui_matrix::sparse::SpMat;
use yui_matrix::MatTrait;

#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Display, Debug, Default)]
#[display("e({},{})", _0, _1)]
pub struct EnumGen(isize, usize);

impl Elem for EnumGen {
    fn math_symbol() -> String {
        "E".into()
    }
}

impl Gen for EnumGen {}

pub type GenericChainComplex<R> = XChainComplexBase<isize, EnumGen, R>;

use crate::{Grid1, XChainComplexBase, XModStr};

impl<R> GenericChainComplex<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn mat<I>(shape: (usize, usize), entries: I) -> SpMat<R>
    where R: From<i32>, I: IntoIterator<Item = i32> { 
        SpMat::from_dense_data(shape, entries.into_iter().map(|x| R::from(x)))
    }

    pub fn from_mats(d_deg: isize, offset: isize, mats: Vec<SpMat<R>>) -> Self { 
        let n = mats.len() as isize;
        let range = offset .. offset + n;
        let range = if d_deg.is_positive() { 
            Either::Left(range) 
        } else { 
            Either::Right(range.rev())
        };

        let summands = Grid1::generate(range.clone(), |i| { 
            let c = (i - offset) as usize;
            let r = mats[c].ncols();
            XModStr::free((0..r).map(|j| EnumGen(i, j)))
        });
        
        Self::new(
            summands.clone(), d_deg, 
            move |i, z| {
                let c = (i - offset) as usize;
                let d = &mats[c];
                let v = summands[i].vectorize(z);
                let dv = d * v;
                let dz = summands[i + d_deg].as_chain(&dv);
                dz
            }
        )
    }
    
    pub fn one() -> GenericChainComplex<R> {
        GenericChainComplex::from_mats(-1, 0,
            vec![
                SpMat::from_dense_data((0, 1), [])
            ]
        )
    }

    pub fn one_one(r: R) -> GenericChainComplex<R> {
        GenericChainComplex::from_mats(-1, 0,
            vec![
                SpMat::from_dense_data((0, 1), []),
                SpMat::from_dense_data((1, 1), [r])
            ]
        )
    }

    pub fn two_one(r1: R, r2: R) -> GenericChainComplex<R> {
        GenericChainComplex::from_mats(-1, 0,
            vec![
                SpMat::from_dense_data((0, 1), []),
                SpMat::from_dense_data((1, 2), [r1, r2])
            ]
        )
    }

    pub fn one_two(r1: R, r2: R) -> GenericChainComplex<R> {
        GenericChainComplex::from_mats(-1, 0,
            vec![
                SpMat::from_dense_data((0, 2), []),
                SpMat::from_dense_data((2, 1), [r1, r2])
            ]
        )
    }

    pub fn d3() -> GenericChainComplex<R>
    where R: From<i32> {
        GenericChainComplex::from_mats(-1, 0,
            vec![
                Self::mat((0, 4), []),
                Self::mat((4, 6), [-1, -1, 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 1] ),
                Self::mat((6, 4), [1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1] ),
                Self::mat((4, 1), [-1, 1, -1, 1]),
            ]
        )
    }

    pub fn s2() -> GenericChainComplex<R>
    where R: From<i32> {
        GenericChainComplex::from_mats(-1, 0,
            vec![
                Self::mat((0, 4), []),
                Self::mat((4, 6), [-1, -1, 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 1]),
                Self::mat((6, 4), [1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1] ),
            ]
        )
    }

    pub fn t2() -> GenericChainComplex<R>
    where R: From<i32> {
        GenericChainComplex::from_mats(-1, 0,
            vec![
                Self::mat((0, 9), []),
                Self::mat((9, 27), [-1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, -1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1]),
                Self::mat((27, 18), [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1]),
            ]
        )
    }

    pub fn rp2() -> GenericChainComplex<R>
    where R: From<i32> { 
        GenericChainComplex::from_mats(-1, 0,
            vec![
                Self::mat((0, 6), []),
                Self::mat((6, 15), [-1, -1, 0, 0, 0, 0, 0, -1, -1, 0, -1, 0, 0, 0, 0, 1, 0, -1, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 1, 1, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1] ),
                Self::mat((15, 10), [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1] ),
            ]
        )
    }
}

#[cfg(test)]
pub(crate) mod tests { 
    use num_traits::Zero;

    use crate::{ChainComplexTrait, RModStr};

    use super::*;

    #[test]
    fn zero() { 
        let c = GenericChainComplex::<i32>::zero();
        assert_eq!(c[0].rank(), 0);
        
        c.check_d_all();

        let h = c.homology(false);
        assert!(h[0].is_zero());
    }

    #[test]
    fn single() { 
        let c = GenericChainComplex::<i32>::one();
        assert_eq!(c[0].rank(), 1);

        c.check_d_all();

        let h = c.homology(false);
        
        assert_eq!(h[0].rank(), 1);
        assert!( h[0].is_free());
        assert!(!h[0].is_zero());
    }

    #[test]
    fn one_to_one() { 
        let c = GenericChainComplex::<i32>::one_one(1);
        let h = c.homology(false);

        assert!(h[0].is_zero());
        assert!(h[1].is_zero());
    }

    #[test]
    fn two_to_one() { 
        let c = GenericChainComplex::<i32>::two_one(1, -1);
        let h = c.homology(false);

        assert!(h[0].is_zero());
        assert_eq!(h[1].rank(), 1);
        assert!(h[1].is_free());
    }

    #[test]
    fn one_to_two() { 
        let c = GenericChainComplex::<i32>::one_two(1, -1);
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());
        assert!(h[1].is_zero());
    }

    #[test]
    fn torsion() { 
        let c = GenericChainComplex::<i32>::one_one(2);
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[0].tors(), &vec![2]);
        assert!(!h[0].is_free());
    }

    #[test]
    fn d3() {
        let c = GenericChainComplex::<i32>::d3();

        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 4);
        assert_eq!(c[3].rank(), 1);

        c.check_d_all();

        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 0);
        assert!(h[1].is_free());

        assert_eq!(h[2].rank(), 0);
        assert!(h[2].is_free());

        assert_eq!(h[3].rank(), 0);
        assert!(h[3].is_free());
    }

    #[test]
    fn s2() {
        let c = GenericChainComplex::<i32>::s2();

        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 4);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();

        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 0);
        assert!(h[1].is_free());

        assert_eq!(h[2].rank(), 1);
        assert!(h[2].is_free());
    }

    #[test]
    fn t2() {
        let c = GenericChainComplex::<i32>::t2();

        assert_eq!(c[0].rank(), 9);
        assert_eq!(c[1].rank(), 27);
        assert_eq!(c[2].rank(), 18);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();

        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 2);
        assert!(h[1].is_free());

        assert_eq!(h[2].rank(), 1);
        assert!(h[2].is_free());
    }

    #[test]
    fn rp2() {
        let c = GenericChainComplex::<i32>::rp2();

        assert_eq!(c[0].rank(), 6);
        assert_eq!(c[1].rank(), 15);
        assert_eq!(c[2].rank(), 10);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();
        
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].tors(), &vec![2]);
        assert!(!h[1].is_free());

        assert_eq!(h[2].rank(), 0);
        assert!(h[2].is_free());
    }

    #[test]
    fn s2_gens() {
        let c = GenericChainComplex::<i32>::s2();
        let h = c.homology(true);

        let z = h[2].gen_chain(0);
        let dz = c.d(2, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());
    }

    #[test]
    fn t2_gens() {
        let c = GenericChainComplex::<i32>::t2();
        let h = c.homology(true);

        let z = h[2].gen_chain(0);
        let dz = c.d(2, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());

        let a = h[1].gen_chain(0);
        let b = h[1].gen_chain(1);
        let da = c.d(1, &a);
        let db = c.d(1, &b);

        assert!(!a.is_zero());
        assert!(!b.is_zero());
        assert!(da.is_zero());
        assert!(db.is_zero());
    }

    #[test]
    fn rp2_gens() {
        let c = GenericChainComplex::<i32>::rp2();
        let h = c.homology(true);

        let z = h[1].gen_chain(0);
        let dz = c.d(1, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());

        assert!(!h[1].vectorize_euc(&z).is_zero());
        assert!(h[1].vectorize_euc(&(z * 2)).is_zero()); // order 2
    }
}