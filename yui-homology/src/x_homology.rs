use yui::{EucRing, EucRingOps};
use yui::lc::Gen;

use crate::{GridDeg, isize2, isize3, Grid, XModStr, XChainComplexBase, GridTrait, ComputeHomology};

pub type XHomologySummand<X, R> = XModStr<X, R>;
pub type XHomologyBase<I, X, R> = Grid<I, XHomologySummand<X, R>>;

pub type XHomology <X, R> = XHomologyBase<isize,  X, R>;
pub type XHomology2<X, R> = XHomologyBase<isize2, X, R>;
pub type XHomology3<X, R> = XHomologyBase<isize3, X, R>;

impl<I, X, R> XChainComplexBase<I, X, R>
where 
    I: GridDeg,
    X: Gen,
    R: EucRing, for<'x> &'x R: EucRingOps<R>
{
    pub fn homology_at(&self, i: I, with_trans: bool) -> XHomologySummand<X, R> {
        let ci = &self[i];
        let hi = self.compute_homology(i, with_trans);
        ci.compose(&hi)
    }

    pub fn homology(&self, with_trans: bool) -> XHomologyBase<I, X, R> {
        XHomologyBase::generate(
            self.support(), 
            |i| self.homology_at(i, with_trans)
        )
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use num_traits::Zero;
    use yui_matrix::sparse::SpVec;

    use crate::{ChainComplex, XChainComplex, RModStr, ChainComplexTrait};
 
    #[test]
    fn zero() { 
        let c = XChainComplex::from(ChainComplex::zero());
        let h = c.homology(false);
        
        assert!(h[0].is_zero());
    }

    #[test]
    fn single() { 
        let c = XChainComplex::from(ChainComplex::one());
        let h = c.homology(false);
        
        assert_eq!(h[0].rank(), 1);
        assert!( h[0].is_free());
        assert!(!h[0].is_zero());
    }

    #[test]
    fn one_to_one() { 
        let c = XChainComplex::from(ChainComplex::one_one(1));
        let h = c.homology(false);

        assert!(h[0].is_zero());
        assert!(h[1].is_zero());
    }

    #[test]
    fn two_to_one() { 
        let c = XChainComplex::from(ChainComplex::two_one(1, -1));
        let h = c.homology(false);

        assert!(h[0].is_zero());
        assert_eq!(h[1].rank(), 1);
        assert!(h[1].is_free());
    }

    #[test]
    fn one_to_two() { 
        let c = XChainComplex::from(ChainComplex::one_two(1, -1));
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());
        assert!(h[1].is_zero());
    }

    #[test]
    fn torsion() { 
        let c = XChainComplex::from(ChainComplex::one_one(2));
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[0].tors(), &vec![2]);
        assert!(!h[0].is_free());
    }

    #[test]
    fn d3() {
        let c = XChainComplex::from(ChainComplex::d3());
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
        let c = XChainComplex::from(ChainComplex::s2());
        let h = c.homology(true);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 0);
        assert!(h[1].is_free());

        assert_eq!(h[2].rank(), 1);
        assert!(h[2].is_free());

        let z = h[2].gen_chain(0);
        let dz = c.d(2, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());
        assert_eq!(h[2].vectorize(&z), SpVec::from(vec![1]));
    }

    #[test]
    fn t2() {
        let c = XChainComplex::from(ChainComplex::t2());
        let h = c.homology(true);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 2);
        assert!(h[1].is_free());

        assert_eq!(h[2].rank(), 1);
        assert!(h[2].is_free());

        let z = h[2].gen_chain(0);
        let dz = c.d(2, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());
        assert_eq!(h[2].vectorize(&z), SpVec::from(vec![1]));

        let a = h[1].gen_chain(0);
        let b = h[1].gen_chain(1);
        let da = c.d(1, &a);
        let db = c.d(1, &b);

        assert!(!a.is_zero());
        assert!(!b.is_zero());
        assert!(da.is_zero());
        assert!(db.is_zero());
        assert_eq!(h[1].vectorize(&a).reduced(), SpVec::from(vec![1, 0]));
        assert_eq!(h[1].vectorize(&b).reduced(), SpVec::from(vec![0, 1]));
    }

    #[test]
    fn rp2() {
        let c = XChainComplex::from(ChainComplex::rp2());
        let h = c.homology(true);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].tors(), &vec![2]);
        assert!(!h[1].is_free());

        assert_eq!(h[2].rank(), 0);
        assert!(h[2].is_free());

        let z = h[1].gen_chain(0);
        let dz = c.d(1, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());
        assert_eq!(h[1].vectorize(&z), SpVec::from(vec![1]));

        let v = SpVec::from(vec![-1,1,-1,1,1,1,1,1,1,-1]);
        let a = c[2].as_chain(&v);
        let da = c.d(2, &a);
        
        assert!(!da.is_zero());
        assert_eq!(h[1].vectorize(&da).to_dense()[0].abs(), 2);
    }

    #[test]
    fn s2_reduced() {
        let c = XChainComplex::from(ChainComplex::s2()).reduced();
        let h = c.homology(true);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 0);
        assert!(h[1].is_free());

        assert_eq!(h[2].rank(), 1);
        assert!(h[2].is_free());

        let z = h[2].gen_chain(0);
        let dz = c.d(2, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());
        assert_eq!(h[2].vectorize(&z), SpVec::from(vec![1]));
    }

    #[test]
    fn t2_reduced() {
        let c = XChainComplex::from(ChainComplex::t2()).reduced();
        let h = c.homology(true);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 2);
        assert!(h[1].is_free());

        assert_eq!(h[2].rank(), 1);
        assert!(h[2].is_free());

        let z = h[2].gen_chain(0);
        let dz = c.d(2, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());
        assert_eq!(h[2].vectorize(&z), SpVec::from(vec![1]));

        let a = h[1].gen_chain(0);
        let b = h[1].gen_chain(1);
        let da = c.d(1, &a);
        let db = c.d(1, &b);

        assert!(!a.is_zero());
        assert!(!b.is_zero());
        assert!(da.is_zero());
        assert!(db.is_zero());
        assert_eq!(h[1].vectorize(&a).reduced(), SpVec::from(vec![1, 0]));
        assert_eq!(h[1].vectorize(&b).reduced(), SpVec::from(vec![0, 1]));
    }

    #[test]
    fn rp2_reduced() {
        let c = XChainComplex::from(ChainComplex::rp2()).reduced();
        let h = c.homology(true);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].tors(), &vec![2]);
        assert!(!h[1].is_free());

        assert_eq!(h[2].rank(), 0);
        assert!(h[2].is_free());

        let z = h[1].gen_chain(0);
        let dz = c.d(1, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());
        assert_eq!(h[1].vectorize(&z), SpVec::from(vec![1]));

        let z = h[1].gen_chain(0);
        let dz = c.d(1, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());
        assert_eq!(h[1].vectorize(&z), SpVec::from(vec![1]));

        assert_eq!(c[2].rank(), 1);
        let a = c[2].gen_chain(0);
        let da = c.d(2, &a);

        assert!(!da.is_zero());
        assert_eq!(h[1].vectorize(&da).to_dense()[0].abs(), 2);
    }
}