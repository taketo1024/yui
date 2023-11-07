use yui_core::{Deg, isize2, isize3, EucRing, EucRingOps};
use yui_lin_comb::Gen;

use crate::{GridBase, XModStr, XChainComplexBase, ComputeHomology, GridTrait, RModStr};

pub type XHomologySummand<X, R> = XModStr<X, R>;
pub type XHomologyBase<I, X, R> = GridBase<I, XHomologySummand<X, R>>;

pub type XHomology <X, R> = XHomologyBase<isize,  X, R>;
pub type XHomology2<X, R> = XHomologyBase<isize2, X, R>;
pub type XHomology3<X, R> = XHomologyBase<isize3, X, R>;

impl<I, X, R> XChainComplexBase<I, X, R>
where 
    I: Deg,
    X: Gen,
    R: EucRing, for<'x> &'x R: EucRingOps<R>
{
    pub fn x_homology(&self) -> XHomologyBase<I, X, R> {
        let h = self.homology(true);

        XHomologyBase::new(h.support(), move |i| { 
            let ci = &self[i];
            let hi = &h[i];

            let t1 = hi.trans().unwrap().clone();
            let t = if let Some(t0) = ci.trans() { 
                t0.compose(&t1)
            } else { 
                t1
            };

            XModStr::new(
                ci.gens().clone(), 
                hi.rank(), 
                hi.tors().clone(), 
                Some(t)
            )
        })
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use num_traits::Zero;
    use yui_matrix::sparse::SpVec;

    use crate::{ChainComplex, XChainComplex, RModStr};
 
    #[test]
    fn zero() { 
        let c = XChainComplex::from(ChainComplex::zero());
        let h = c.x_homology();
        
        assert!(h[0].is_zero());
    }

    #[test]
    fn single() { 
        let c = XChainComplex::from(ChainComplex::one());
        let h = c.x_homology();
        
        assert_eq!(h[0].rank(), 1);
        assert!( h[0].is_free());
        assert!(!h[0].is_zero());
    }

    #[test]
    fn one_to_one() { 
        let c = XChainComplex::from(ChainComplex::one_one(1));
        let h = c.x_homology();

        assert!(h[0].is_zero());
        assert!(h[1].is_zero());
    }

    #[test]
    fn two_to_one() { 
        let c = XChainComplex::from(ChainComplex::two_one(1, -1));
        let h = c.x_homology();

        assert!(h[0].is_zero());
        assert_eq!(h[1].rank(), 1);
        assert!(h[1].is_free());
    }

    #[test]
    fn one_to_two() { 
        let c = XChainComplex::from(ChainComplex::one_two(1, -1));
        let h = c.x_homology();

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());
        assert!(h[1].is_zero());
    }

    #[test]
    fn torsion() { 
        let c = XChainComplex::from(ChainComplex::one_one(2));
        let h = c.x_homology();

        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[0].tors(), &vec![2]);
        assert!(!h[0].is_free());
    }

    #[test]
    fn d3() {
        let c = XChainComplex::from(ChainComplex::d3());
        let h = c.x_homology();

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[2].is_free(), true);

        assert_eq!(h[3].rank(), 0);
        assert_eq!(h[3].is_free(), true);
    }

    #[test]
    fn s2() {
        let c = XChainComplex::from(ChainComplex::s2());
        let h = c.x_homology();

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);

        let h2 = &h[2];
        let z = h2.gen_chain(0);

        assert!(!z.is_zero());
        assert!(c.d(2, &z).is_zero());
        assert_eq!(h2.vectorize(&z), SpVec::from(vec![1]));
    }

    #[test]
    fn t2() {
        let c = XChainComplex::from(ChainComplex::t2());
        let h = c.x_homology();

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);

        let h2 = &h[2];
        let z = h2.gen_chain(0);

        assert!(!z.is_zero());
        assert!(c.d(2, &z).is_zero());
        assert_eq!(h2.vectorize(&z), SpVec::from(vec![1]));

        let h1 = &h[1];
        let a = h1.gen_chain(0);
        let b = h1.gen_chain(1);

        assert!(!a.is_zero());
        assert!(!b.is_zero());
        assert!(c.d(1, &a).is_zero());
        assert!(c.d(1, &b).is_zero());
        assert_eq!(h1.vectorize(&a).reduced(), SpVec::from(vec![1, 0]));
        assert_eq!(h1.vectorize(&b).reduced(), SpVec::from(vec![0, 1]));
    }

    #[test]
    fn rp2() {
        let c = XChainComplex::from(ChainComplex::rp2());
        let h = c.x_homology();

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].tors(), &vec![2]);
        assert_eq!(h[1].is_free(), false);

        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[2].is_free(), true);

        let h1 = &h[1];
        let z = h1.gen_chain(0);

        assert!(!z.is_zero());
        assert!(c.d(1, &z).is_zero());
        assert_eq!(h1.vectorize(&z), SpVec::from(vec![1]));
    }
}