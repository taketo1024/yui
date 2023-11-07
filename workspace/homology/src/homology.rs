use yui_core::{EucRing, EucRingOps, Deg, isize2, isize3};

use crate::utils::HomologyCalc;
use crate::{GridBase, ChainComplexTrait, RModStr, SimpleRModStr};

pub type HomologySummand<R> = SimpleRModStr<R>;
pub type HomologyBase<I, R> = GridBase<I, HomologySummand<R>>;

pub type Homology<R>  = HomologyBase<isize,  R>;
pub type Homology2<R> = HomologyBase<isize2, R>;
pub type Homology3<R> = HomologyBase<isize3, R>;

pub trait ComputeHomology<H> {
    fn homology(&self, with_trans: bool) -> H;
}

impl<I, R, C> ComputeHomology<HomologyBase<I, R>> for C 
where 
    I: Deg, 
    R: EucRing, for<'x> &'x R: EucRingOps<R>,
    C: ChainComplexTrait<I, R = R>,
    C::E: RModStr<R = R>
{
    fn homology(&self, with_trans: bool) -> HomologyBase<I, R> {
        HomologyBase::new(
            self.support(), 
            |i| {
                let i0 = i - self.d_deg();
                let d0 = self.d_matrix(i0);
                let d1 = self.d_matrix(i );
                HomologyCalc::calculate(&d0, &d1, with_trans)
            }
        )
    }
}

#[cfg(test)]
mod tests { 
    use yui_matrix::sparse::SpVec;
    use crate::homology::ComputeHomology;
    use crate::{ChainComplex, RModStr};

    #[test]
    fn zero() { 
        let c = ChainComplex::<i32>::zero();
        let h = c.homology(false);
        
        assert!(h[0].is_zero());
    }

    #[test]
    fn single() { 
        let c = ChainComplex::<i32>::one();
        let h = c.homology(false);
        
        assert_eq!(h[0].rank(), 1);
        assert!( h[0].is_free());
        assert!(!h[0].is_zero());
    }

    #[test]
    fn one_to_one() { 
        let c = ChainComplex::<i32>::one_one(1);
        let h = c.homology(false);

        assert!(h[0].is_zero());
        assert!(h[1].is_zero());
    }

    #[test]
    fn two_to_one() { 
        let c = ChainComplex::<i32>::two_one(1, -1);
        let h = c.homology(false);

        assert!(h[0].is_zero());
        assert_eq!(h[1].rank(), 1);
        assert!(h[1].is_free());
    }

    #[test]
    fn one_to_two() { 
        let c = ChainComplex::<i32>::one_two(1, -1);
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());
        assert!(h[1].is_zero());
    }

    #[test]
    fn torsion() { 
        let c = ChainComplex::<i32>::one_one(2);
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[0].tors(), &vec![2]);
        assert!(!h[0].is_free());
    }

    #[test]
    fn d3() {
        let c = ChainComplex::<i32>::d3();
        let h = c.homology(false);

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
        let c = ChainComplex::<i32>::s2();
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn t2() {
        let c = ChainComplex::<i32>::t2();
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[1].is_free(), true);

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn rp2() {
        let c = ChainComplex::<i32>::rp2();
        let h = c.homology(false);

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].tors(), &vec![2]);
        assert_eq!(h[1].is_free(), false);

        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[2].is_free(), true);
    }

    #[test]
    fn s2_gens() {
        let c = ChainComplex::<i32>::s2();
        let h = c.homology(true);

        let h2 = &h[2];
        let t = h2.trans().unwrap();
        let z = h2.gen_vec(0).unwrap();

        assert!(!z.is_zero());
        assert!(c.d(2, &z).is_zero());
        assert_eq!(t.forward(&z), SpVec::from(vec![1]));
    }

    #[test]
    fn t2_gens() {
        let c = ChainComplex::<i32>::t2();
        let h = c.homology(true);

        let h2 = &h[2];
        let t2 = h2.trans().unwrap();
        let z = h2.gen_vec(0).unwrap();

        assert!(!z.is_zero());
        assert!(c.d(2, &z).is_zero());
        assert_eq!(t2.forward(&z), SpVec::from(vec![1]));

        let h1 = &h[1];
        let t1 = h1.trans().unwrap();
        let a = h1.gen_vec(0).unwrap();
        let b = h1.gen_vec(1).unwrap();

        assert!(!a.is_zero());
        assert!(!b.is_zero());
        assert!(c.d(1, &a).is_zero());
        assert!(c.d(1, &b).is_zero());
        assert_eq!(t1.forward(&a), SpVec::from(vec![1, 0]));
        assert_eq!(t1.forward(&b), SpVec::from(vec![0, 1]));
    }

    #[test]
    fn rp2_gens() {
        let c = ChainComplex::<i32>::rp2();
        let h = c.homology(true);

        let h1 = &h[1];
        let t = h1.trans().unwrap();
        let z = h1.gen_vec(0).unwrap();

        assert!(!z.is_zero());
        assert!(c.d(1, &z).is_zero());
        assert_eq!(t.forward(&z), SpVec::from(vec![1]));
    }
}