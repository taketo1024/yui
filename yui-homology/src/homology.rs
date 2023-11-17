use yui::{EucRing, EucRingOps};

use crate::utils::HomologyCalc;
use crate::{GridDeg, isize2, isize3};
use crate::{Grid, GridTrait, ChainComplexTrait, SimpleRModStr, ChainComplexBase};

pub type HomologySummand<R> = SimpleRModStr<R>;
pub type HomologyBase<I, R> = Grid<I, HomologySummand<R>>;

pub type Homology<R>  = HomologyBase<isize,  R>;
pub type Homology2<R> = HomologyBase<isize2, R>;
pub type Homology3<R> = HomologyBase<isize3, R>;

pub trait ComputeHomology<I> {
    type Output;
    fn compute_homology(&self, i: I, with_trans: bool) -> Self::Output;
}

impl<I, R, C> ComputeHomology<I> for C 
where 
    I: GridDeg, 
    R: EucRing, for<'x> &'x R: EucRingOps<R>,
    C: ChainComplexTrait<I, R = R>
{
    type Output = HomologySummand<R>;

    fn compute_homology(&self, i: I, with_trans: bool) -> HomologySummand<R> {
        let i0 = i - self.d_deg();
        let d0 = self.d_matrix(i0);
        let d1 = self.d_matrix(i );
        HomologyCalc::calculate(d0, d1, with_trans)
    }
}

impl<I, R> ChainComplexBase<I, R> 
where I: GridDeg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology_at(&self, i: I, with_trans: bool) -> HomologySummand<R> {
        let c = &self[i];
        let h = self.compute_homology(i, with_trans);
        c.compose(&h)
    }

    pub fn homology(&self, with_trans: bool) -> HomologyBase<I, R> {
        HomologyBase::generate(
            self.support(), 
            |i| self.homology_at(i, with_trans)
        )
    }
}

#[cfg(test)]
mod tests { 
    use yui_matrix::sparse::SpVec;
    use crate::{ChainComplex, RModStr, ChainComplexTrait};

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
        let c = ChainComplex::<i32>::s2();
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
        let c = ChainComplex::<i32>::t2();
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
        let c = ChainComplex::<i32>::rp2();
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
        let c = ChainComplex::<i32>::s2();
        let h = c.homology(true);

        let t = h[2].trans().unwrap();
        let z = h[2].gen_vec(0);
        let dz = c.d(2, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());
        assert_eq!(t.forward(&z), SpVec::from(vec![1]));
    }

    #[test]
    fn t2_gens() {
        let c = ChainComplex::<i32>::t2();
        let h = c.homology(true);

        let t2 = h[2].trans().unwrap();
        let z = h[2].gen_vec(0);
        let dz = c.d(2, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());
        assert_eq!(t2.forward(&z), SpVec::from(vec![1]));

        let t1 = h[1].trans().unwrap();
        let a = h[1].gen_vec(0);
        let b = h[1].gen_vec(1);
        let da = c.d(1, &a);
        let db = c.d(1, &b);

        assert!(!a.is_zero());
        assert!(!b.is_zero());
        assert!(da.is_zero());
        assert!(db.is_zero());
        assert_eq!(t1.forward(&a), SpVec::from(vec![1, 0]));
        assert_eq!(t1.forward(&b), SpVec::from(vec![0, 1]));
    }

    #[test]
    fn rp2_gens() {
        let c = ChainComplex::<i32>::rp2();
        let h = c.homology(true);

        let t = h[1].trans().unwrap();
        let z = h[1].gen_vec(0);
        let dz = c.d(1, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());
        assert_eq!(t.forward(&z), SpVec::from(vec![1]));

        let v = SpVec::from(vec![-1,1,-1,1,1,1,1,1,1,-1]);
        let dv = c.d(2, &v);
        let w = h[1].trans().unwrap().forward(&dv);

        assert!(!dv.is_zero());
        assert_eq!(w.to_dense()[0].abs(), 2);
    }
}