use std::ops::Index;

use delegate::delegate;
use itertools::Itertools;
use yui_core::{EucRing, EucRingOps, Ring, RingOps, Deg, isize2, isize3};
use yui_matrix::sparse::{SpVec, Trans};

use crate::{DisplayAt, GridBase, GridIter};

use super::grid::GridTrait;
use super::complex::ChainComplexBase;

pub type Homology<R>  = HomologyBase<isize,  R>;
pub type Homology2<R> = HomologyBase<isize2, R>;
pub type Homology3<R> = HomologyBase<isize3, R>;

const ERR_NO_TRANS: &'static str = "not computed with trans.";

#[derive(Debug, Clone)]
pub struct HomologySummand<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    rank: usize,
    tors: Vec<R>,
    trans: Option<Trans<R>>
}

impl<R> HomologySummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(rank: usize, tors: Vec<R>, trans: Option<Trans<R>>) -> Self { 
        Self { rank, tors, trans }
    }

    pub fn zero() -> Self { 
        Self::new(0, vec![], None)
    }

    pub fn rank(&self) -> usize { 
        self.rank
    }

    pub fn tors(&self) -> &Vec<R> {
        &self.tors
    }

    pub fn is_zero(&self) -> bool { 
        self.rank() == 0 && self.is_free()
    }

    pub fn is_free(&self) -> bool { 
        self.tors().is_empty()
    }

    pub fn ngens(&self) -> usize { 
        let t = self.trans.as_ref().expect(ERR_NO_TRANS);
        t.tgt_dim()
    }

    pub fn gen(&self, i: usize) -> SpVec<R> {
        assert!(i < self.ngens());

        let t = self.trans.as_ref().expect(ERR_NO_TRANS);
        t.backward_mat().col_vec(i)
    }

    pub fn trans_forward(&self, v: &SpVec<R>) -> SpVec<R> { 
        let t = self.trans.as_ref().expect(ERR_NO_TRANS);
        assert_eq!(t.src_dim(), v.dim());

        t.forward(v)
    }

    pub fn trans_backward(&self, v: &SpVec<R>) -> SpVec<R> { 
        let t = self.trans.as_ref().expect(ERR_NO_TRANS);
        assert_eq!(t.tgt_dim(), v.dim());
        
        t.backward(v)
    }

    pub fn display(&self) -> Option<String> {
        use yui_utils::superscript;

        let rank = self.rank();
        let tors = self.tors().iter()
            .into_group_map_by(|r| r.to_string())
            .into_iter().map(|(k, list)| (k, list.len()))
            .collect_vec();

        if rank == 0 && tors.is_empty() { 
            return None
        }
    
        let mut res = vec![];
        let symbol = R::math_symbol();
    
        if rank > 1 {
            let str = format!("{}{}", symbol, superscript(rank as isize));
            res.push(str);
        } else if rank == 1 { 
            let str = format!("{}", symbol);
            res.push(str);
        }
        
        for (t, r) in tors.iter() { 
            let str = if r > &1 { 
                format!("({}/{}){}", symbol, t, superscript(*r as isize))
            } else { 
                format!("({}/{})", symbol, t)
            };
            res.push(str);
        }
    
        let str = res.join(" ⊕ ");
        Some(str)
    }
}

impl<R> Default for HomologySummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn default() -> Self {
        Self::zero()
    }
}

pub struct HomologyBase<I, R>
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    summands: GridBase<I, HomologySummand<R>>
}

impl<I, R> HomologyBase<I, R>
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(complex: &ChainComplexBase<I, R>, with_trans: bool) -> Self { 
        let summands = GridBase::new(
            complex.support(), 
            |i| complex.homology_at(i, with_trans)
        );
        Self { summands }
    }

    delegate! { 
        to self.summands { 
            fn get(&self, i: I) -> &HomologySummand<R>;    
        }
    }
}

impl<R> Index<isize> for HomologyBase<isize, R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = HomologySummand<R>;
    delegate! { 
        to self.summands { 
            fn index(&self, i: isize) -> &Self::Output;
        }
    }
}

impl<R> Index<(isize, isize)> for HomologyBase<isize2, R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = HomologySummand<R>;
    delegate! { 
        to self.summands { 
            fn index(&self, i: (isize, isize)) -> &Self::Output;
        }
    }
}

impl<R> Index<(isize, isize, isize)> for HomologyBase<isize3, R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = HomologySummand<R>;
    delegate! { 
        to self.summands { 
            fn index(&self, i: (isize, isize, isize)) -> &Self::Output;
        }
    }
}

impl<I, R> GridTrait<I> for HomologyBase<I, R>
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Itr = GridIter<I>;
    delegate! { 
        to self.summands { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: I) -> bool;
        }
    }
}

impl<I, R> DisplayAt<I> for HomologyBase<I, R>
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    fn display_at(&self, i: I) -> Option<String> {
        self.get(i).display()
    }
}

#[cfg(test)]
mod tests { 
    use yui_matrix::sparse::SpVec;
    use crate::{ChainComplexTrait, ChainComplex};

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
        assert_eq!(h2.ngens(), 1);

        let z = h2.gen(0);
        assert!(!z.is_zero());
        assert!(c.is_cycle(2, &z));
        assert_eq!(h2.trans_forward(&z), SpVec::from(vec![1]));
    }

    #[test]
    fn t2_gens() {
        let c = ChainComplex::<i32>::t2();
        let h = c.homology(true);

        let h2 = &h[2];
        assert_eq!(h2.ngens(), 1);

        let z = h2.gen(0);
        assert!(!z.is_zero());
        assert!(c.is_cycle(2, &z));
        assert_eq!(h2.trans_forward(&z), SpVec::from(vec![1]));
        assert_eq!(h2.ngens(), 1);

        let h1 = &h[1];
        assert_eq!(h1.ngens(), 2);

        let a = h1.gen(0);
        let b = h1.gen(1);

        assert!(!a.is_zero());
        assert!(!b.is_zero());
        assert!(c.is_cycle(1, &a));
        assert!(c.is_cycle(1, &b));
        assert_eq!(h1.trans_forward(&a), SpVec::from(vec![1, 0]));
        assert_eq!(h1.trans_forward(&b), SpVec::from(vec![0, 1]));
    }

    #[test]
    fn rp2_gens() {
        let c = ChainComplex::<i32>::rp2();
        let h = c.homology(true);

        let h1 = &h[1];
        assert_eq!(h1.ngens(), 1);

        let z = h1.gen(0);

        assert!(!z.is_zero());
        assert!(c.is_cycle(1, &z));
        assert_eq!(h1.trans_forward(&z), SpVec::from(vec![1]));
    }
}