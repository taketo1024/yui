use std::ops::Index;
use std::collections::HashMap;
use std::fmt::Display;

use itertools::Itertools;
use yui_core::{EucRing, EucRingOps, Ring, RingOps, Deg, isize2, isize3};
use yui_matrix::sparse::SpVec;

use super::trans::Trans;
use super::graded::Graded;
use super::complex::ChainComplexBase;

pub type Homology<R>  = HomologyBase<isize,  R>;
pub type Homology2<R> = HomologyBase<isize2, R>;
pub type Homology3<R> = HomologyBase<isize3, R>;

#[derive(Debug, Clone)]
pub struct HomologySummand<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    rank: usize,
    tors: Vec<R>,
    trans: Option<Trans<R>>
}

impl<R> HomologySummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(rank: usize, tors: Vec<R>, trans: Trans<R>) -> Self { 
        let trans = Some(trans);
        Self { rank, tors, trans }
    }

    pub fn new_no_trans(rank: usize, tors: Vec<R>) -> Self { 
        Self { rank, tors, trans: None }
    }

    fn zero() -> Self { 
        Self::new_no_trans(0, vec![])
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

    pub fn c2h(&self, v: &SpVec<R>) -> SpVec<R> { 
        let t = self.trans.as_ref().expect("not computed with gens.");
        assert_eq!(t.src_dim(), v.dim());

        t.forward(v)
    }

    pub fn h2c(&self, v: &SpVec<R>) -> SpVec<R> { 
        let t = self.trans.as_ref().expect("not computed with gens.");
        assert_eq!(t.tgt_dim(), v.dim());
        
        t.backward(v)
    }

    pub fn ngens(&self) -> usize { 
        let t = self.trans.as_ref().expect("not computed with gens.");
        t.tgt_dim()
    }

    pub fn gen(&self, i: usize) -> SpVec<R> {
        assert!(i < self.ngens());

        let t = self.trans.as_ref().expect("not computed with gens.");
        t.b_mat().col_vec(i)
    }

    pub fn print(&self) {
        println!("{}", self.to_string())
    }
}

impl<R> Display for HomologySummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use yui_utils::superscript;

        let rank = self.rank();
        let tors = self.tors().iter()
            .into_group_map_by(|r| r.to_string())
            .into_iter().map(|(k, list)| (k, list.len()))
            .collect_vec();

        if rank == 0 && tors.is_empty() { 
            return f.write_str(".")
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
        f.write_str(&str) 
    }
}

pub struct HomologyBase<I, R>
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    support: Vec<I>,
    summands: HashMap<I, HomologySummand<R>>,
    zero_smd: HomologySummand<R>
}

impl<I, R> HomologyBase<I, R>
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(complex: &ChainComplexBase<I, R>, with_gens: bool) -> Self { 
        let support = complex.support().collect_vec();
        let summands = HashMap::from_iter(
            support.iter().map(|&i| {
                let h_i = complex.homology_at(i, with_gens);
                (i, h_i)
            })
        );
        let zero_smd = HomologySummand::zero();
        Self { support, summands, zero_smd }
    }

    pub fn get(&self, i: I) -> &HomologySummand<R> {
        self.summands.get(&i).unwrap_or(&self.zero_smd)
    }
}

impl<R> Index<isize> for HomologyBase<isize, R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = HomologySummand<R>;
    fn index(&self, i: isize) -> &Self::Output {
        self.get(i)
    }
}

impl<R> Index<(isize, isize)> for HomologyBase<isize2, R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = HomologySummand<R>;
    fn index(&self, i: (isize, isize)) -> &Self::Output {
        self.get(isize2(i.0, i.1))
    }
}

impl<R> Index<(isize, isize, isize)> for HomologyBase<isize3, R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = HomologySummand<R>;
    fn index(&self, i: (isize, isize, isize)) -> &Self::Output {
        self.get(isize3(i.0, i.1, i.2))
    }
}

impl<I, R> Graded<I> for HomologyBase<I, R>
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Itr = std::vec::IntoIter<I>;
    fn support(&self) -> Self::Itr {
        self.support.clone().into_iter()
    }

    fn display(&self, i: I) -> String {
        self.get(i).to_string()
    }
}

#[cfg(test)]
mod tests { 
    use yui_matrix::sparse::SpMat;
    use yui_matrix::sparse::SpVec;

    use super::super::complex::ChainComplex;
    use super::super::complex::tests::*;

    #[test]
    fn singleton() { 
        let c = ChainComplex::<i32>::from_mats(-1, 0,
            vec![
                SpMat::from_vec((0, 1), vec![]),
                SpMat::from_vec((1, 0), vec![])
            ]
        );        
        let h = c.homology(false);
        
        assert_eq!(h[0].rank(), 1);
        assert!( h[0].is_free());
        assert!(!h[0].is_zero());
    }

    #[test]
    fn torsion() { 
        let c = ChainComplex::<i32>::from_mats(-1, 0,
            vec![
                SpMat::from_vec((0, 1), vec![]),
                SpMat::from_vec((1, 1), vec![2]),
                SpMat::from_vec((1, 0), vec![]),
            ]
        );
        let h = c.homology(false);
        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[0].tors(), &vec![2]);
        assert!(!h[0].is_free());

    }

    #[test]
    fn d3() {
        let c = Samples::<i32>::d3();
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
        let c = Samples::<i32>::s2();
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
        let c = Samples::<i32>::t2();
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
        let c = Samples::<i32>::rp2();
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
        let c = Samples::<i32>::s2();
        let h = c.homology(true);

        let h2 = &h[2];
        assert_eq!(h2.ngens(), 1);

        let z = h2.gen(0);
        assert!(!z.is_zero());
        assert!(c.is_cycle(2, &z));
        assert_eq!(h2.c2h(&z), SpVec::from(vec![1]));
    }

    #[test]
    fn t2_gens() {
        let c = Samples::<i32>::t2();
        let h = c.homology(true);

        let h2 = &h[2];
        assert_eq!(h2.ngens(), 1);

        let z = h2.gen(0);
        assert!(!z.is_zero());
        assert!(c.is_cycle(2, &z));
        assert_eq!(h2.c2h(&z), SpVec::from(vec![1]));
        assert_eq!(h2.ngens(), 1);

        let h1 = &h[1];
        assert_eq!(h1.ngens(), 2);

        let a = h1.gen(0);
        let b = h1.gen(1);

        assert!(!a.is_zero());
        assert!(!b.is_zero());
        assert!(c.is_cycle(1, &a));
        assert!(c.is_cycle(1, &b));
        assert_eq!(h1.c2h(&a), SpVec::from(vec![1, 0]));
        assert_eq!(h1.c2h(&b), SpVec::from(vec![0, 1]));
    }

    #[test]
    fn rp2_gens() {
        let c = Samples::<i32>::rp2();
        let h = c.homology(true);

        let h1 = &h[1];
        assert_eq!(h1.ngens(), 1);

        let z = h1.gen(0);

        assert!(!z.is_zero());
        assert!(c.is_cycle(1, &z));
        assert_eq!(h1.c2h(&z), SpVec::from(vec![1]));
    }
}