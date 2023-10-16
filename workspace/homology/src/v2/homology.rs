use std::ops::Index;
use std::cell::OnceCell;
use std::collections::HashMap;
use std::fmt::Display;

use itertools::Itertools;
use yui_core::{EucRing, EucRingOps, Ring, RingOps};
use yui_matrix::sparse::{SpVec, SpMat, MatType};

use crate::RModStr;
use crate::utils::HomologyCalc;

use super::complex::ChainComplexBase;
use super::deg::{Deg, isize2, isize3};
use super::graded::Graded;

pub type Homology<R>  = HomologyBase<isize,  R>;
pub type Homology2<R> = HomologyBase<isize2, R>;
pub type Homology3<R> = HomologyBase<isize3, R>;

#[derive(Debug, Clone)]
pub struct HomologySummand<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    rank: usize,
    tors: Vec<R>,
    gens: Option<Vec<SpVec<R>>>,
    trans: Option<SpMat<R>>
}

impl<R> HomologySummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn new(rank: usize, tors: Vec<R>, gens: Option<Vec<SpVec<R>>>, trans: Option<SpMat<R>>) -> Self { 
        Self { rank, tors, gens, trans }
    }

    fn new_gens(rank: usize, tors: Vec<R>, gens: Vec<SpVec<R>>, trans: SpMat<R>) -> Self { 
        Self::new(rank, tors, Some(gens), Some(trans))
    }

    fn new_no_gens(rank: usize, tors: Vec<R>) -> Self { 
        Self::new(rank, tors, None, None)
    }

    fn zero() -> Self { 
        Self::new_no_gens(0, vec![])
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

    pub fn gens(&self) -> &Vec<SpVec<R>> {
        self.gens.as_ref().expect("not computed with gens.")
    }

    pub fn gen(&self, i: usize) -> &SpVec<R> {
        &self.gens()[i]
    }

    pub fn vectorize(&self, z: &SpVec<R>) -> SpVec<R> { 
        let p = self.trans.as_ref().expect("not computed with gens.");
        assert_eq!(p.shape().1, z.dim());
        p * z
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
            return f.write_str("0")
        }
    
        let mut res = vec![];
        let symbol = R::set_symbol();
    
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
    with_gens: bool,
    complex: ChainComplexBase<I, R>,
    cache: HashMap<I, OnceCell<HomologySummand<R>>>,
    zero_smd: HomologySummand<R>
}

impl<I, R> HomologyBase<I, R>
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(complex: ChainComplexBase<I, R>, with_gens: bool) -> Self { 
        let cache = HashMap::from_iter(
            complex.support().map(|i| (i, OnceCell::new()) )
        );
        let zero_smd = HomologySummand::zero();
        Self { with_gens, complex, cache, zero_smd }
    }

    pub fn complex(&self) -> &ChainComplexBase<I, R> {
        &self.complex
    }

    pub fn get(&self, i: I) -> &HomologySummand<R> {
        if self.complex.is_supported(i) { 
            self.cache[&i].get_or_init(|| 
                self.calc(i)
            )
        } else { 
            &self.zero_smd
        }
    }

    fn calc(&self, i: I) -> HomologySummand<R> {
        let i0 = i - self.complex.d_deg();
        let d0 = self.complex.d_matrix(i0);
        let d1 = self.complex.d_matrix(i);

        if self.with_gens { 
            let (res, p, q) = HomologyCalc::calculate_with_trans(&d0, &d1);
            let l = q.shape().1;
            let gens = (0..l).map(|j| q.col_vec(j)).collect();
            HomologySummand::new_gens(res.rank(), res.tors().to_owned(), gens, p)
        } else { 
            let res = HomologyCalc::calculate(&d0, &d1);
            HomologySummand::new_no_gens(res.rank(), res.tors().to_owned())
        }
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
        self.complex.support()
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
        let c = ChainComplex::<i32>::from_mats(-1,
            vec![
                SpMat::from_vec((0, 1), vec![]),
                SpMat::from_vec((1, 0), vec![])
            ]
        );        
        let h = c.homology();
        
        assert_eq!(h[0].rank(), 1);
        assert!( h[0].is_free());
        assert!(!h[0].is_zero());
    }

    #[test]
    fn torsion() { 
        let c = ChainComplex::<i32>::from_mats(-1,
            vec![
                SpMat::from_vec((0, 1), vec![]),
                SpMat::from_vec((1, 1), vec![2]),
                SpMat::from_vec((1, 0), vec![]),
            ]
        );
        let h = c.homology();
        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[0].tors(), &vec![2]);
        assert!(!h[0].is_free());

    }

    #[test]
    fn d3() {
        let c = Samples::<i32>::d3();
        let h = c.homology();

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
        let h = c.homology();

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
        let h = c.homology();

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
        let h = c.homology();

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
        let h = c.homology_with_gens();
        let c = &h.complex;

        let h2 = &h[2];
        assert_eq!(h2.gens().len(), 1);

        let z = h2.gen(0);
        assert!(!z.is_zero());
        assert!(c.differentiate(2, z).is_zero());
        assert_eq!(h2.vectorize(z), SpVec::from(vec![1]));
    }

    #[test]
    fn t2_gens() {
        let c = Samples::<i32>::t2();
        let h = c.homology_with_gens();
        let c = &h.complex;

        let h2 = &h[2];
        assert_eq!(h2.gens().len(), 1);

        let z = h2.gen(0);
        assert!(!z.is_zero());
        assert!(c.differentiate(2, z).is_zero());
        assert_eq!(h2.vectorize(z), SpVec::from(vec![1]));
        assert_eq!(h2.gens().len(), 1);

        let h1 = &h[1];
        assert_eq!(h1.gens().len(), 2);

        let a = h1.gen(0);
        let b = h1.gen(1);

        assert!(!a.is_zero());
        assert!(!b.is_zero());
        assert!(c.differentiate(1, a).is_zero());
        assert!(c.differentiate(1, b).is_zero());
        assert_eq!(h1.vectorize(a), SpVec::from(vec![1, 0]));
        assert_eq!(h1.vectorize(b), SpVec::from(vec![0, 1]));
    }

    #[test]
    fn rp2_gens() {
        let c = Samples::<i32>::rp2();
        let h = c.homology_with_gens();
        let c = &h.complex;

        let h1 = &h[1];
        assert_eq!(h1.gens().len(), 1);

        let z = h1.gen(0);

        assert!(!z.is_zero());
        assert!(c.differentiate(1, z).is_zero());
        assert_eq!(h1.vectorize(z), SpVec::from(vec![1]));
    }
}