use std::ops::Index;
use std::cell::OnceCell;
use std::collections::HashMap;
use std::fmt::Display;

use itertools::Itertools;
use yui_core::{EucRing, EucRingOps, Ring, RingOps};

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
    tors: Vec<R>
}

impl<R> HomologySummand<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn new(rank: usize, tors: Vec<R>) -> Self { 
        Self { rank, tors }
    }

    fn zero() -> Self { 
        Self::new(0, vec![])
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
    complex: ChainComplexBase<I, R>,
    cache: HashMap<I, OnceCell<HomologySummand<R>>>,
    zero_smd: HomologySummand<R>
}

impl<I, R> HomologyBase<I, R>
where I: Deg, R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(complex: ChainComplexBase<I, R>) -> Self { 
        let cache = HashMap::from_iter(
            complex.support().map(|i| (i, OnceCell::new()) )
        );
        let zero_smd = HomologySummand::zero();
        Self { complex, cache, zero_smd }
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
        let res = HomologyCalc::calculate(&d0, &d1);
        HomologySummand::new(res.rank(), res.tors().clone())
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
}