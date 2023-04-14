use std::fmt::Display;
use std::ops::{RangeInclusive, Index};
use std::rc::Rc;
use std::vec::IntoIter;

use yui_core::{Ring, RingOps};
use yui_matrix::sparse::SpMat;
use yui_link::Link;
use yui_homology::{Idx2, Idx2Iter, Grid, ChainComplex, FreeRModStr, FreeChainComplex};

use crate::{KhAlgStr, KhEnhState, KhCube, KhChain};

pub type KhComplexSummand<R> = FreeRModStr<KhEnhState, R>;
pub struct KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    link: Link,
    str: KhAlgStr<R>,
    complex: FreeChainComplex<KhEnhState, R, RangeInclusive<isize>>,
    reduced: bool
}

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(link: Link, h: R, t: R, reduced: bool) -> Self { 
        let cube = KhCube::new_ht(&link, h, t);
        let str = cube.structure().clone();

        let i0 = Self::deg_shift_for(&link, reduced).0;
        let range = cube.h_range().shift(i0);

        let cube0 = Rc::new(cube);
        let cube1 = cube0.clone();

        let complex = FreeChainComplex::new(range, 1, 
            |i| {
                let i = i - i0;
                let gens = if reduced {
                    let e = link.first_edge().unwrap();
                    cube0.reduced_generators(i, e)
                } else { 
                    cube0.generators(i) 
                };
                gens.into_iter().cloned().collect()
            },
            move |x| { 
                cube1.differentiate(x)
            }
        );

        KhComplex { link, str, complex, reduced }
    }

    pub fn unreduced(l: Link) -> Self { 
        Self::new(l, R::zero(), R::zero(), false)
    }

    pub fn unreduced_ht(l: Link, h: R, t: R) -> Self { 
        Self::new(l, h, t, false)
    }

    pub fn reduced(l: Link) -> Self { 
        Self::new(l, R::zero(), R::zero(), true)
    }

    pub fn reduced_ht(l: Link, h: R, t: R) -> Self { 
        Self::new(l, h, t, true)
    }

    pub fn link(&self) -> &Link { 
        &self.link
    }

    pub fn structure(&self) -> &KhAlgStr<R> {
        &self.str
    }

    pub fn deg_shift(&self) -> (isize, isize) { 
        Self::deg_shift_for(&self.link, self.reduced)
    }

    pub fn is_reduced(&self) -> bool { 
        self.reduced
    }

    pub fn differentiate_x(&self, x: &KhEnhState) -> Vec<(KhEnhState, R)> {
        self.complex.differentiate_x(x)
    }

    pub fn differetiate(&self, z: &KhChain<R>) -> KhChain<R> { 
        self.complex.differetiate(z)
    }

    pub fn deg_shift_for(l: &Link, reduced: bool) -> (isize, isize) {
        let (n_pos, n_neg) = l.signed_crossing_nums();
        let (n_pos, n_neg) = (n_pos as isize, n_neg as isize);
        let h = -n_neg;
        let q = n_pos - 2 * n_neg + (if reduced { 1 } else { 0 });
        (h, q)
    }
}

impl<R> Display for KhComplex<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.complex.fmt(f)
    }
}

impl<R> Index<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    type Output = KhComplexSummand<R>;
    
    fn index(&self, index: isize) -> &Self::Output {
        &self.complex[index]
    }
}

impl<R> Grid for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    type Idx = isize;
    type IdxIter = RangeInclusive<isize>;
    type Output = KhComplexSummand<R>;

    fn contains_idx(&self, k: Self::Idx) -> bool {
        self.complex.contains_idx(k)
    }

    fn indices(&self) -> Self::IdxIter {
        self.complex.indices()
    }

    fn get(&self, i: Self::Idx) -> Option<&Self::Output> {
        self.complex.get(i)
    }
}

impl<R> ChainComplex for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    fn d_degree(&self) -> Self::Idx {
        self.complex.d_degree()
    }

    fn d_matrix(&self, k: Self::Idx) -> SpMat<Self::R> {
        self.complex.d_matrix(k)
    }
}

pub struct KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    complex: FreeChainComplex<KhEnhState, R, Idx2Iter>
}

impl<R> KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(l: Link, reduced: bool) -> Self { 
        let cube = KhCube::new(&l);
        let (i0, j0) = KhComplex::deg_shift_for(&l, reduced);
        let h_range = cube.h_range().shift(i0);
        let q_range = cube.q_range().shift(j0);
        
        let start = Idx2(*h_range.start(), *q_range.start());
        let end   = Idx2(*h_range.end(),   *q_range.end());
        let range = start.iter_rect(end, (1, 2));

        let cube0 = Rc::new(cube);
        let cube1 = cube0.clone();

        let complex = FreeChainComplex::new(range, Idx2(1, 0), 
            |idx| {
                let (i, j) = idx.as_tuple();
                let i = i - i0;
                let j = j - j0;

                let gens = if reduced {
                    let e = l.first_edge().unwrap();
                    cube0.reduced_generators(i, e)
                } else { 
                    cube0.generators(i)
                };
                gens.into_iter().filter(|x| { 
                    x.q_deg() == j
                }).cloned().collect()
            },
            move |x| { 
                cube1.differentiate(x)
            }
        );

        Self { complex }
    }

    pub fn unreduced(l: Link) -> Self { 
        Self::new(l, false)
    }

    pub fn reduced(l: Link) -> Self { 
        Self::new(l, true)
    }
}

impl<R> Index<Idx2> for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = KhComplexSummand<R>;

    fn index(&self, index: Idx2) -> &Self::Output {
        &self.complex[index]
    }
}

impl<R> Grid for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Idx = Idx2;
    type IdxIter = IntoIter<Idx2>;
    type Output = KhComplexSummand<R>;

    fn contains_idx(&self, k: Self::Idx) -> bool {
        self.complex.contains_idx(k)
    }

    fn indices(&self) -> Self::IdxIter {
        self.complex.indices()
    }

    fn get(&self, i: Self::Idx) -> Option<&Self::Output> {
        self.complex.get(i)
    }
}

impl<R> ChainComplex for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn d_degree(&self) -> Self::Idx {
        self.complex.d_degree()
    }

    fn d_matrix(&self, idx: Self::Idx) -> SpMat<Self::R> {
        self.complex.d_matrix(idx)
    }
}

trait Shift { 
    fn shift(self, i: isize) -> Self;
}

impl Shift for RangeInclusive<isize> { 
    fn shift(self, i: isize) -> Self {
        RangeInclusive::new(
            self.start() + i, 
            self.end() + i
        )
    }
}

#[cfg(test)]
mod tests {
    use yui_link::Link;
    use yui_homology::Grid;
    use yui_homology::test::ChainComplexValidation;
    use super::KhComplex;

    #[test]
    fn kh_empty() {
        let l = Link::empty();
        let c = KhComplex::<i32>::unreduced(l);

        assert_eq!(c.indices(), 0..=0);
        assert_eq!(c.deg_shift(), (0, 0));

        c.check_d_all();
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let c = KhComplex::<i32>::unreduced(l);

        assert_eq!(c.indices(), 0..=0);
        assert_eq!(c.deg_shift(), (0, 0));
        
        c.check_d_all();
    }

    #[test]
    fn kh_unknot_twist() {
        let l = Link::from(&[[0, 0, 1, 1]]);
        let c = KhComplex::<i32>::unreduced(l);

        assert_eq!(c.indices(), 0..=1);
        assert_eq!(c.deg_shift(), (0, 1));
        
        c.check_d_all();
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let c = KhComplex::<i32>::unreduced(l);

        assert_eq!(c.indices(), -3..=0);
        assert_eq!(c.deg_shift(), (-3, -6));

        assert_eq!(c[-3].generators().len(), 8);
        assert_eq!(c[-2].generators().len(), 12);
        assert_eq!(c[-1].generators().len(), 6);
        assert_eq!(c[ 0].generators().len(), 4);

        c.check_d_all();
    }

    #[test]
    fn kh_figure8() {
        let l = Link::figure8();
        let c = KhComplex::<i32>::unreduced(l);

        assert_eq!(c.indices(), -2..=2);
        assert_eq!(c.deg_shift(), (-2, -2));

        assert_eq!(c[-2].generators().len(), 8);
        assert_eq!(c[-1].generators().len(), 16);
        assert_eq!(c[ 0].generators().len(), 18);
        assert_eq!(c[ 1].generators().len(), 16);
        assert_eq!(c[ 2].generators().len(), 8);

        c.check_d_all();
    }
}