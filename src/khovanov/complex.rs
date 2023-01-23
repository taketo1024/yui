use std::collections::HashMap;
use std::ops::{RangeInclusive, Index};
use std::vec::IntoIter;

use crate::math::homology::base::{GradedRModStr, RModGrid};
use crate::math::homology::free::FreeRModStr;
use crate::math::traits::{Ring, RingOps};
use crate::math::homology::complex::ChainComplex;
use crate::links::Link;
use crate::utils::misc::{Idx2, Idx2Range};
use super::algebra::{KhAlgStr, KhEnhState, KhChain};
use super::cube::KhCube;

pub type KhComplexSummand<R> = FreeRModStr<KhEnhState, R>;
pub struct KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    link: Link,
    cube: KhCube<R>,
    grid: RModGrid<KhComplexSummand<R>, RangeInclusive<isize>>,
    reduced: bool
}

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(link: Link, h: R, t: R, reduced: bool) -> Self { 
        let cube = KhCube::new_ht(&link, h, t);

        let i0 = Self::calc_deg_shift(&link, reduced).0;
        let range = cube.h_range().shift(i0);

        let grid = RModGrid::new(range, |i| {
            let i = i - i0;
            let gens = if reduced {
                let e = link.first_edge().unwrap();
                cube.reduced_generators(i, e)
            } else { 
                cube.generators(i) 
            };
            let s = FreeRModStr::new(gens);
            Some(s)
        });

        KhComplex { link, cube, grid, reduced }
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
        self.cube.structure()
    }

    pub fn deg_shift(&self) -> (isize, isize) { 
        Self::calc_deg_shift(&self.link, self.reduced)
    }

    pub fn is_reduced(&self) -> bool { 
        self.reduced
    }

    // TODO abstract as FreeChainComplex
    pub fn differentiate_x(&self, x: &KhEnhState) -> Vec<(KhEnhState, R)> {
        self.cube.differentiate(x)
    }

    pub fn differetiate(&self, z: &KhChain<R>) -> KhChain<R> { 
        z.apply(|x| self.differentiate_x(x))
    }

    fn calc_deg_shift(l: &Link, reduced: bool) -> (isize, isize) {
        let (n_pos, n_neg) = l.signed_crossing_nums();
        let (n_pos, n_neg) = (n_pos as isize, n_neg as isize);
        let h = -n_neg;
        let q = n_pos - 2 * n_neg + (if reduced { 1 } else { 0 });
        (h, q)
    }
}

impl<R> Index<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    type Output = KhComplexSummand<R>;
    
    fn index(&self, index: isize) -> &Self::Output {
        &self.grid[index]
    }
}

impl<R> GradedRModStr for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    type R = R;
    type Index = isize;
    type IndexRange = RangeInclusive<isize>;

    fn in_range(&self, k: Self::Index) -> bool {
        self.grid.in_range(k)
    }

    fn range(&self) -> Self::IndexRange {
        self.grid.range()
    }
}

impl<R> ChainComplex for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    fn d_degree(&self) -> Self::Index {
        1
    }

    fn d_matrix(&self, k: Self::Index) -> sprs::CsMat<Self::R> {
        let c1 = &self.grid[k];
        let c2 = &self.grid[k + 1];
        c1.make_matrix(c2, |x| self.cube.differentiate(x))
    }
}

pub struct KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    cube: KhCube<R>,
    grid: RModGrid<KhComplexSummand<R>, Idx2Range>
}

impl<R> KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(l: Link, reduced: bool) -> Self { 
        use grouping_by::GroupingBy;

        let cube = KhCube::new(&l);
        let (i0, j0) = KhComplex::calc_deg_shift(&l, reduced);
        let h_range = cube.h_range().shift(i0);
        let q_range = cube.q_range().shift(j0);
        
        let range = Idx2::iterate(
            Idx2(*h_range.start(), *q_range.start()), 
            Idx2(*h_range.end(),   *q_range.end()), 
            (1, 2)
        );

        let mut gens: HashMap<_, _> = h_range.clone().map(|i| { 
            let gens = if reduced {
                let e = l.first_edge().unwrap();
                cube.reduced_generators(i - i0, e)
            } else { 
                cube.generators(i - i0)
            };

            let set = gens.into_iter().grouping_by(|x| 
                x.q_deg() + j0
            );
            (i, set)
        }).collect();

        let grid = RModGrid::new(range, |idx| {
            let (i, j) = idx.as_tuple();
            let set = gens.get_mut(&i).unwrap();
            if let Some(g) = set.remove(&j) {
                let s = FreeRModStr::new(g);
                Some(s)
            } else { 
                None
            }
        });

        Self { cube, grid }
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
        &self.grid[index]
    }
}

impl<R> GradedRModStr for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;
    type Index = Idx2;
    type IndexRange = IntoIter<Idx2>;

    fn in_range(&self, k: Self::Index) -> bool {
        self.grid.in_range(k)
    }

    fn range(&self) -> Self::IndexRange {
        self.grid.range()
    }
}

impl<R> ChainComplex for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn d_degree(&self) -> Self::Index {
        Idx2(1, 0)
    }

    fn d_matrix(&self, idx: Self::Index) -> sprs::CsMat<Self::R> {
        let c1 = &self[idx];
        let c2 = &self[idx + self.d_degree()];
        c1.make_matrix(c2, |x| self.cube.differentiate(x))
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
    use crate::links::Link;
    use crate::math::homology::base::GradedRModStr;
    use crate::math::homology::complex::ChainComplexValidation;
    use super::KhComplex;

    #[test]
    fn kh_empty() {
        let l = Link::empty();
        let c = KhComplex::<i32>::unreduced(l);

        assert_eq!(c.range(), 0..=0);
        assert_eq!(c.deg_shift(), (0, 0));

        c.check_d_all();
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let c = KhComplex::<i32>::unreduced(l);

        assert_eq!(c.range(), 0..=0);
        assert_eq!(c.deg_shift(), (0, 0));
        
        c.check_d_all();
    }

    #[test]
    fn kh_unknot_twist() {
        let l = Link::from(&[[0, 0, 1, 1]]);
        let c = KhComplex::<i32>::unreduced(l);

        assert_eq!(c.range(), 0..=1);
        assert_eq!(c.deg_shift(), (0, 1));
        
        c.check_d_all();
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let c = KhComplex::<i32>::unreduced(l);

        assert_eq!(c.range(), -3..=0);
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

        assert_eq!(c.range(), -2..=2);
        assert_eq!(c.deg_shift(), (-2, -2));

        assert_eq!(c[-2].generators().len(), 8);
        assert_eq!(c[-1].generators().len(), 16);
        assert_eq!(c[ 0].generators().len(), 18);
        assert_eq!(c[ 1].generators().len(), 16);
        assert_eq!(c[ 2].generators().len(), 8);

        c.check_d_all();
    }
}