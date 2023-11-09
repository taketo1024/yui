use std::ops::{RangeInclusive, Index};
use cartesian::cartesian;

use delegate::delegate;
use yui_core::{Ring, RingOps, EucRing, EucRingOps, isize2};
use yui_link::Link;
use yui_homology::{ChainComplexTrait, XChainComplex, XChainComplex2, GridTrait, XChainComplexSummand, ChainComplexDisplay, Grid2, XModStr};
use yui_matrix::sparse::SpMat;

use crate::{KhEnhState, KhChain, KhHomology, KhHomologyBigraded};

pub type KhComplexSummand<R> = XChainComplexSummand<KhEnhState, R>;

pub struct KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    inner: XChainComplex<KhEnhState, R>,
    canon_cycles: Vec<KhChain<R>>,
    reduced: bool,
    deg_shift: (isize, isize)
}

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(link: &Link, h: &R, t: &R, reduced: bool) -> Self { 
        Self::new_v2(link, h, t, reduced)
    }

    pub(crate) fn _new(inner: XChainComplex<KhEnhState, R>, canon_cycles: Vec<KhChain<R>>, reduced: bool, deg_shift: (isize, isize)) -> Self { 
        KhComplex { inner, canon_cycles, reduced, deg_shift }
    }

    pub fn canon_cycles(&self) -> &Vec<KhChain<R>> { 
        &self.canon_cycles
    }

    pub fn canon_cycle(&self, i: usize) -> &KhChain<R> { 
        &self.canon_cycles[i]
    }

    pub fn is_reduced(&self) -> bool { 
        self.reduced
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        let h_min = self.support().min().unwrap_or(0);
        let h_max = self.support().max().unwrap_or(0);
        h_min ..= h_max
    }

    pub fn q_range(&self) -> RangeInclusive<isize> { 
        let q_min = self.support().filter_map(|i| self[i].gens().iter().map(|x| x.q_deg()).min()).min().unwrap_or(0);
        let q_max = self.support().filter_map(|i| self[i].gens().iter().map(|x| x.q_deg()).max()).max().unwrap_or(0);
        let q0 = self.deg_shift.1;
        (q_min + q0) ..= (q_max + q0)
    }

    pub fn inner(&self) -> &XChainComplex<KhEnhState, R> {
        &self.inner
    }

    pub fn into_bigraded(self) -> KhComplexBigraded<R> {
        let reduced = self.reduced;
        let deg_shift = self.deg_shift;

        let h_range = self.h_range();
        let q_range = self.q_range().step_by(2);
        let support = cartesian!(h_range, q_range.clone()).map(|(i, j)| 
            isize2(i, j)
        );

        let summands = Grid2::generate(support, |idx| { 
            let isize2(i, j) = idx;
            let q = j - deg_shift.1;
            let gens = self[i].gens().iter().filter(|x| { 
                x.q_deg() == q
            }).cloned();
            XModStr::free(gens)
        });

        let inner = XChainComplex2::new(summands, isize2(1, 0), move |idx, x| { 
            let i = idx.0;
            let x = KhChain::from(x.clone());
            let dx = self.d(i, &x);
            dx.into_iter().collect()
        });

        KhComplexBigraded { inner, reduced }
    }

    pub fn deg_shift_for(l: &Link, reduced: bool) -> (isize, isize) {
        let (n_pos, n_neg) = l.signed_crossing_nums();
        let (n_pos, n_neg) = (n_pos as isize, n_neg as isize);
        let h = -n_neg;
        let q = n_pos - 2 * n_neg + (if reduced { 1 } else { 0 });
        (h, q)
    }

    delegate! { 
        to self.inner { 
            pub fn d(&self, i: isize, z: &KhChain<R>) -> KhChain<R>;
        }
    }
}

impl<R> Index<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = KhComplexSummand<R>;

    delegate! { 
        to self.inner {
            fn index(&self, index: isize) -> &Self::Output;
        }
    }
}

impl<R> GridTrait<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Itr = std::vec::IntoIter<isize>;
    type E = KhComplexSummand<R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: isize) -> bool;
            fn get(&self, i: isize) -> &Self::E;
        }
    }
}

impl<R> ChainComplexTrait<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    delegate! { 
        to self.inner { 
            fn d_deg(&self) -> isize;
            fn d_matrix(&self, i: isize) -> SpMat<R>;
        }
    }
}

impl<R> ChainComplexDisplay<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {}

impl<R> KhComplex<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology(&self, with_trans: bool) -> KhHomology<R> {
        let h = self.inner.homology(with_trans);
        KhHomology::_new(h)
    }
}

pub struct KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    inner: XChainComplex2<KhEnhState, R>,
    reduced: bool,
}

impl<R> KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(l: Link, reduced: bool) -> Self { 
        Self::new_v2(l, reduced)
    }

    pub fn is_reduced(&self) -> bool { 
        self.reduced
    }

    pub fn inner(&self) -> &XChainComplex2<KhEnhState, R> {
        &self.inner
    }
}

impl<R> Index<(isize, isize)> for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = KhComplexSummand<R>;

    delegate! { 
        to self.inner { 
            fn index(&self, index: (isize, isize)) -> &Self::Output;
        }
    }
}

impl<R> GridTrait<isize2> for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Itr = std::vec::IntoIter<isize2>;
    type E = KhComplexSummand<R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Itr;
            fn is_supported(&self, i: isize2) -> bool;
            fn get(&self, i: isize2) -> &Self::E;
        }
    }
}

impl<R> ChainComplexTrait<isize2> for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    delegate! { 
        to self.inner { 
            fn d_deg(&self) -> isize2;
            fn d_matrix(&self, i: isize2) -> SpMat<Self::R>;
        }
    }
}

impl<R> ChainComplexDisplay<isize2> for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {}

impl<R> KhComplexBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn homology(&self, with_trans: bool) -> KhHomologyBigraded<R> {
        let h = self.inner.homology(with_trans);
        KhHomologyBigraded::_new(h)
    }
}