use std::ops::RangeInclusive;
use cartesian::cartesian;

use delegate::delegate;
use yui_core::{Ring, RingOps, EucRing, EucRingOps, isize2};
use yui_matrix::sparse::{SpMat, SpVec};
use yui_link::Link;
use yui_homology::{ChainComplexTrait, XChainComplex, XChainComplex2, Graded};

use crate::{KhEnhState, KhChain, KhHomology, KhHomologyBigraded};

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
        let q_min = self.support().filter_map(|i| self.gens(i).map(|x| x.q_deg()).min()).min().unwrap_or(0);
        let q_max = self.support().filter_map(|i| self.gens(i).map(|x| x.q_deg()).max()).max().unwrap_or(0);
        let q0 = self.deg_shift.1;
        (q_min + q0) ..= (q_max + q0)
    }

    pub fn inner(&self) -> &XChainComplex<KhEnhState, R> {
        &self.inner
    }

    pub fn as_bigraded(self) -> KhComplexBigraded<R> {
        let reduced = self.reduced;
        let deg_shift = self.deg_shift;

        let h_range = self.h_range();
        let q_range = self.q_range().step_by(2);
        let support = cartesian!(h_range, q_range.clone()).map(|(i, j)| 
            isize2(i, j)
        );

        let complex = XChainComplex2::new(support, isize2(1, 0), 
            |idx| {
                let isize2(i, j) = idx;
                let q = j - deg_shift.1;

                self.gens(i).filter(|x| { 
                    x.q_deg() == q
                }).cloned().collect()
            },
            |idx, x| { 
                let i = idx.0;
                let x = KhChain::from(x.clone());
                let dx = self.d(i, &x);
                dx.into_iter().collect()
            }
        );

        KhComplexBigraded { inner: complex, reduced }
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
            pub fn gens(&self, i: isize) -> impl Iterator<Item = &KhEnhState>;
            pub fn vectorize(&self, i: isize, z: &KhChain<R>) -> SpVec<R>;
            pub fn d(&self, i: isize, z: &KhChain<R>) -> KhChain<R>;
        }
    }
}

impl<R> Graded<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Itr = std::vec::IntoIter<isize>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Itr;
            fn display(&self, i: isize) -> String;
        }
    }
}

impl<R> ChainComplexTrait<isize> for KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    delegate! { 
        to self.inner { 
            fn rank(&self, i: isize) -> usize;
            fn d_deg(&self) -> isize;
            fn d_matrix(&self, i: isize) -> &SpMat<Self::R>;
        }
    }
}

impl<R> KhComplex<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    pub fn homology(self) -> KhHomology<R> { 
        KhHomology::from(self)
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

impl<R> Graded<isize2> for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Itr = std::vec::IntoIter<isize2>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Itr;
            fn display(&self, i: isize2) -> String;
        }
    }
}

impl<R> ChainComplexTrait<isize2> for KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    delegate! { 
        to self.inner { 
            fn rank(&self, i: isize2) -> usize;
            fn d_deg(&self) -> isize2;
            fn d_matrix(&self, i: isize2) -> &SpMat<Self::R>;
        }
    }
}


impl<R> KhComplexBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    pub fn homology(self) -> KhHomologyBigraded<R> { 
        KhHomologyBigraded::from(self)
    }
}