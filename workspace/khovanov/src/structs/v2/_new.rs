use itertools::Itertools;
use num_traits::Zero;
use yui_core::{Ring, RingOps, EucRing, EucRingOps};
use yui_homology::ChainComplexTrait;
use yui_link::Link;

use crate::{KhComplex, KhComplexBigraded, KhHomology, KhHomologyBigraded};

use super::builder::TngComplexBuilder;

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new_v2(l: &Link, h: &R, t: &R, reduced: bool) -> Self { 
        let mut b = TngComplexBuilder::new(l, h, t, reduced);

        if t.is_zero() && l.is_knot() {
            b.make_canon_cycles();
        }
        
        b.process();

        let canon_cycles = b.canon_cycles().iter().map(|z| 
            z.eval(h, t)
        ).collect_vec();
        let complex = b.into_complex().eval(h, t);

        assert!(canon_cycles.iter().all(|z| 
            complex.d(0, z).is_zero()
        ));

        let deg_shift = Self::deg_shift_for(l, reduced);
        Self::_new(complex, canon_cycles, reduced, deg_shift)
    }        
}

impl<R> KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new_v2(l: Link, reduced: bool) -> Self { 
        let c = KhComplex::new_v2(&l, &R::zero(), &R::zero(), reduced);
        c.into_bigraded()
    }
}

impl<R> KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new_v2(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        KhComplex::new_v2(l, h, t, reduced).homology(false)
    }
}

impl<R> KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new_v2(l: Link, reduced: bool) -> Self {
        KhComplexBigraded::new_v2(l, reduced).homology(false)
    }
}

#[cfg(test)]
mod tests {
    use yui_homology::{RModStr, DisplaySeq};

    use super::*;
 
    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let c = KhComplex::new_v2(&l, &0, &0, false);
        let h: KhHomology<_> = c.homology(false);

        assert_eq!(h.h_range(), -3..=0);

        assert_eq!(h[-3].rank(), 1);
        assert_eq!(h[-3].is_free(), true);

        assert_eq!(h[-2].rank(), 1);
        assert_eq!(h[-2].tors(), &vec![2]);

        assert_eq!(h[-1].is_zero(), true);

        assert_eq!(h[ 0].rank(), 2);
        assert_eq!(h[ 0].is_free(), true);

        h.print_seq();
    }
}