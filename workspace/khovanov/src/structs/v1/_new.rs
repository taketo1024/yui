use yui_core::{Ring, RingOps, EucRing, EucRingOps};
use yui_link::Link;

use super::canon_cycle::CanonCycles;
use super::cube::KhCube;

use crate::{KhComplex, KhChain, KhHomology, KhComplexBigraded, KhHomologyBigraded};

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new_v1(l: &Link, h: &R, t: &R, reduced: bool) -> Self { 
        let deg_shift = Self::deg_shift_for(l, reduced);
        let canon_cycles = if t.is_zero() && l.is_knot() {
            let ori = if reduced { vec![true] } else { vec![true, false] };
            ori.into_iter().map(|o| 
                KhChain::canon_cycle(l, &R::zero(), h, o)
            ).collect()
        } else { 
            vec![]
        };
        let cube = KhCube::new(l, h, t);
        let complex = cube.as_complex(deg_shift.0, reduced);

        KhComplex::_new(complex, canon_cycles, reduced, deg_shift)
    }        
}

impl<R> KhComplexBigraded<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new_v1(l: Link, reduced: bool) -> Self { 
        let c = KhComplex::new_v1(&l, &R::zero(), &R::zero(), reduced);
        c.as_bigraded()
    }
}

impl<R> KhHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new_v1(l: &Link, h: &R, t: &R, reduced: bool) -> Self {
        let c = KhComplex::new_v1(l, h, t, reduced);
        Self::from(c)
    }
}

impl<R> KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new_v1(l: Link, reduced: bool) -> Self {
        let c = KhComplexBigraded::new_v1(l, reduced);
        Self::from(c)
    }
}

#[cfg(test)]
mod tests {
    use yui_homology::v2::ChainComplexTrait;
    use yui_link::Link;
    use super::KhComplex;

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let c = KhComplex::new_v1(&l, &0, &0, false);

        assert_eq!(c.h_range(), -3..=0);

        assert_eq!(c.rank(-3), 8);
        assert_eq!(c.rank(-2), 12);
        assert_eq!(c.rank(-1), 6);
        assert_eq!(c.rank( 0), 4);

        c.check_d_all();

        let h = c.homology();
    
        assert_eq!(h.h_range(), -3..=0);

        assert_eq!(h[-3].rank(), 1);
        assert_eq!(h[-3].is_free(), true);

        assert_eq!(h[-2].rank(), 1);
        assert_eq!(h[-2].tors(), &vec![2]);

        assert_eq!(h[-1].is_zero(), true);

        assert_eq!(h[ 0].rank(), 2);
        assert_eq!(h[ 0].is_free(), true);
    }
}
