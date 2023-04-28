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
        let complex = KhComplex::new_v1(l, h, t, reduced);
        Self::from(&complex)
    }
}

impl<R> KhHomologyBigraded<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new_v1(l: Link, reduced: bool) -> Self {
        let c = KhComplexBigraded::new_v1(l, reduced);
        Self::from(&c)
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
        let c = KhComplex::new_v1(&l, &0, &0, false);

        assert_eq!(c.indices(), 0..=0);

        c.check_d_all();
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let c = KhComplex::new_v1(&l, &0, &0, false);

        assert_eq!(c.indices(), 0..=0);
        
        c.check_d_all();
    }

    #[test]
    fn kh_unknot_twist() {
        let l = Link::from_pd_code([[0, 0, 1, 1]]);
        let c = KhComplex::new_v1(&l, &0, &0, false);

        assert_eq!(c.indices(), 0..=1);
        
        c.check_d_all();
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let c = KhComplex::new_v1(&l, &0, &0, false);

        assert_eq!(c.indices(), -3..=0);

        assert_eq!(c[-3].generators().len(), 8);
        assert_eq!(c[-2].generators().len(), 12);
        assert_eq!(c[-1].generators().len(), 6);
        assert_eq!(c[ 0].generators().len(), 4);

        c.check_d_all();
    }

    #[test]
    fn kh_figure8() {
        let l = Link::figure8();
        let c = KhComplex::new_v1(&l, &0, &0, false);

        assert_eq!(c.indices(), -2..=2);

        assert_eq!(c[-2].generators().len(), 8);
        assert_eq!(c[-1].generators().len(), 16);
        assert_eq!(c[ 0].generators().len(), 18);
        assert_eq!(c[ 1].generators().len(), 16);
        assert_eq!(c[ 2].generators().len(), 8);

        c.check_d_all();
    }
}
