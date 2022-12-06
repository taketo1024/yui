use crate::math::homology::chain_complex::ChainComplex;
use crate::math::homology::homology::{HomologyComputable, SimpleHomology};
use crate::math::matrix::CsMatElem;
use crate::math::traits::{Ring, RingOps, EucRing, EucRingOps};
use crate::links::Link;
use super::algebra::KhAlgStr;
use super::cube::{KhEnhState, KhCube};

pub struct KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    cube: KhCube<R>
}

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(l: &Link, str: KhAlgStr<R>) -> Self { 
        let cube = KhCube::new(l, str);
        KhComplex { cube }
    }
}

impl<R> ChainComplex for KhComplex<R>
where R: Ring + CsMatElem, for<'x> &'x R: RingOps<R> { 
    type R = R;
    type Generator = KhEnhState;

    fn hdeg_range(&self) -> std::ops::RangeInclusive<isize> {
        let n = self.cube.dim() as isize;
        0 ..= n
    }

    fn generators(&self, k: isize) -> Vec<&Self::Generator> {
        let k = k as usize;
        self.cube.generators(k)
    }

    fn differentiate(&self, _k: isize, x:&Self::Generator) -> Vec<(Self::Generator, Self::R)> {
        self.cube.differentiate(x)
    }
}

impl<R> HomologyComputable for KhComplex<R>
where R: EucRing + CsMatElem, for<'x> &'x R: EucRingOps<R> { 
    type Homology = SimpleHomology<R>;
    fn homology(&self) -> Self::Homology {
        SimpleHomology::from(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::{links::Link, khovanov::algebra::KhAlgStr, math::homology::homology::{HomologyComputable, Homology}};
    use super::KhComplex;

    #[test]
    fn kh_empty() {
        let l = Link::empty();
        let a = KhAlgStr::new(0, 0);
        let c = KhComplex::new(&l, a);
        let h = c.homology();

        assert_eq!(h.range(), 0..=0);
        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);
    }

    #[test]
    fn kh_unknot() {
        let l = Link::unknot();
        let a = KhAlgStr::new(0, 0);
        let c = KhComplex::new(&l, a);
        let h = c.homology();

        assert_eq!(h.range(), 0..=0);
        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[0].is_free(), true);
    }

    #[test]
    fn kh_trefoil() {
        let l = Link::trefoil();
        let a = KhAlgStr::new(0, 0);
        let c = KhComplex::new(&l, a);
        let h = c.homology();

        assert_eq!(h.range(), 0..=3);
        for i in h.range() { 
            dbg!(&h[i]);
        }
    }
}