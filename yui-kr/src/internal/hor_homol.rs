use std::ops::{Index, RangeInclusive};
use std::sync::Arc;

use delegate::delegate;
use num_traits::Zero;
#[cfg(feature = "multithread")]
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use yui::{EucRing, EucRingOps};
use yui_homology::{Grid1, GridTrait, Homology, Summand, SummandTrait};
use yui_matrix::sparse::SpVec;
use yui::bitseq::BitSeq;

// use log::info;
// use yui_homology::DisplaySeq;

use super::base::{KRGen, KRChain, extend_ends_bounded};
use super::data::KRCubeData;
use super::hor_cpx::KRHorComplex;
use super::hor_excl::KRHorExcl;

pub struct KRHorHomol<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    q: isize,
    excl: Arc<KRHorExcl<R>>,
    inner: Homology<KRGen, R>
} 

impl<R> KRHorHomol<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    pub fn new(data: Arc<KRCubeData<R>>, q: isize, v_coords: BitSeq) -> Self { 
        let n = data.dim() as isize;
        Self::new_restr(data, q, v_coords, 0..=n)
    }

    pub fn new_restr(data: Arc<KRCubeData<R>>, q: isize, v_coords: BitSeq, h_range: RangeInclusive<isize>) -> Self { 
        let n = data.dim() as isize;
        let h_range_ext = extend_ends_bounded(h_range.clone(), 1, 0..=n);

        let excl = data.excl(v_coords);
        let complex = KRHorComplex::new_restr(data.clone(), q, v_coords, h_range_ext).reduced();
        
        // info!("H_hor (q: {}, v: {:?})..", q, v_coords);

        let inner = Grid1::generate(h_range, |i|
            complex.homology_at(i)
        );

        // info!("H_hor (q: {}, v: {:?})\n{}", q, v_coords, inner.display_seq("h"));

        Self { q, excl, inner }
    }

    pub fn zero(data: Arc<KRCubeData<R>>, q: isize, v_coords: BitSeq) -> Self { 
        let excl = data.excl(v_coords);
        let inner = Grid1::default();
        Self { q, excl, inner }
    }

    pub fn q_deg(&self) -> isize { 
        self.q
    }

    pub fn v_coords(&self) -> BitSeq { 
        self.excl.v_coords()
    }

    pub fn rank(&self, i: isize) -> usize {
        self.inner[i].rank()
    }

    pub fn gens(&self, i: isize) -> Vec<KRChain<R>> { 
        let r = self.rank(i);
        
        cfg_if::cfg_if! { 
            if #[cfg(feature = "multithread")] { 
                let itr = (0..r).into_par_iter();
            } else { 
                let itr = (0..r).into_iter();
            }
        };

        itr.map(|k| { 
            let v = SpVec::unit(r, k);
            self.as_chain(i, &v)
        }).collect()
    }

    #[inline(never)] // for profilability
    pub fn vectorize(&self, i: isize, z: &KRChain<R>) -> SpVec<R> {
        debug_assert!(z.iter().all(|(x, _)| 
            x.0.weight() as isize == i && 
            x.1 == self.v_coords()
        ));

        if self.rank(i).is_zero() { 
            return SpVec::zero(0)
        }

        let z_exc = self.excl.forward(z);
        let h = &self.inner[i];
        h.vectorize(&z_exc)
    }

    pub fn as_chain(&self, i: isize, v_hml: &SpVec<R>) -> KRChain<R> {
        let h = &self.inner[i];
        let z_exc = h.devectorize(v_hml);
        self.excl.backward(&z_exc, true)
    }
}

impl<R> GridTrait<isize> for KRHorHomol<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Support = std::vec::IntoIter<isize>;
    type Item = Summand<KRGen, R>;

    delegate! { 
        to self.inner { 
            fn support(&self) -> Self::Support;
            fn is_supported(&self, i: isize) -> bool;
            fn get(&self, i: isize) -> &Self::Item;
            fn get_default(&self) -> &Self::Item;
        }
    }
}

impl<R> Index<isize> for KRHorHomol<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = Summand<KRGen, R>;
    
    fn index(&self, index: isize) -> &Self::Output { 
        self.get(index)
    }
}

#[cfg(test)]
mod tests { 
    use yui_link::Link;
    use yui::Ratio;
    use super::*;

    type R = Ratio<i64>;

    fn make_hml(l: &Link, q: isize, v: BitSeq) -> KRHorHomol<R> {
        let data = KRCubeData::<R>::new(l);
        let rc = Arc::new(data);
        KRHorHomol::new(rc, q, v)
    }

    #[test]
    fn rank() { 
        let l = Link::trefoil();
        let v = BitSeq::from([1,0,0]);
        let q = 0;
        let hml = make_hml(&l, q, v);

        assert_eq!(hml[0].rank(), 0);
        assert_eq!(hml[1].rank(), 0);
        assert_eq!(hml[2].rank(), 1);
        assert_eq!(hml[3].rank(), 1);
    }


    #[test]
    fn rank2() { 
        let l = Link::from_pd_code([[1,4,2,5],[5,2,6,3],[3,6,4,1]]); // trefoil
        let v = BitSeq::from([1,0,0]);
        let q = -4;
        let hml = make_hml(&l, q, v);

        assert_eq!(hml[0].rank(), 0);
        assert_eq!(hml[1].rank(), 0);
        assert_eq!(hml[2].rank(), 0);
        assert_eq!(hml[3].rank(), 1);
    }

    #[test]
    fn vectorize() { 
        let l = Link::trefoil().mirror();
        let v = BitSeq::from([1,0,0]);
        let q = 1;
        let hml = make_hml(&l, q, v);

        let zs = hml.gens(2);
        assert_eq!(zs.len(), 2);

        let z0 = &zs[0];
        let z1 = &zs[1];
        let w = z0 * R::from(3) + z1 * R::from(-1);

        assert_eq!(hml.vectorize(2, z0), SpVec::from(vec![R::from(1), R::from(0)]));
        assert_eq!(hml.vectorize(2, z1), SpVec::from(vec![R::from(0), R::from(1)]));
        assert_eq!(hml.vectorize(2, &w), SpVec::from(vec![R::from(3), R::from(-1)]));
    }
}