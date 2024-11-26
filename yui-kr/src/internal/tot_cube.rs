use std::ops::RangeInclusive;
use std::sync::Arc;

use itertools::Itertools;
use log::info;
use num_traits::One;
#[cfg(feature = "multithread")]
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use yui::util::sync::SyncCounter;
use yui::{EucRing, EucRingOps};
use yui::bitseq::BitSeq;

use super::base::{KRGen, KRPoly, KRMono, KRChain, KRPolyChain, sign_between, combine, decombine};
use super::data::KRCubeData;
use super::hor_homol::KRHorHomol;

pub struct KRTotCube<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    data: Arc<KRCubeData<R>>,
    q: isize,
    verts: Vec<KRHorHomol<R>> // serialized for fast access
} 

impl<R> KRTotCube<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(data: Arc<KRCubeData<R>>, q: isize) -> Self {
        let n = data.dim() as isize;
        Self::new_restr(data, q, (0..=n, 0..=n))
    }

    pub fn new_restr(data: Arc<KRCubeData<R>>, q: isize, range: (RangeInclusive<isize>, RangeInclusive<isize>)) -> Self {
        let n = data.dim();
        let (h_range, v_range) = range;

        let total = if log::log_enabled!(log::Level::Info) { 
            BitSeq::generate(n).filter(|v| 
                v_range.contains(&(v.weight() as isize))
            ).count()
        } else { 
            0
        };
        let counter = SyncCounter::new();

        info!("setup cube.. (total: {total})");

        let vs = BitSeq::generate(n);

        cfg_if::cfg_if! { 
            if #[cfg(feature = "multithread")] { 
                let itr = vs.collect_vec().into_par_iter();
            } else { 
                let itr = vs;
            }
        };

        let verts = itr.map(|v| {
            if v_range.contains(&(v.weight() as isize)) { 
                if log::log_enabled!(log::Level::Info) { 
                    let c = counter.incr();
                    if c % 3000 == 0 { 
                        info!("  {c}/{total}");
                    }
                }
    
                KRHorHomol::new_restr(
                    data.clone(), 
                    q, 
                    v, 
                    h_range.clone()
                )
            } else { 
                KRHorHomol::zero(data.clone(), q, v)
            }
        }).collect();

        Self { data, q, verts }
    }

    pub fn dim(&self) -> usize { 
        self.data.dim()
    }

    pub fn q_deg(&self) -> isize { 
        self.q
    }

    pub fn vert(&self, v_coords: BitSeq) -> &KRHorHomol<R> {
        let i = v_coords.as_u64() as usize;
        &self.verts[i]
    }

    pub fn edge_poly(&self, h_coords: BitSeq, i: usize) -> KRPoly<R> {
        self.data.ver_edge_poly(h_coords, i)
    }

    #[inline(never)] // for profilability
    pub fn d(&self, z: &KRChain<R>) -> KRChain<R> { 
        let z = combine(z.clone());
        let n = self.dim();

        let dz = z.iter().flat_map(|(v, f)| { 
            (0..n).filter(|&i|
                v.1[i].is_zero()
            ).map(move |i| {
                let w = KRGen(v.0, v.1.edit(|b| b.set_1(i)), KRMono::one());
                let e = R::from_sign( sign_between(v.1, w.1) );
                let g = self.edge_poly(v.0, i);
                let h = f * g * e;
                (w, h)
            })
        }).collect::<KRPolyChain<_>>();

        decombine(dz)
    }
}

#[cfg(test)]
mod tests {
    use num_traits::One;
    use yui::Ratio;
    use yui_link::Link;
    use super::*;

    type R = Ratio<i64>;
    type P = KRPoly<R>;

    fn make_cube(l: &Link, q: isize) -> KRTotCube<R> {
        let data = Arc::new( KRCubeData::<R>::new(l) );
        KRTotCube::new(data, q)
    }

    #[test]
    fn edge_poly() { 
        let x = P::variable;
        let l = Link::from_pd_code([[1,4,2,5],[5,2,6,3],[3,6,4,1]]); // trefoil
        let cube = make_cube(&l, 0);

        // p: neg, h: 0, v: 0 -> 1
        let h  = BitSeq::from([0,0,0]);
        let p = cube.edge_poly(h, 0); // 1
        assert_eq!(p, P::one());

        // p: neg, h: 1, v: 0 -> 1
        let h  = BitSeq::from([1,0,0]);
        let p = cube.edge_poly(h, 0); // x_bc
        assert_eq!(p, x(0));

        let l = l.mirror(); // trefoil
        let cube = make_cube(&l, 0);

        // p: pos, h: 0, v: 0 -> 1
        let h  = BitSeq::from([0,0,0]);
        let p = cube.edge_poly(h, 0); // x_bc
        assert_eq!(p, x(0));

        // p: pos, h: 1, v: 0 -> 1
        let h  = BitSeq::from([1,0,0]);
        let p = cube.edge_poly(h, 0); // x_bc
        assert_eq!(p, P::one());
    }
}