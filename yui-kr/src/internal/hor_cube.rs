use std::sync::Arc;

use num_traits::One;
use yui::{Ring, RingOps};
use yui_homology::{ChainComplex, Grid1, Summand};
use yui::bitseq::BitSeq;

use super::base::{KRPoly, KRMono, KRGen, KRChain, KRPolyChain, sign_between, combine, decombine};
use super::data::KRCubeData;

pub struct KRHorCube<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    data: Arc<KRCubeData<R>>,
    q: isize,
    v_coords: BitSeq,
} 

impl<R> KRHorCube<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(data: Arc<KRCubeData<R>>, q: isize, v_coords: BitSeq) -> Self {
        assert_eq!(v_coords.len(), data.dim());
        Self {
            data,
            q,
            v_coords,
        }
    }

    pub fn dim(&self) -> usize { 
        self.data.dim()
    }

    pub fn q_deg(&self) -> isize { 
        self.q
    }

    pub fn data(&self) -> &KRCubeData<R> { 
        self.data.as_ref()
    }

    fn mon_deg(&self, h_coords: BitSeq) -> isize { 
        let i0 = self.data.root_grad().0;
        let i  = self.data.grad_at(h_coords, self.v_coords).0; // <= i0
        let h = h_coords.weight() as isize;

        self.q + h + (i0 - i) / 2
    }

    pub fn vert_gens(&self, h_coords: BitSeq) -> impl Iterator<Item = KRGen> + '_ {
        let deg = self.mon_deg(h_coords);
        
        let gens = if deg >= 0 {
            let n = self.dim();
            KRMono::generate(n, deg as usize)
        } else { 
            // a trick to make empty iter
            let mut itr = KRMono::generate(1, 1);
            itr.next();
            itr
        };

        gens.map(move |x| 
            KRGen(h_coords, self.v_coords, x)
        )
    }

    pub fn gens(&self, i: usize) -> impl Iterator<Item = KRGen> + '_ {
        self.data.verts_of_weight(i).iter().flat_map(|&v| 
            self.vert_gens(v)
        )
    }

    pub fn edge_poly(&self, i: usize) -> KRPoly<R> {
        self.data.hor_edge_poly(self.v_coords, i)
    }

    pub fn d(&self, z: &KRChain<R>) -> KRChain<R> { 
        let z = combine(z.clone());
        let n = self.dim();

        let dz = z.iter().flat_map(|(v, f)| { 
            (0..n).filter(|&i|
                v.0[i].is_zero()
            ).map(move |i| {
                let w = KRGen(v.0.edit(|b| b.set_1(i)), v.1, KRMono::one());
                let e = R::from_sign( sign_between(v.0, w.0) );
                let g = self.edge_poly(i);
                let h = f * g * e;
                (w, h)
            })
        }).collect::<KRPolyChain<_>>();

        decombine(dz)
    }

    #[allow(unused)]
    fn d_x(&self, e: &KRGen) -> KRChain<R> { 
        self.d(&KRChain::from(e.clone()))
    }

    pub fn into_complex(self) -> ChainComplex<KRGen, R> {
        let n = self.dim() as isize;
        let summands = Grid1::generate(0..=n, |i| { 
            let gens = self.gens(i as usize);
            Summand::from_raw_gens(gens)
        });
        ChainComplex::new(summands, 1, move |_, z| {
            self.d(z)
        })
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use num_traits::{One, Zero};
    use yui_homology::{SummandTrait, ChainComplexTrait};
    use yui::Ratio;
    use yui_link::Link;
    use yui::util::macros::hashmap;
    use super::*;

    type R = Ratio<i64>;
    type P = KRPoly<R>;

    fn make_cube(l: &Link, q: isize, v: BitSeq) -> KRHorCube<R> {
        let data = KRCubeData::<R>::new_no_excl(l);
        let rc = Arc::new(data);
        
        KRHorCube::new(rc, q, v)
    }

    #[test]
    fn mon_deg() { 
        let l = Link::trefoil();

        let v = BitSeq::from([0,0,0]);
        let h = BitSeq::from([0,0,0]);
        let q = 0;

        let cube = make_cube(&l, q, v);
        let deg = cube.mon_deg(h);
        assert_eq!(deg, 0);

        let v = BitSeq::from([0,0,0]);
        let h = BitSeq::from([0,0,0]);
        let q = 1;

        let cube = make_cube(&l, q, v);
        let deg = cube.mon_deg(h);
        assert_eq!(deg, 1);

        let v = BitSeq::from([1,0,0]);
        let h = BitSeq::from([0,0,0]);
        let q = 0;

        let cube = make_cube(&l, q, v);
        let deg = cube.mon_deg(h);
        assert_eq!(deg, 0);

        let v = BitSeq::from([0,0,0]);
        let h = BitSeq::from([1,0,0]);
        let q = 0;

        let cube = make_cube(&l, q, v);
        let deg = cube.mon_deg(h);
        assert_eq!(deg, 1);

        let v = BitSeq::from([1,0,0]);
        let h = BitSeq::from([1,0,0]);
        let q = 0;

        let cube = make_cube(&l, q, v);
        let deg = cube.mon_deg(h);
        assert_eq!(deg, 2);
    }

    #[test]
    fn edge_poly() { 
        let x = P::variable;

        // p: neg, v: 0, h: 0 -> 1
        let l = Link::from_pd_code([[1,4,2,5],[5,2,6,3],[3,6,4,1]]); // trefoil
        let v = BitSeq::from([0,0,0]);
        let q = 0;
        let cube = make_cube(&l, q, v);

        let p = cube.edge_poly(0); // x_ac
        assert_eq!(p, -&x(1) - &x(2));

        // p: neg, v: 1, h: 0 -> 1
        let v = BitSeq::from([1,0,0]);
        let q = 0;
        let cube = make_cube(&l, q, v);

        let p = cube.edge_poly(0); // x_ac * x_bc
        assert_eq!(p, (-&x(1) - &x(2)) * &x(0));
    }

    #[test]
    fn differentiate() { 
        let one = KRMono::one();
        let x = |i| KRMono::from((i, 1));

        let l = Link::from_pd_code([[1,4,2,5],[5,2,6,3],[3,6,4,1]]); // trefoil
        let v = BitSeq::from([0,0,0]);
        let q = 0;
        let cube = make_cube(&l, q, v);

        let h = BitSeq::from([0,0,0]);
        let z = KRGen(h, v, one.clone());
        let ys = cube.d_x(&z);

        assert_eq!(ys, hashmap!{
            KRGen(BitSeq::from([1,0,0]), v, x(1)) => -R::one(),
            KRGen(BitSeq::from([1,0,0]), v, x(2)) => -R::one(),
            KRGen(BitSeq::from([0,1,0]), v, x(1)) =>  R::one(),
            KRGen(BitSeq::from([0,0,1]), v, x(2)) =>  R::one()
        }.into());

        let h = BitSeq::from([0,1,0]);
        let z = KRGen(h, v, one.clone());
        let ys = cube.d_x(&z);

        assert_eq!(ys, hashmap! {
            KRGen(BitSeq::from([1,1,0]), v, x(1)) => -R::one(),
            KRGen(BitSeq::from([1,1,0]), v, x(2)) => -R::one(),
            KRGen(BitSeq::from([0,1,1]), v, x(2)) => -R::one()
        }.into());

        let h = BitSeq::from([1,1,1]);
        let z = KRGen(h, v, one.clone());
        let ys = cube.d_x(&z);

        assert!(ys.is_zero());
    }

    #[test]
    fn vert_gens() {
        let l = Link::trefoil();
        let v = BitSeq::from([0,0,0]);
        let q = 0;

        let cube = make_cube(&l, q, v);
        let h = BitSeq::from([0, 0, 0]);
        let gens = cube.vert_gens(h).collect_vec();

        assert_eq!(gens.len(), 1);
        assert!(gens.iter().all(|x| x.0 == h));

        let h = BitSeq::from([1, 0, 0]);
        let gens = cube.vert_gens(h).collect_vec();

        assert_eq!(gens.len(), 3);
        assert!(gens.iter().all(|x| x.0 == h));
    }

    #[test]
    fn vert_gens_empty() {
        let l = Link::trefoil();
        let v = BitSeq::from([0,0,0]);
        let q = -4;

        let cube = make_cube(&l, q, v);
        let h = BitSeq::from([0, 0, 0]);
        let gens = cube.vert_gens(h).collect_vec();

        assert_eq!(gens.len(), 0);
    }

    #[test]
    fn gens() {
        let l = Link::trefoil();
        let v = BitSeq::from([0,0,0]);
        let q = 0;
        let cube = make_cube(&l, q, v);

        let gens = cube.gens(0).collect_vec();
        assert_eq!(gens.len(), 1);

        let gens = cube.gens(1).collect_vec();
        assert_eq!(gens.len(), 9);
    }

    #[test]
    fn as_complex() { 
        let l = Link::trefoil();
        let v = BitSeq::from([1,0,0]);
        let q = 0;
        let cube = make_cube(&l, q, v);
        let c = cube.into_complex();

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 12);
        assert_eq!(c[2].rank(), 26);
        assert_eq!(c[3].rank(), 15);

        c.check_d_all();

        let v = BitSeq::from([1,0,0]);
        let q = 1;
        let cube = make_cube(&l.mirror(), q, v);
        let c = cube.into_complex();
        
        assert_eq!(c[0].rank(), 6);
        assert_eq!(c[1].rank(), 40);
        assert_eq!(c[2].rank(), 70);
        assert_eq!(c[3].rank(), 36);
        
        c.check_d_all();
    }

    #[test]
    fn homology() { 
        let l = Link::trefoil();
        let v = BitSeq::from([1,0,0]);
        let q = 0;
        let cube = make_cube(&l, q, v);

        let c = cube.into_complex();
        let h = c.homology();

        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[3].rank(), 1);

        let l = Link::trefoil().mirror();
        let v = BitSeq::from([1,0,0]);
        let q = 1;
        let cube = make_cube(&l, q, v);
        let c = cube.into_complex();
        let h = c.homology();
        
        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 2);
    }
}