use std::collections::{HashSet, HashMap};

use log::trace;
use num_traits::{Zero, One};
use yui::{Ring, RingOps, IndexList};
use yui_homology::utils::make_matrix;
use yui_matrix::sparse::Trans;
use yui::poly::Mono;
use yui::bitseq::BitSeq;

use super::base::{KRMono, KRPoly, KRGen, KRChain, KRPolyChain, sign_between, combine, decombine};
use super::data::KRCubeData;

#[derive(Debug, Clone)]
struct Process<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    dir: usize,
    var: usize,
    deg: usize,
    divisor: KRPoly<R>,
    edge_polys: HashMap<usize, KRPoly<R>>
}

impl<R> Process<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn divisor(&self) -> (&KRPoly<R>, usize, usize) { 
        (&self.divisor, self.var, self.deg)
    }
}

#[derive(Debug, Clone)]
pub struct KRHorExcl<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    v_coords: BitSeq,
    edge_polys: HashMap<usize, KRPoly<R>>,
    exc_dirs: HashSet<usize>,
    fst_exc_vars: HashSet<usize>,
    snd_exc_vars: HashSet<usize>,
    remain_vars: HashSet<usize>,
    process: Vec<Process<R>>,
}

impl<R> KRHorExcl<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn new(v_coords: BitSeq, edge_polys: HashMap<usize, KRPoly<R>>) -> Self { 
        let n = v_coords.len();
        Self { 
            v_coords,
            edge_polys,
            exc_dirs: HashSet::new(), 
            fst_exc_vars: HashSet::new(), 
            snd_exc_vars: HashSet::new(),
            remain_vars: (0..n).collect(),
            process: vec![], 
        }
    }

    pub fn from(data: &KRCubeData<R>, v_coords: BitSeq, level: usize) -> Self { 
        assert!(level <= 2);

        let n = data.dim();
        let edge_polys = (0..n).map(|i| 
            (i, data.hor_edge_poly(v_coords, i))
        ).collect();

        let mut res = Self::new(v_coords, edge_polys);

        trace!("start excl for v: {v_coords}");
        trace!("initial: {:#?}", res.edge_polys);

        for d in 1..=level { 
            res.excl_all(d);
        }

        trace!("excluded: {:?}", res.excl_vars());
        trace!("result: {:#?}", res.edge_polys);

        res
    }

    pub fn v_coords(&self) -> BitSeq { 
        self.v_coords
    }

    fn excl_all(&mut self, deg: usize) {
        while let Some((i, k)) = self.find_excl_var(deg) { 
            self.perform_excl(i, k, deg);
        }
    }

    fn find_excl_var(&self, deg: usize) -> Option<(usize, usize)> { 
        if self.remain_vars.is_empty() { 
            return None
        }
        
        let cands = self.edge_polys.keys().filter_map(|&i| { 
            let p = &self.edge_polys[&i];
            self.find_excl_var_in(p, deg).map(|k| (i, k))
        });

        cands.max_by(|(i1, _), (i2, _)| { 
            let p1 = self.edge_polys[&i1].nterms();
            let p2 = self.edge_polys[&i2].nterms();
            // prefer smaller poly.
            Ord::cmp(&p1, &p2).reverse().then(
                // prefer smaller index.
                Ord::cmp(&i1, &i2).reverse()
            )
        })
    }

    fn find_excl_var_in(&self, p: &KRPoly<R>, deg: usize) -> Option<usize> { 
        assert!(deg > 0);

        // all vars consisting p must be contained in remain_vars. 
        if !self.is_free(p) { 
            return None;
        }

        // collect terms of the form a * (x_k)^deg.
        let cands = p.iter().filter_map(|(x, a)| {
            if x.total_deg() == deg && x.deg().ninds() == 1 { 
                let k = x.deg().min_index().unwrap();
                Some((k, a))
            } else {
                None
            }
        });
        
        // choose best candidate
        cands.max_by(|(k1, a1), (k2, a2)| { 
            // prefer coeff ±1
            Ord::cmp(&a1.is_pm_one(), &a2.is_pm_one()).then( 
                // prefer smaller index
                Ord::cmp(&k1, &k2).reverse() 
            )
        }).map(|(k, _)| k)
    }

    fn perform_excl(&mut self, i: usize, k: usize, deg: usize) {
        assert!(deg == 1 || deg == 2);
        
        // take divisor and edge-polys.
        let p = self.edge_polys.remove(&i).unwrap();
        let edge_polys = std::mem::take(&mut self.edge_polys);

        debug_assert!(p.lead_term_for(k).unwrap().0.deg_for(k) == deg);
        debug_assert!(self.is_free(&p));

        // update edge-polys.
        self.edge_polys = edge_polys.clone().into_iter().map(|(j, f)| {
            let f_rem = rem(f, &p, k);
            (j, f_rem)
        }).collect();

        // update other infos.
        self.exc_dirs.insert(i);
        if deg == 1 { 
            self.fst_exc_vars.insert(k);
        } else { 
            self.snd_exc_vars.insert(k);
        }

        // update indep vars.
        if deg == 1 { 
            // vars except x_k remain indep.
            self.remain_vars.remove(&k);
        } else { 
            // vars appearing in p are no longer indep. 
            for x in p.iter().map(|(x, _)| x) {
                for l in x.deg().iter().map(|(l, _)| l) {
                    self.remain_vars.remove(l);
                }
            }
        }

        trace!("excl step: {}", self.process.len());
        trace!("direction: {i}");
        trace!("chosen: {} in {p}", KRMono::from((k, deg)));
        trace!("edge-poly: {:#?}", self.edge_polys);
        trace!("remain: {:?}", self.remain_vars);
        trace!("fst_exc: {:?}", self.fst_exc_vars);
        trace!("snd_exc: {:?}", self.snd_exc_vars);

        // insert new process.
        let d = Process {
            dir: i,
            var: k,
            deg,
            divisor: p,
            edge_polys,
        };
        
        self.process.push(d);
    }

    fn is_free(&self, p: &KRPoly<R>) -> bool { 
        p.iter().all(|(x, _)| 
            x.deg().iter().all(|(k, _)| 
                self.remain_vars.contains(&k)
            )
        )
    }

    pub fn excl_vars(&self) -> Vec<KRMono> { 
        Iterator::chain(
            self.fst_exc_vars.iter().map(|&k| KRMono::from((k, 1))),
            self.snd_exc_vars.iter().map(|&k| KRMono::from((k, 2)))
        ).collect()
    }

    pub fn should_remain(&self, v: &KRGen) -> bool {
        !self.should_drop(v) && !self.should_reduce(&v.2)
    }

    fn should_drop(&self, v: &KRGen) -> bool { 
        let h = &v.0;
        self.exc_dirs.iter().any(|&i| 
            h[i].is_zero()
        )
    }

    fn should_reduce(&self, x: &KRMono) -> bool { 
        x.deg().iter().any(|(k, &d)| 
            self.fst_exc_vars.contains(k) || 
            (self.snd_exc_vars.contains(k) && d >= 2)
        )
    }

    fn is_reduced(&self, p: &KRPoly<R>) -> bool { 
        p.iter().all(|(x, _)| 
            !self.should_reduce(x)
        )
    }

    fn is_reduced_upto(&self, p: &KRPoly<R>, step: usize) -> bool { 
        (0..step).all(|i| { 
            let (_, k, deg) = self.process[i].divisor();
            if let Some((x, _)) = p.lead_term_for(k) { 
                x.deg_for(k) < deg
            } else {
                true
            }
        })
    }

    fn reduce(&self, f: KRPoly<R>) -> KRPoly<R> {
        if self.is_reduced(&f) { 
            return f
        }
        let l = self.process.len();
        self._reduce_upto(f, l)
    }

    fn reduce_upto(&self, f: KRPoly<R>, step: usize) -> KRPoly<R> {
        if self.is_reduced_upto(&f, step) { 
            return f
        }
        self._reduce_upto(f, step)
    }

    fn _reduce_upto(&self, f: KRPoly<R>, step: usize) -> KRPoly<R> {
        (0..step).fold(f, |f, i| { 
            let (p, k, _) = self.process[i].divisor();
            rem(f, p, k)
        })
    }

    pub fn d(&self, z: &KRChain<R>) -> KRChain<R> { 
        let z = combine(z.clone());

        let dz = z.iter().flat_map(|(v, f)| { 
            self.edge_polys.iter().filter(|(&i, _)|
                v.0[i].is_zero()
            ).map(move |(&i, g)| {
                let w = KRGen(v.0.edit(|b| b.set_1(i)), v.1, KRMono::one());
                let e = R::from_sign( sign_between(v.0, w.0) );
                let h = self.reduce(f * g * e);
                (w, h)
            })
        }).collect::<KRPolyChain<_>>();

        decombine(dz)
    }

    fn d_step(&self, z: &KRPolyChain<R>, step: usize, after_excl: bool) -> KRPolyChain<R> { 
        let proc = &self.process[step];
        z.iter().flat_map(|(v, f)| { 
            proc.edge_polys.iter().filter(|(&i, _)|
                v.0[i].is_zero()
            ).map(move |(&i, g)| {
                let w = KRGen(v.0.edit(|b| b.set_1(i)), v.1, KRMono::one());
                let e = R::from_sign( sign_between(v.0, w.0) );
                let h = if after_excl { 
                    let (p, k, _) = proc.divisor();
                    let g = rem(g.clone(), &p, k);
                    self.reduce_upto(f * g * e, step + 1)
                } else { 
                    self.reduce_upto(f * g * e, step)
                };
                (w, h)
            })
        }).collect()
    }

    pub fn forward(&self, z: &KRChain<R>) -> KRChain<R> {
        // keep only non-vanishing gens. 
        let z = z.filter_gens(|v| 
            !self.should_drop(v)
        );

        // reduce monomials. 
        if z.iter().any(|(v, _)| self.should_reduce(&v.2)) { 
            let z = combine(z);
            let res = z.into_map_coeffs::<KRPoly<R>, _>(|f| self.reduce(f));
            decombine(res)
        } else { 
            z
        }
    }

    fn forward_x(&self, v: &KRGen) -> KRChain<R> {
        self.forward(&KRChain::from(v.clone()))
    }

    pub fn backward(&self, z: &KRChain<R>, is_cycle: bool) -> KRChain<R> {
        if is_cycle { 
            debug_assert_eq!(
                self.d(&z), 
                KRChain::zero()
            );
        }
        
        let z = combine(z.clone());

        let l = self.process.len();
        let res = (0..l).rev().fold(z, |z, step|
            self.backward_step(z, step, is_cycle)
        );

        decombine(res)
    }

    fn backward_step(&self, z: KRPolyChain<R>, step: usize, is_cycle: bool) -> KRPolyChain<R> {
        let d = &self.process[step];
        let i = d.dir;
        let (p, k, _) = d.divisor();

        //    iz < . . . . . . z
        //    |                |
        //  d |                | d'
        //    V                V
        //   diz + id'z < . . d'z
        //   ---------
        //       p
        //
        //  ψ(z) := ((diz + id'z) / p, z).

        let iz = self.send_back(&z, i);
        let diz = self.d_step(&iz, step, false);

        let idz = if is_cycle {
            debug_assert_eq!(
                self.d_step(&z, step, true), 
                KRPolyChain::zero()
            );
            KRPolyChain::zero()
        } else {
            let dz = self.d_step(&z, step, true);
            self.send_back(&dz, i)
        };

        let w = (diz + idz).into_map_coeffs::<KRPoly<R>, _>(|f| {
            div(f, p, k)
        });
        
        w + z
    }

    fn backward_x(&self, w: &KRGen) -> KRChain<R> {
        let z = KRChain::from(w.clone());
        self.backward(&z, false)
    }

    fn send_back(&self, z: &KRPolyChain<R>, dir: usize) -> KRPolyChain<R> {
        let i = dir;
        z.map::<_, KRPoly<R>, _>(|v, f| { 
            debug_assert!(v.0[i].is_one());
            let u = KRGen(v.0.edit(|b| b.set_0(i)), v.1, v.2.clone());
            let e = R::from_sign( sign_between(u.0, v.0) );
            (u, f * e)
        })
    }

    pub fn trans_for(&self, from: &IndexList<KRGen>, to: &IndexList<KRGen>) -> Trans<R> {
        let fwd = make_matrix(from, to, |x| self.forward_x(x));
        let bwd = make_matrix(to, from, |x| self.backward_x(x));
        Trans::new(fwd, bwd)
    }
}

fn div_rem<R>(f: KRPoly<R>, g: &KRPoly<R>, k: usize) -> (KRPoly<R>, KRPoly<R>)
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (e0, a0) = g.lead_term_for(k).unwrap();

    assert!(e0.deg_for(k) > 0);
    assert!(a0.is_unit());

    let a0_inv = a0.inv().unwrap();

    let mut q = KRPoly::zero();
    let mut r = f;
    
    while let Some((e1, a1)) = r.lead_term_for(k) { 
        if !e0.divides(e1) {
            break
        }

        let b = a1 * &a0_inv;
        let x = KRPoly::from((e1 / e0, b));

        r -= &x * g;
        q += x;
    }
    
    (q, r)
}

fn div<R>(f: KRPoly<R>, p: &KRPoly<R>, k: usize) -> KRPoly<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (q, r) = div_rem(f, p, k);
    assert_eq!(r, KRPoly::zero());
    q
}

fn rem<R>(f: KRPoly<R>, p: &KRPoly<R>, k: usize) -> KRPoly<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    div_rem(f, p, k).1
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use itertools::Itertools;
    use yui_link::Link;
    use yui_matrix::MatTrait;
    use yui::poly::MultiVar;
    use yui::Ratio;
    use yui::util::macros::hashmap;
    use yui::bitseq::BitSeq;

    use crate::internal::data::KRCubeData;
    use crate::internal::hor_cube::KRHorCube;

    use super::*;

    type R = Ratio<i64>;
    type P = KRPoly<R>;

    fn gen<const N: usize, const M: usize>(h: [usize; N], v: [usize; N], m: [usize; M]) -> KRGen {
        KRGen(BitSeq::from(h), BitSeq::from(v), MultiVar::from(m))
    }

    fn reduce<R>(excl: &KRHorExcl<R>, gens: &Vec<KRGen>) -> Vec<KRGen>
    where R: Ring, for<'x> &'x R: RingOps<R> { 
        gens.iter().filter(|&v| excl.should_remain(v)).cloned().collect_vec()
    }

    #[test]
    fn init() { 
        let l = Link::trefoil();
        let v = BitSeq::from([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let excl = KRHorExcl::from(&data, v, 0);

        assert!(excl.exc_dirs.is_empty());
        assert!(excl.fst_exc_vars.is_empty());
        assert!(excl.snd_exc_vars.is_empty());
        assert_eq!(excl.remain_vars, [0,1,2].into());
        assert!(excl.process.is_empty());

        let x = P::variable;
        assert_eq!(excl.edge_polys, hashmap!{ 
            0 => -x(1) - x(2),
            1 => x(1),
            2 => x(2) * (x(0) + x(1))
        });
    }

    #[test]
    fn find_excl_var() { 
        let l = Link::trefoil();
        let v = BitSeq::from([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let excl = KRHorExcl::from(&data, v, 0);
        let p = &excl.edge_polys;

        assert_eq!(excl.find_excl_var_in(&p[&0], 1), Some(1));
        assert_eq!(excl.find_excl_var_in(&p[&1], 1), Some(1));
        assert_eq!(excl.find_excl_var_in(&p[&2], 1), None);
        
        assert_eq!(excl.find_excl_var_in(&p[&0], 2), None);
        assert_eq!(excl.find_excl_var_in(&p[&1], 2), None);
        assert_eq!(excl.find_excl_var_in(&p[&2], 1), None);
    }

    #[test]
    fn perform_excl() { 
        let l = Link::trefoil();
        let v = BitSeq::from([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);
        let x = P::variable;

        excl.perform_excl(0, 1, 1); // x1 -> -x2

        assert_eq!(excl.edge_polys, hashmap!{ 
            1 => -x(2),
            2 => x(2) * (x(0) - x(2))
        });

        assert_eq!(excl.exc_dirs, [0].into());
        assert_eq!(excl.fst_exc_vars, [1].into());
        assert_eq!(excl.snd_exc_vars, [].into());
        assert_eq!(excl.remain_vars, [0, 2].into());
        assert_eq!(excl.process.len(), 1);

        let proc = &excl.process[0];

        assert_eq!(proc.dir, 0);
        assert_eq!(proc.var, 1);
        assert_eq!(proc.divisor, -x(1) - x(2));
        assert_eq!(proc.edge_polys, hashmap!{ 
            1 => x(1),
            2 => x(2) * (x(0) + x(1))
        });
    }

    #[test]
    fn perform_excl_2() { 
        let l = Link::trefoil();
        let v = BitSeq::from([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);
        let x = P::variable;

        excl.perform_excl(0, 1, 1); // x1 -> -x2
        excl.perform_excl(1, 2, 1); // x2 -> 0

        assert_eq!(excl.edge_polys, hashmap!{ 
            2 => P::zero()
        });

        assert_eq!(excl.exc_dirs, [0,1].into());
        assert_eq!(excl.fst_exc_vars, [1,2].into());
        assert_eq!(excl.snd_exc_vars, [].into());
        assert_eq!(excl.remain_vars, [0].into());
        assert_eq!(excl.process.len(), 2);

        let proc = &excl.process[1];
        
        assert_eq!(proc.dir, 1);
        assert_eq!(proc.var, 2);
        assert_eq!(proc.divisor, -x(2));
        assert_eq!(proc.edge_polys, hashmap!{ 
            2 => x(2) * (x(0) - x(2))
        });
    }

    #[test]
    fn perform_excl_3() { 
        let l = Link::trefoil();
        let v = BitSeq::from([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);
        let x = P::variable;

        excl.perform_excl(0, 1, 1); // x1 -> -x2

        let p = &excl.edge_polys;
        assert_eq!(excl.find_excl_var_in(&p[&1], 2), None);
        assert_eq!(excl.find_excl_var_in(&p[&2], 2), Some(2));

        excl.perform_excl(2, 2, 2); // x2 * x2 -> x0 * x2

        assert_eq!(excl.edge_polys, hashmap!{ 
            1 => -x(2)
        });

        assert_eq!(excl.exc_dirs, [0,2].into());
        assert_eq!(excl.fst_exc_vars, [1].into());
        assert_eq!(excl.snd_exc_vars, [2].into());
        assert_eq!(excl.remain_vars, [].into()); // no indep vars.
        assert_eq!(excl.process.len(), 2);

        let proc = &excl.process[1];

        assert_eq!(proc.dir, 2);
        assert_eq!(proc.var, 2);
        assert_eq!(proc.divisor, x(2) * (x(0) - x(2)));
        assert_eq!(proc.edge_polys, hashmap!{ 
            1 => -x(2)
        });
    }

    #[test]
    fn should_vanish_before() {
        let l = Link::trefoil();
        let q = 0;
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let excl = KRHorExcl::from(&data, v, 0);
        let cube = KRHorCube::new(Arc::new(data), q, v);

        assert!((0..=3).all(|i| 
            cube.gens(i).all(|v| 
                !excl.should_drop(&v)
            )
        ));
    }

    #[test]
    fn should_vanish() {
        let l = Link::trefoil();
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);

        excl.perform_excl(0, 1, 1);

        assert!( excl.should_drop(&gen([0,1,0], [0,0,1], [1,2,3]))); // h[0] = 0
        assert!(!excl.should_drop(&gen([1,0,0], [0,0,1], [1,2,3])));
    }

    #[test]
    fn should_reduce_before() {
        let l = Link::trefoil();
        let q = 0;
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let excl = KRHorExcl::from(&data, v, 0);
        let cube = KRHorCube::new(Arc::new(data), q, v);

        assert!((0..=3).all(|i| 
            cube.gens(i).all(|v| 
                !excl.should_reduce(&v.2)
            )
        ));
    }

    #[test]
    fn should_reduce() {
        let l = Link::trefoil();
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);

        excl.perform_excl(0, 1, 1); // x1 -> -x2

        assert!(!excl.should_reduce(&MultiVar::from([0,0,0])));
        assert!(!excl.should_reduce(&MultiVar::from([1,0,0])));
        assert!( excl.should_reduce(&MultiVar::from([0,1,0])));
        assert!(!excl.should_reduce(&MultiVar::from([0,0,1])));
    }

    #[test]
    fn should_reduce_2() {
        let l = Link::trefoil();
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);

        excl.perform_excl(0, 1, 1); // x1 -> -x2
        excl.perform_excl(2, 2, 2); // x2 * x2 -> x0 * x2

        assert!(!excl.should_reduce(&MultiVar::from([0,0,0])));
        assert!(!excl.should_reduce(&MultiVar::from([1,0,0])));
        assert!( excl.should_reduce(&MultiVar::from([0,1,0])));
        assert!(!excl.should_reduce(&MultiVar::from([0,0,1])));
        assert!( excl.should_reduce(&MultiVar::from([0,0,2])));
    }

    #[test]
    fn reduce_gens() {
        let l = Link::trefoil();
        let q = 0;
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);
        let cube = KRHorCube::new(Arc::new(data), q, v);

        let g = (0..=3).map(|i| 
            cube.gens(i).collect_vec()
        ).collect_vec();

        assert_eq!(g[0].len(),  1);
        assert_eq!(g[1].len(), 12);
        assert_eq!(g[2].len(), 26);
        assert_eq!(g[3].len(), 15);

        excl.perform_excl(0, 1, 1); // x1 -> -x2

        let rg = g.iter().map(|g|
            reduce(&excl, g)
        ).collect_vec();

        assert_eq!(rg[0].len(), 0);
        assert_eq!(rg[1].len(), 2);
        assert_eq!(rg[2].len(), 7);
        assert_eq!(rg[3].len(), 5);
    }

    #[test]
    fn reduce_gens_2() {
        let l = Link::trefoil();
        let q = 0;
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);
        let cube = KRHorCube::new(Arc::new(data), q, v);

        let g = (0..=3).map(|i| 
            cube.gens(i).collect_vec()
        ).collect_vec();

        assert_eq!(g[0].len(),  1);
        assert_eq!(g[1].len(), 12);
        assert_eq!(g[2].len(), 26);
        assert_eq!(g[3].len(), 15);

        excl.perform_excl(0, 1, 1); // x1 -> -x2
        excl.perform_excl(2, 2, 2); // x2 * x2 -> x0 * x2

        let rg = g.iter().map(|g|
            reduce(&excl, g)
        ).collect_vec();

        assert_eq!(rg[0].len(), 0);
        assert_eq!(rg[1].len(), 0);
        assert_eq!(rg[2].len(), 2);
        assert_eq!(rg[3].len(), 2);
    }

    #[test]
    fn forward() {
        let l = Link::trefoil();
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);

        excl.perform_excl(0, 1, 1); // x1 -> -x2

        assert_eq!(
            excl.forward_x(&gen([0,0,0], [0,0,1], [0,0,0])), 
            KRChain::zero()
        ); // vanish
        assert_eq!(
            excl.forward_x(&gen([1,0,0], [0,0,1], [0,0,0])), 
              KRChain::from(gen([1,0,0], [0,0,1], [0,0,0])) // 1: id
        );
        assert_eq!(
            excl.forward_x(&gen([1,0,0], [0,0,1], [1,0,0])), 
              KRChain::from(gen([1,0,0], [0,0,1], [1,0,0])) // x0: id
        );
        assert_eq!(
            excl.forward_x(&gen([1,0,0], [0,0,1], [0,1,0])), 
              KRChain::from((gen([1,0,0], [0,0,1], [0,0,1]), -R::one())) // x1 -> -x2
        );
        assert_eq!(
            excl.forward_x(&gen([1,0,0], [0,0,1], [0,0,1])), 
              KRChain::from(gen([1,0,0], [0,0,1], [0,0,1])) // x2: id
        );
    }

    #[test]
    fn forward_2() {
        let l = Link::trefoil();
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);

        excl.perform_excl(0, 1, 1); // x1 -> -x2
        excl.perform_excl(2, 2, 2); // x2 * x2 -> x0 * x2

        assert_eq!(
            excl.forward_x(&gen([0,1,1], [0,0,1], [0,0,0])), 
            KRChain::zero()
        ); // vanish
        assert_eq!(
            excl.forward_x(&gen([1,1,0], [0,0,1], [0,0,0])), 
            KRChain::zero()
        ); // vanish

        assert_eq!(
            excl.forward_x(&gen([1,0,1], [0,0,1], [0,0,0])), 
              KRChain::from(gen([1,0,1], [0,0,1], [0,0,0])) // 1: id
        );
        assert_eq!(
            excl.forward_x(&gen([1,0,1], [0,0,1], [1,0,0])), 
              KRChain::from(gen([1,0,1], [0,0,1], [1,0,0])) // x1: id
        );
        assert_eq!(
            excl.forward_x(&gen([1,0,1], [0,0,1], [0,1,0])), 
              KRChain::from((gen([1,0,1], [0,0,1], [0,0,1]), -R::one())) // x1 -> -x2
        );
        assert_eq!(
            excl.forward_x(&gen([1,0,1], [0,0,1], [0,0,1])), 
              KRChain::from(gen([1,0,1], [0,0,1], [0,0,1])) // x2: id
        );
        assert_eq!(
            excl.forward_x(&gen([1,0,1], [0,0,1], [0,0,2])), 
              KRChain::from(gen([1,0,1], [0,0,1], [1,0,1])) // x2^2 -> x0 x2
        );
        assert_eq!(
            excl.forward_x(&gen([1,0,1], [0,0,1], [1,1,1])), 
              KRChain::from((gen([1,0,1], [0,0,1], [2,0,1]), -R::one())) // x0x1x2 -> -x0^2 x2
        );
    }

    #[test]
    fn backward() {
        let l = Link::trefoil();
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);
        
        excl.perform_excl(0, 1, 1); // x1 -> -x2

        assert_eq!(
            excl.backward_x(&gen([1,0,0], [0,0,1], [0,0,0])), // 1 at [1,0,0]
            KRChain::from_iter([
                (gen([1,0,0], [0,0,1], [0,0,0]),  R::one()), //   1 at [1,0,0]
                (gen([0,1,0], [0,0,1], [0,0,0]), -R::one()), //  -1 at [0,1,0]
                (gen([0,0,1], [0,0,1], [0,0,1]), -R::one())  // -x2 at [0,0,1]
            ])
        );

        assert_eq!(
            excl.backward_x(&gen([1,1,0], [0,0,1], [0,0,0])), // 1 at [1,1,0]
            KRChain::from_iter([
                (gen([1,1,0], [0,0,1], [0,0,0]), R::one()), //   1 at [1,1,0]
                (gen([0,1,1], [0,0,1], [0,0,1]), R::one())  //  x2 at [0,1,1]
            ])
        );

        assert_eq!(
            excl.backward_x(&gen([1,0,1], [0,0,1], [0,0,0])), 
            KRChain::from_iter([
                (gen([1,0,1], [0,0,1], [0,0,0]),  R::one()), //  1 at [1,0,1]
                (gen([0,1,1], [0,0,1], [0,0,0]), -R::one()), // -1 at [0,1,1]
            ])
        );

        assert_eq!(
            excl.backward_x(&gen([1,1,1], [0,0,1], [0,0,0])), 
               KRChain::from(gen([1,1,1], [0,0,1], [0,0,0])), //  id at [1,1,1]
        );
    }

    #[test]
    fn backward_2() {
        let l = Link::trefoil();
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);
        
        excl.perform_excl(0, 1, 1); // x1 -> -x2
        excl.perform_excl(1, 2, 1); // x2 -> 0

        assert_eq!(
            excl.backward_x(&gen([1,1,0], [0,0,1], [0,0,0])), 
            KRChain::from_iter([
                (gen([1,1,0], [0,0,1], [0,0,0]),  R::one()), //     id at [1,1,0]
                (gen([1,0,1], [0,0,1], [1,0,0]), -R::one()), //    -x0
                (gen([1,0,1], [0,0,1], [0,0,1]),  R::one()), //    +x2 at [1,0,1]
                (gen([0,1,1], [0,0,1], [1,0,0]),  R::one()), //     x0 at [0,1,1]
            ])
        );
        assert_eq!(
            excl.backward_x(&gen([1,1,1], [0,0,1], [0,0,0])), 
               KRChain::from(gen([1,1,1], [0,0,1], [0,0,0]))
        );
    }

    #[test]
    fn backward_3() {
        let l = Link::trefoil();
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);
        
        excl.perform_excl(0, 1, 1); // x1 -> -x2
        excl.perform_excl(2, 2, 2); // x2^2 -> x0 x2

        assert_eq!(
            excl.backward_x(&gen([1,0,1], [0,0,1], [0,0,0])), 
            KRChain::from_iter([
                (gen([1,0,1], [0,0,1], [0,0,0]),  R::one()),
                (gen([0,1,1], [0,0,1], [0,0,0]), -R::one()),
            ])
        );
        assert_eq!(
            excl.backward_x(&gen([1,1,1], [0,0,1], [0,0,0])), 
               KRChain::from(gen([1,1,1], [0,0,1], [0,0,0]))
        );
    }

    fn make_trans<R>(excl: &KRHorExcl<R>, gens: &Vec<KRGen>) -> Trans<R>
    where R: Ring, for<'x> &'x R: RingOps<R> { 
        let from = IndexList::from_iter( gens.iter().cloned() );
        let to = IndexList::from_iter( reduce(excl, gens) );
        excl.trans_for(&from, &to)
    }

    #[test]
    fn trans_before() { 
        let l = Link::trefoil();
        let q = 0;
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let excl = KRHorExcl::from(&data, v, 0);
        let cube = KRHorCube::new(Arc::new(data), q, v);

        let gens = cube.gens(2).collect_vec();
        let t = make_trans(&excl, &gens);

        assert_eq!(gens.len(), 26);
        assert_eq!(t.forward_mat().shape(), (26, 26));
        assert!(t.forward_mat().is_id());
        assert!(t.backward_mat().is_id());
    }

    #[test]
    fn trans() { 
        let l = Link::trefoil();
        let q = 0;
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);
        let cube = KRHorCube::new(Arc::new(data), q, v);

        let gens = cube.gens(2).collect_vec();
        excl.perform_excl(0, 1, 1); // x1 -> -x2

        let t = make_trans(&excl, &gens);

        assert_eq!(t.forward_mat().shape(), (7, 26));
        assert!((t.forward_mat() * t.backward_mat()).is_id());
    }

    #[test]
    fn trans_2() { 
        let l = Link::trefoil();
        let q = 0;
        let v = BitSeq::from_iter([0,0,1]);

        let data = KRCubeData::<R>::new_no_excl(&l);
        let mut excl = KRHorExcl::from(&data, v, 0);
        let cube = KRHorCube::new(Arc::new(data), q, v);

        let gens = cube.gens(2).collect_vec();
        excl.perform_excl(0, 1, 1); // x1 -> -x2
        excl.perform_excl(2, 2, 2); // x2 * x2 -> x0 * x2

        let t = make_trans(&excl, &gens);

        assert_eq!(t.forward_mat().shape(), (2, 26));
        assert!((t.forward_mat() * t.backward_mat()).is_id());
    }
}