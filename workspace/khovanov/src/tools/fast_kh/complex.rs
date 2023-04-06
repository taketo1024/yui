use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use std::ops::RangeInclusive;

use itertools::{Itertools, join};
use num_traits::Zero;
use cartesian::cartesian;
use yui_core::{Ring, RingOps};
use yui_homology::{GenericChainComplex, FreeRModStr};
use yui_lin_comb::LinComb;
use yui_link::{Crossing, Resolution, Link};

use crate::{KhAlgLabel, KhComplex, KhGen};
use super::cob::{Cob, Dot, Bottom, CobComp};
use super::tng::{Tng, TngComp};

// MEMO remove 
#[cfg(not(test))] 
use log::info; 
#[cfg(test)]
use std::println as info;

type Mor = LinComb<Cob, i32>; // Z-linear combination of cobordisms.

trait MorTrait: Sized {
    fn is_invertible(&self) -> bool;
    fn inv(&self) -> Option<Self>;
    fn map_cob<F>(self, f: F) -> Self where F: Fn(&mut Cob);
    fn connect(self, c: &Cob) -> Self;
    fn connect_comp(self, c: &CobComp) -> Self;
    fn cap_off(self, b: Bottom, c: &TngComp, dot: Dot) -> Self;
    fn eval<R>(&self, h: &R, t: &R) -> R
    where R: Ring + From<i32>, for<'x> &'x R: RingOps<R>;
}

impl MorTrait for Mor {
    fn is_invertible(&self) -> bool { 
        self.len() == 1 && 
        self.iter().next().map(|(c, a)| 
            c.is_invertible() && a.is_unit()
        ).unwrap_or(false)
    }

    fn inv(&self) -> Option<Self> { 
        if let Some((Some(cinv), Some(ainv))) = self.iter().next().map(|(c, a)| 
            (c.inv(), a.inv())
        ) { 
            let inv = Mor::from((cinv, ainv));
            Some(inv)
        } else { 
            None
        }
    }

    fn map_cob<F>(self, f: F) -> Self 
    where F: Fn(&mut Cob) {
        self.into_map(|mut cob, r| { 
            f(&mut cob);
            if cob.is_zero() { 
                (cob, 0)
            } else { 
                (cob, r)
            }
        })
    }

    fn connect(self, c: &Cob) -> Self {
        self.map_cob(|cob| cob.connect(c.clone()) )
    }

    fn connect_comp(self, c: &CobComp) -> Self {
        self.map_cob(|cob| cob.connect_comp(c.clone()) )
    }

    fn cap_off(self, b: Bottom, c: &TngComp, dot: Dot) -> Self {
        self.map_cob(|cob| cob.cap_off(b, c, dot) )
    }

    fn eval<R>(&self, h: &R, t: &R) -> R
    where R: Ring + From<i32>, for<'x> &'x R: RingOps<R> {
        self.iter().map(|(c, &a)| { 
            R::from(a) * c.eval(h, t)
        }).sum()
    }
}

#[derive(Clone, Debug)]
pub struct TngVertex { 
    key: KhGen,
    tng: Tng,
    in_edges: HashSet<KhGen>,
    out_edges: HashMap<KhGen, Mor>
}

impl TngVertex { 
    pub fn init() -> Self { 
        let key = KhGen::init();
        let tng = Tng::empty();
        let in_edges = HashSet::new();
        let out_edges = HashMap::new();
        Self { key, tng, in_edges, out_edges }
    }
}

impl Display for TngVertex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let key = {
            let state = join(self.key.state().values().iter(), "");
            if self.key.label().is_empty() { 
                state
            } else { 
                let label = join(self.key.label().iter(), "");
                format!("{}-{}", state, label)
            }
        };

        write!(f, "({}; {})", key, self.tng)
    }
}

pub struct TngComplex {
    vertices: HashMap<KhGen, TngVertex>,
    len: usize,
    deg_shift: (isize, isize)
}

impl TngComplex {
    pub fn new(deg_shift: (isize, isize)) -> Self { 
        let mut vertices = HashMap::new();
        let k0 = KhGen::init();
        let v0 = TngVertex::init();
        vertices.insert(k0, v0);

        let len = 1;

        TngComplex{ vertices, len, deg_shift }
    }

    pub fn len(&self) -> usize { 
        self.len
    }

    pub fn rank(&self, i: isize) -> usize { 
        let i0 = self.deg_shift.0;
        self.vertices.keys().filter(|k| {
            let w = k.state().weight() as isize;
            w == i - i0
        }).count()
    }

    pub fn vertex(&self, v: &KhGen) -> &TngVertex { 
        &self.vertices[v]
    }

    pub fn nverts(&self) -> usize { 
        self.vertices.len()
    }

    pub fn iter_verts(&self) -> impl Iterator<Item = (&KhGen, &TngVertex)> {
        self.vertices.iter().sorted_by(|(k0, _), (k1, _)| k0.cmp(k1))
    }

    //                          d0
    //                 C[i]#x0 ---> C[i+1]#x0
    //                   |             :
    //                 f |             :
    //            -d1    V             :
    //  C[i-1]#x1 ---> C[i]#x1 .... C[i+1]#x1 
    //
    //  C'[i] = C[i]#x0 ⊕ C[i-1]#x1
    //     d' = [d0    ]
    //          [f  -d1]

    pub fn append(&mut self, x: &Crossing) {
        info!("({}) append: {x}", self.nverts());

        let sdl = CobComp::from(x);

        let c0 = std::mem::take(&mut self.vertices);
        let c1 = c0.clone();

        let c0 = Self::modify_vertices(c0, &sdl, sdl.src(), Resolution::Res0);
        let c1 = Self::modify_vertices(c1, &sdl, sdl.tgt(), Resolution::Res1);

        self.vertices.extend(c0);
        self.vertices.extend(c1);
        self.len += 1;
    }

    fn modify_vertices(vertices: HashMap<KhGen, TngVertex>, sdl: &CobComp, t: &Tng, r: Resolution) -> HashMap<KhGen, TngVertex> {
        let c = Cob::id(t); // t × I

        vertices.into_iter().map(|(k, mut v)| { 
            let v_in  = std::mem::take(&mut v.in_edges);
            v.in_edges = v_in.into_iter().map(|mut k| {
                k.append_state(r);
                k
            }).collect();
    
            let v_out = std::mem::take(&mut v.out_edges);
            v.out_edges = v_out.into_iter().map(|(mut k, mut v)| { 
                k.append_state(r);
                v = v.connect(&c);
                
                if r == Resolution::Res0 { 
                    (k, v)
                } else { 
                    (k, -v)
                }
            }).collect();

            if r == Resolution::Res0 { 
                let mut k1 = k.clone();
                k1.append_state(Resolution::Res1);

                let mut f = Cob::id(&v.tng);
                f.connect_comp(sdl.clone());
                v.out_edges.insert(k1.clone(), Mor::from(f));

            } else { 
                let mut k0 = k.clone();
                k0.append_state(Resolution::Res0);
                v.in_edges.insert(k0);
            }
    
            v.key.append_state(r);
            v.tng.connect(t.clone());

            (v.key.clone(), v)
        }).collect()
    }

    pub fn deloop(&mut self, simplify: bool) { 
        while let Some((k, r)) = self.find_loop() { 
            let (k0, k1) = self.deloop_at(&k, r);
            if simplify { 
                self.simplify_on_deloop(&k0);
                self.simplify_on_deloop(&k1);
            }
        }
    }

    fn find_loop(&self) -> Option<(KhGen, usize)> { 
        for (k, v) in self.iter_verts() { 
            if let Some(r) = v.tng.find_loop() { 
                return Some((k.clone(), r))
            }
        }
        None
    }

    fn deloop_at(&mut self, k: &KhGen, r: usize) -> (KhGen, KhGen) { 
        info!("({}) deloop {} at {r}", self.nverts(), &self.vertices[k]);

        let mut v0 = self.vertices.remove(k).unwrap();
        let c = v0.tng.remove_at(r);
        let mut v1 = v0.clone();

        v0.key.append_label(KhAlgLabel::X);
        v1.key.append_label(KhAlgLabel::I);

        let k0 = &v0.key;
        let k1 = &v1.key;

        let v0_in = v0.in_edges.iter().cloned().collect_vec();
        for j in v0_in.iter() { 
            let u = self.vertices.get_mut(j).unwrap();

            let f = u.out_edges.remove(&k).unwrap();
            v0.in_edges.remove(&j);
            v1.in_edges.remove(&j);

            let f0 = f.clone().cap_off(Bottom::Tgt, &c, Dot::None);
            let f1 = f.cap_off(Bottom::Tgt, &c, Dot::Y);

            if !f0.is_zero() {
                u.out_edges.insert(k0.clone(), f0);
                v0.in_edges.insert(j.clone());
            }

            if !f1.is_zero() {
                u.out_edges.insert(k1.clone(), f1);
                v1.in_edges.insert(j.clone());
            }
        }
        
        let v0_out = v0.out_edges.keys().cloned().collect_vec();
        for l in v0_out.iter() { 
            let w = self.vertices.get_mut(l).unwrap();

            let f0 = v0.out_edges.remove(&l).unwrap();
            let f1 = v1.out_edges.remove(&l).unwrap();
            w.in_edges.remove(&k);

            let f0 = f0.cap_off(Bottom::Src, &c, Dot::X);
            let f1 = f1.cap_off(Bottom::Src, &c, Dot::None);

            if !f0.is_zero() { 
                v0.out_edges.insert(l.clone(), f0);
                w.in_edges.insert(k0.clone());
            }

            if !f1.is_zero() { 
                v1.out_edges.insert(l.clone(), f1);
                w.in_edges.insert(k1.clone());
            }
        }
        
        let k0 = k0.clone();
        let k1 = k1.clone();

        self.vertices.insert(k0.clone(), v0);
        self.vertices.insert(k1.clone(), v1);

        (k0, k1)
    }

    fn simplify_on_deloop(&mut self, k: &KhGen) { 
        let v = self.vertex(k);
        
        for k_from in v.in_edges.iter() { 
            let w = &self.vertices[k_from];
            let f = &w.out_edges[k];
            if f.is_invertible() { 
                self.eliminate(&k_from.clone(), k);
                return;
            }
        }

        for (k_to, f) in v.out_edges.iter() { 
            if f.is_invertible() { 
                self.eliminate(k, &k_to.clone());
                return;
            }
        }
    }

    fn eliminate(&mut self, k0: &KhGen, l0: &KhGen) {

        let v0 = self.vertex(k0);
        let w0 = self.vertex(l0);
        let a = &v0.out_edges[l0];

        println!("({}) eliminate {}: {} -> {}", self.nverts(), a, &v0, &w0);

        let Some(ainv) = a.inv() else { 
            panic!()
        };
        
        //       a
        //  v0 - - -> w0
        //     \   / b
        //       / 
        //     /   \ c
        //  v1 -----> w1
        //       d

        let w0_in  = w0.in_edges.iter().cloned().collect_vec();
        let v0_out = v0.out_edges.keys().cloned().collect_vec();
        
        for (k1, l1) in cartesian!(w0_in.iter(), v0_out.iter()) {
            if k1 == k0 || l1 == l0 { 
                continue
            }

            let b = &self.vertex(k1).out_edges[l0];
            let c = &self.vertex(k0).out_edges[l1];
            
            let cab = c * &ainv * b;
            let v1 = self.vertices.get_mut(&k1).unwrap();

            let d = if let Some(d) = v1.out_edges.get(&l1) { 
                d - cab
            } else { 
                -cab
            };

            if d.is_zero() { 
                self.remove_edge(k1, l1);
            } else { 
                self.add_edge(k1, l1, d);
            }
        }

        // remove edges into w0
        for k1 in w0_in.iter() {
            self.remove_edge(k1, l0);
        }

        // remove edges out from v0
        for l1 in v0_out.iter() {
            self.remove_edge(k0, l1);
        }

        let v0_in  = self.vertex(k0).in_edges.iter().cloned().collect_vec();
        let w0_out = self.vertex(l0).out_edges.keys().cloned().collect_vec();

        // remove edges into v0
        for j in v0_in.iter() { 
            self.remove_edge(j, k0);
        }

        // remove edges out from w0
        for j in w0_out.iter() { 
            self.remove_edge(l0, j);
        }

        // remove v0 and w0.
        self.vertices.remove(k0);
        self.vertices.remove(l0);
    }

    fn add_edge(&mut self, k: &KhGen, l: &KhGen, f: Mor) { 
        let v = self.vertices.get_mut(k).unwrap();
        v.out_edges.insert(l.clone(), f);

        let w = self.vertices.get_mut(l).unwrap();
        w.in_edges.insert(k.clone());
    }

    fn remove_edge(&mut self, k: &KhGen, l: &KhGen) { 
        let v = self.vertices.get_mut(k).unwrap();
        v.out_edges.remove(l);

        let w = self.vertices.get_mut(l).unwrap();
        w.in_edges.remove(k);
    }

    pub fn is_completely_delooped(&self) -> bool { 
        self.vertices.iter().all(|(_, v)|
            v.tng.is_empty()
        )
    }

    pub fn as_generic<R>(&self, h: R, t: R) -> GenericChainComplex<R, RangeInclusive<isize>> 
    where R: Ring + From<i32>, for<'x> &'x R: RingOps<R> {
        assert!(self.is_completely_delooped());

        // TODO improve construction. 

        let n = self.len();
        let c = (0..n).map(|i| {
            let gens = self.vertices.keys().filter(|k| k.state().weight() == i).sorted().cloned().collect();
            FreeRModStr::new(gens)
        }).collect_vec();

        let i0 = self.deg_shift.0;

        GenericChainComplex::ascending_from(i0, (0..n-1).map( |i| {
            c[i].make_matrix(&c[i+1], |x| { 
                let v = self.vertex(x);
                v.out_edges.iter().map(|(y, f)| 
                    (y.clone(), f.eval(&h, &t))
                ).collect()
            })
        }).collect())
    }

    pub fn describe(&self) { 
        let mut str = "".to_string();
        for k0 in self.vertices.keys().sorted() { 
            let v = &self.vertices[&k0];
            str += &format!("{k0}: {}\n", v.tng);

            for k1 in v.out_edges.keys().sorted() { 
                let f = &v.out_edges[&k1];
                str += &format!(" -> {k1}: {f}\n");
            }
            str += "\n";
        }
        println!("{str}");
    }

    fn validate_edges(&self) {
        for (k, v) in self.vertices.iter() { 
            for j in v.in_edges.iter() {
                assert!(
                    self.vertices.contains_key(j),
                    "no vertex for in-edge {j} -> {k}"
                );
                let u = self.vertex(j);
                assert!(
                    u.out_edges.contains_key(k),
                    "no out-edge {j} -> {k}"
                );
            }
            for (l, _f) in v.out_edges.iter() {
                assert!(
                    self.vertices.contains_key(l),
                    "no vertex for out-edge {k} -> {l}"
                );
                let w = self.vertex(l);
                assert!(
                    w.in_edges.contains(k),
                    "no in-edge {k} -> {l}"
                );
                // assert!(!f.is_zero());
            }
        }
    }
}

impl From<&Link> for TngComplex {
    fn from(l: &Link) -> Self {
        info!("construct TngComplex.");

        let deg_shift = KhComplex::<i64>::deg_shift_for(l, false);
        let mut c = TngComplex::new(deg_shift);
        for x in l.data() {
            c.append(x);
            c.deloop(true);
        }

        info!("result: #v = {}", c.vertices.len());

        c
    }
}

#[test]
fn test() { 
    let l = Link::from(&[[1, 30, 2, 31], [3, 27, 4, 26], [5, 12, 6, 13], [7, 29, 8, 28], [9, 35, 10, 34], [11, 4, 12, 5], [13, 8, 14, 9], [14, 36, 15, 35], [17, 24, 18, 25], [18, 32, 19, 31], [21, 3, 22, 2], [23, 16, 24, 17], [25, 22, 26, 23], [27, 7, 28, 6], [29, 21, 30, 20], [32, 16, 33, 15], [33, 11, 34, 10], [36, 20, 1, 19]]);
    let c = TngComplex::from(&l);
    info!("{:?}", c.iter_verts().map(|x| x.0.state().to_string()).join(", "));
}

#[cfg(test)]
mod tests { 
    use yui_link::*;
    use yui_homology::*;
    use yui_homology::test::ChainComplexValidation;
    use super::*;

    #[test]
    fn mor_inv() { 
        let c = Cob::id(&Tng::new(vec![
            TngComp::arc(0,1),
            TngComp::arc(2,3)
        ]));
        let f = Mor::from((c.clone(), -1));

        assert!(f.is_invertible());
        assert_eq!(f.inv(), Some(f.clone()));

        let f = Mor::from((c.clone(), 2));
        assert_eq!(f.is_invertible(), false);
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn empty() { 
        let c = TngComplex::new((0, 0));

        assert_eq!(c.len(), 1);
        assert_eq!(c.rank(0), 1);
    }

    #[test]
    fn single_x() { 
        let mut c = TngComplex::new((0, 0));
        let x = Crossing::new(CrossingType::Xn, [0,1,2,3]);
        c.append(&x);

        assert_eq!(c.len(), 2);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 1);
    }

    #[test]
    fn two_x_disj() { 
        let mut c = TngComplex::new((0, 0));
        let x0 = Crossing::new(CrossingType::Xn, [0,1,2,3]);
        let x1 = Crossing::new(CrossingType::Xn, [4,5,6,7]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.len(), 3);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);
    }

    #[test]
    fn two_x() { 
        let mut c = TngComplex::new((0, 0));
        let x0 = Crossing::new(CrossingType::Xn, [0,4,1,5]);
        let x1 = Crossing::new(CrossingType::Xn, [3,1,4,2]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.len(), 3);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);
    }

    #[test]
    fn deloop_one() { 
        let mut c = TngComplex::new((0, 0));
        let x0 = Crossing::new(CrossingType::Xp, [1,4,2,5]);
        let x1 = Crossing::new(CrossingType::Xn, [3,6,4,1]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.len(), 3);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);

        c.deloop(false);

        assert_eq!(c.len(), 3);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 3); // delooped here
        assert_eq!(c.rank(2), 1);
    }

    #[test]
    fn trefoil_no_deloop() { 
        let mut c = TngComplex::new((0, 0));
        let l = Link::trefoil();

        for x in l.data() {
            c.append(x);
        }

        assert_eq!(c.len(), 4);

        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 3);
        assert_eq!(c.rank(2), 3);
        assert_eq!(c.rank(3), 1);
    }

    #[test]
    fn trefoil_deloop() { 
        let mut c = TngComplex::new((0, 0));
        let l = Link::trefoil();

        for x in l.data() {
            c.append(x);
            c.deloop(false);
        }

        assert_eq!(c.len(), 4);

        assert_eq!(c.rank(0), 8);
        assert_eq!(c.rank(1), 12);
        assert_eq!(c.rank(2), 6);
        assert_eq!(c.rank(3), 4);

        let c = c.as_generic(0, 0);
        let h = c.homology(); // TODO should shift degree

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 1);
        assert_eq!(h[1].tors(), &vec![2]);
        
        assert_eq!(h[2].is_zero(), true);
        assert_eq!(h[2].is_free(), true);
        
        assert_eq!(h[3].rank(), 2);
        assert_eq!(h[3].is_free(), true);
    }

    #[test]
    fn trefoil_deloop_simplify() { 
        let mut c = TngComplex::new((0, 0));
        let l = Link::trefoil();

        for x in l.data() {
            c.append(x);
            c.deloop(true);
        }

        assert_eq!(c.len(), 4);

        assert_eq!(c.rank(0), 2);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 0);
        assert_eq!(c.rank(3), 2);

        let c = c.as_generic(0, 0);
        let h = c.homology(); // TODO should shift degree

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 1);
        assert_eq!(h[1].tors(), &vec![2]);
        
        assert_eq!(h[2].is_zero(), true);
        assert_eq!(h[2].is_free(), true);
        
        assert_eq!(h[3].rank(), 2);
        assert_eq!(h[3].is_free(), true);
    }

    #[test]
    fn test_unknot_rm1() {
        let l = Link::from(&[[0,0,1,1]]);
        let c = TngComplex::from(&l);
        let c = c.as_generic(0, 0);

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);

        let l = Link::from(&[[0,1,1,0]]);
        let c = TngComplex::from(&l);
        let c = c.as_generic(0, 0);

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[0].rank(), 2);
    }

    #[test]
    fn test_unknot_rm2() {
        let l = Link::from(&[[1,4,2,1],[2,4,3,3]]);
        let c = TngComplex::from(&l);
        let c = c.as_generic(0, 0);

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_unlink_2() {
        let pd_code = [[1,2,3,4], [3,2,1,4]];
        let l = Link::from(&pd_code);
        let c = TngComplex::from(&l);
        let c = c.as_generic(0, 0);

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_hopf_link() {
        let l = Link::hopf_link();
        let c = TngComplex::from(&l);
        let c = c.as_generic(0, 0);

        assert_eq!(c[-2].rank(), 2);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[0].rank(), 2);
    }

    #[test]
    fn test_8_19() {
        let l = Link::from(&[[4,2,5,1],[8,4,9,3],[9,15,10,14],[5,13,6,12],[13,7,14,6],[11,1,12,16],[15,11,16,10],[2,8,3,7]]);
        let c = TngComplex::from(&l);
        let c = c.as_generic(0, 0);

        c.check_d_all();

        let h = c.homology();

        for i in [1,6,7,8] {
            assert_eq!(h[i].rank(), 0);
            assert_eq!(h[i].is_free(), true);
        }

        for i in [0,4,5] {
            assert_eq!(h[i].rank(), 2);
            assert_eq!(h[i].is_free(), true);
        }

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);

        assert_eq!(h[3].rank(), 1);
        assert_eq!(h[3].tors(), &vec![2]);
    }
}
