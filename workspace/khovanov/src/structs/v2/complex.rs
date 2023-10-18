use std::collections::{HashMap, HashSet};
use std::fmt::Display;

use log::info; 
use itertools::Itertools;
use num_traits::Zero;
use cartesian::cartesian;
use yui_core::{Ring, RingOps};
use yui_homology::v2::XChainComplex;
use yui_link::{Crossing, Resolution, Edge};

use crate::{KhAlgGen, KhEnhState};
use super::cob::{Cob, Dot, Bottom, CobComp};
use super::tng::{Tng, TngComp};
use super::mor::{Mor, MorTrait};

#[derive(Clone, Debug)]
pub struct TngVertex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    key: KhEnhState,
    tng: Tng,
    in_edges: HashSet<KhEnhState>,
    out_edges: HashMap<KhEnhState, Mor<R>>
}

impl<R> TngVertex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn init() -> Self { 
        let key = KhEnhState::init();
        let tng = Tng::empty();
        let in_edges = HashSet::new();
        let out_edges = HashMap::new();
        Self { key, tng, in_edges, out_edges }
    }

    pub fn tng(&self) -> &Tng { 
        &self.tng
    }

    pub fn out_edges(&self) -> &HashMap<KhEnhState, Mor<R>> {
        &self.out_edges
    }
}

impl<R> Display for TngVertex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}; {})", self.key, self.tng)
    }
}

pub struct TngComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    vertices: HashMap<KhEnhState, TngVertex<R>>,
    len: usize,
    h: R,
    t: R,
    deg_shift: (isize, isize)
}

impl<R> TngComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(h: &R, t: &R, deg_shift: (isize, isize)) -> Self { 
        let mut vertices = HashMap::new();
        let k0 = KhEnhState::init();
        let v0 = TngVertex::init();
        vertices.insert(k0, v0);

        let len = 1;
        let h = h.clone();
        let t = t.clone();

        TngComplex{ vertices, len, h, t, deg_shift }
    }

    pub fn len(&self) -> usize { 
        self.len
    }

    pub fn rank(&self, i: isize) -> usize { 
        let i0 = self.deg_shift.0;
        self.vertices.keys().filter(|k| {
            let w = k.state.weight() as isize;
            w == i - i0
        }).count()
    }

    pub fn vertex(&self, v: &KhEnhState) -> &TngVertex<R> { 
        &self.vertices[v]
    }

    pub fn nverts(&self) -> usize { 
        self.vertices.len()
    }

    pub fn iter_verts(&self) -> impl Iterator<Item = (&KhEnhState, &TngVertex<R>)> {
        self.vertices.iter().sorted_by(|(k0, _), (k1, _)| k0.cmp(k1))
    }

    pub fn edge(&self, k: &KhEnhState, l: &KhEnhState) -> &Mor<R> {
        &self.vertices[k].out_edges[l]
    }

    fn add_edge(&mut self, k: &KhEnhState, l: &KhEnhState, f: Mor<R>) { 
        let v = self.vertices.get_mut(k).unwrap();
        v.out_edges.insert(*l, f);

        let w = self.vertices.get_mut(l).unwrap();
        w.in_edges.insert(*k);
    }

    fn remove_edge(&mut self, k: &KhEnhState, l: &KhEnhState) { 
        let v = self.vertices.get_mut(k).unwrap();
        v.out_edges.remove(l);

        let w = self.vertices.get_mut(l).unwrap();
        w.in_edges.remove(k);
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

        if !x.is_resolved() { 
            self.append_x(x)
        } else { 
            self.append_a(x)
        }
    }

    fn append_x(&mut self, x: &Crossing) {
        assert!(!x.is_resolved());

        let verts = std::mem::take(&mut self.vertices);

        let sdl = CobComp::from(x);
        let c0 = Cob::id(sdl.src());
        let c1 = Cob::id(sdl.tgt());

        let mut rmv = HashSet::new();

        for (_, v) in verts.into_iter() { 
            let vs = Self::make_cone(v, &c0, &c1, &sdl);
            for v in vs { 
                for (l, f) in v.out_edges.iter() { 
                    if f.is_zero() { 
                        rmv.insert((v.key, *l));
                    }
                }
                self.vertices.insert(v.key, v);
            }
        }

        for (k, l) in rmv { 
            self.remove_edge(&k, &l);
        }

        self.len += 1;

        debug_assert!(self.validate_edges());
    }

    fn append_a(&mut self, x: &Crossing) {
        assert!(x.is_resolved());

        let verts = std::mem::take(&mut self.vertices);
        let c = Cob::id(&Tng::from_a(x));

        self.vertices = verts.into_iter().map(|(k, mut v)| { 
            v.tng.connect(c.src());

            modify!(v.out_edges, |edges: HashMap<KhEnhState, Mor<R>>| { 
                edges.into_iter().map(|(k, f)| { 
                    (k, f.connect(&c))
                }).collect()
            });
    
            (k, v)
        }).collect();

        debug_assert!(self.validate_edges());
    }

    fn make_cone(v: TngVertex<R>, c0: &Cob, c1: &Cob, sdl: &CobComp) -> [TngVertex<R>; 2] {
        use Resolution::{Res0, Res1};

        let mut v0 = v.clone();
        let mut v1 = v;

        let mut c = Cob::id(&v0.tng);
        c.connect_comp(sdl.clone());
        
        Self::extend_by(&mut v0, c0, Res0);
        Self::extend_by(&mut v1, c1, Res1);
        Self::insert_sdl(&mut v0, &mut v1, c);

        [v0, v1]
    }

    fn extend_by(v: &mut TngVertex<R>, c: &Cob, r: Resolution) {
        v.key.state.push(r);
        v.tng.connect(c.src());

        // append 0 / 1 to in_edges.
        modify!(v.in_edges, |edges: HashSet<KhEnhState>| { 
            edges.into_iter().map(|mut k| { 
                k.state.push(r);
                k
            }).collect()
        });

        // append 0 / 1 to out_edges with id-cob connected.
        let e = if r.is_zero() { R::one() } else { -R::one() };
        modify!(v.out_edges, |edges: HashMap<KhEnhState, Mor<R>>| { 
            edges.into_iter().map(|(mut k, f)| { 
                k.state.push(r);
                (k, f.connect(c) * &e)
            }).collect()
        });
    }

    fn insert_sdl(v0: &mut TngVertex<R>, v1: &mut TngVertex<R>, sdl: Cob) { 
        let f = Mor::from_gen(sdl);
        v0.out_edges.insert(v1.key, f);
        v1.in_edges.insert(v0.key);
    }

    pub fn find_loop(&self, exclude: Option<Edge>) -> Option<(KhEnhState, usize, &TngComp)> { 
        for (k, v) in self.iter_verts() { 
            if let Some((r, c)) = v.tng.find_loop(exclude) { 
                return Some((*k, r, c))
            }
        }
        None
    }

    pub fn deloop(&mut self, k: &KhEnhState, r: usize, reduced: bool) -> Vec<KhEnhState> { 
        info!("({}) deloop {} at {r}", self.nverts(), &self.vertices[k]);

        let mut v = self.vertices.remove(k).unwrap();
        let c = v.tng.remove_at(r);

        if reduced { 
            self.deloop_in(&mut v, &c, KhAlgGen::X, Dot::X, Dot::None, true);

            let k_new = v.key;
            self.vertices.insert(k_new, v);

            debug_assert!(self.validate_edges());

            vec![k_new]
        } else { 
            let mut v0 = v;
            let mut v1 = v0.clone();

            self.deloop_in(&mut v0, &c, KhAlgGen::X, Dot::X, Dot::None, false);
            self.deloop_in(&mut v1, &c, KhAlgGen::I, Dot::None, Dot::Y, true);

            let k0 = v0.key;
            let k1 = v1.key;

            self.vertices.insert(k0, v0);
            self.vertices.insert(k1, v1);

            debug_assert!(self.validate_edges());

            vec![k0, k1]
        }
    }

    fn deloop_in(&mut self, v: &mut TngVertex<R>, c: &TngComp, label: KhAlgGen, birth_dot: Dot, death_dot: Dot, remove: bool) { 
        assert!(c.is_circle());

        let k_old = v.key;

        v.key.label.push(label);
        let k_new = v.key;

        let v_in = v.in_edges.iter().cloned().collect_vec();
        for j in v_in.iter() { 
            v.in_edges.remove(&j);
            let u = self.vertices.get_mut(j).unwrap();

            let f_old = if remove {
                u.out_edges.remove(&k_old).unwrap()
            } else { 
                u.out_edges[&k_old].clone()
            };

            let (h, t) = (&self.h, &self.t);
            let f_new = f_old.cap_off(Bottom::Tgt, &c, death_dot).part_eval(h, t);

            if !f_new.is_zero() {
                u.out_edges.insert(k_new, f_new);
                v.in_edges.insert(*j);
            }
        }
        
        let v_out = v.out_edges.keys().cloned().collect_vec();
        for l in v_out.iter() { 
            let w = self.vertices.get_mut(l).unwrap();
            if remove {
                w.in_edges.remove(&k_old);
            }

            let f_old = v.out_edges.remove(&l).unwrap();

            let (h, t) = (&self.h, &self.t);
            let f_new = f_old.cap_off(Bottom::Src, &c, birth_dot).part_eval(h, t);

            if !f_new.is_zero() { 
                v.out_edges.insert(*l, f_new);
                w.in_edges.insert(k_new);
            }
        }
    }

    pub fn find_inv_edge(&self, k: &KhEnhState) -> Option<(KhEnhState, KhEnhState)> { 
        let mut cand = None;
        let mut cand_s = usize::MAX;

        fn score<_R>(v: &TngVertex<_R>, w: &TngVertex<_R>) -> usize
        where _R: Ring, for<'x> &'x _R: RingOps<_R> { 
            let v_out = v.out_edges.len();
            let w_in  = w.in_edges.len();
            (v_out - 1) * (w_in - 1)
        }

        // collect candidate edges into v, and out from v. 

        let v = self.vertex(k);
        let edges = v.in_edges.iter().map(|j| (j, k)).chain(
            v.out_edges.keys().map(|l| (k, l))
        );

        for (k, l) in edges { 
            let f = self.edge(k, l);
            if f.is_invertible() { 
                let v = self.vertex(k);
                let w = self.vertex(l);
                let s = score(v, w);

                if s == 0 {
                    info!("({}) good pivot {}: {} -> {}", self.nverts(), f, v, w);
                    return Some((*k, *l));
                } else if s < cand_s { 
                    cand = Some((k, l));
                    cand_s = s;
                }
            }
        }

        if let Some((k, l)) = cand { 
            info!("({}) pivot (score: {cand_s}) {}: {} -> {}", self.nverts(), self.edge(k, l), self.vertex(k), self.vertex(l));
            Some((*k, *l))
        } else { 
            None
        }
    }

    pub fn eliminate(&mut self, k0: &KhEnhState, l0: &KhEnhState) {
        let v0 = self.vertex(k0);
        let w0 = self.vertex(l0);
        let a = &v0.out_edges[l0];

        info!("({}) eliminate {}: {} -> {}", self.nverts(), a, &v0, &w0);

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

            let (h, t) = (&self.h, &self.t);
            let d = d.part_eval(h, t);

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

        debug_assert!(self.validate_edges());
    }

    pub fn eval(&self, h: &R, t: &R) -> XChainComplex<KhEnhState, R> {
        debug_assert!(self.is_evalable());

        let (h, t) = (h.clone(), t.clone());

        let n = self.len();
        let i0 = self.deg_shift.0;
        let i1 = i0 + (n as isize) - 1;

        let all_gens = self.vertices.keys().cloned().collect_vec();
        let vertices = self.vertices.clone();

        XChainComplex::new(i0..=i1, 1, 
            |i| { 
                let w = (i - i0) as usize;
                all_gens.iter().filter(|x| 
                    x.state.weight() == w
                ).sorted().cloned().collect()
            }, 
            |_i, x| { 
                let v = &vertices[&x];
                v.out_edges.iter().map(|(y, f)| 
                    (*y, f.eval(&h, &t))
                ).collect()
            }
        )
    }

    fn is_evalable(&self) -> bool { 
        self.vertices.iter().all(|(_, v)|
            v.tng.is_empty()
        )
    }

    pub fn desc_d(&self) -> String { 
        let mut str = "".to_string();
        for k0 in self.vertices.keys().sorted() { 
            let v = &self.vertices[&k0];
            str += &format!("{k0}: {}", v.tng);

            for k1 in v.out_edges.keys().sorted() { 
                let f = &v.out_edges[&k1];
                str += &format!("\n  -> {k1}: {f}");
            }
            str += "\n";
        }
        str
    }

    pub fn print_d(&self) { 
        println!("{}", self.desc_d());
    }

    pub fn validate_edges(&self) -> bool {
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
            for (l, f) in v.out_edges.iter() {
                assert!(
                    self.vertices.contains_key(l),
                    "no vertex for out-edge {k} -> {l}"
                );
                let w = self.vertex(l);
                assert!(
                    w.in_edges.contains(k),
                    "no in-edge {k} -> {l}"
                );
                assert!(!f.is_zero());
            }
        }
        true
    }
}

macro_rules! modify {
    ($e: expr, $f: expr) => {{
        let val = std::mem::take(&mut $e);
        $e = $f(val);
    }};
}

pub(self) use modify;

#[cfg(test)]
mod tests { 
    use yui_link::*;
    use crate::KhLabel;
    use super::super::tng::TngComp;

    use super::*;

    #[test]
    fn mor_inv() { 
        let c = Cob::id(&Tng::new(vec![
            TngComp::arc(0,1),
            TngComp::arc(2,3)
        ]));
        let f = Mor::from_pair(c.clone(), -1);

        assert!(f.is_invertible());
        assert_eq!(f.inv(), Some(f.clone()));

        let f = Mor::from_pair(c.clone(), 2);
        assert_eq!(f.is_invertible(), false);
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn empty() { 
        let c = TngComplex::new(&0, &0, (0, 0));

        assert_eq!(c.len(), 1);
        assert_eq!(c.rank(0), 1);
    }

    #[test]
    fn single_x() { 
        let mut c = TngComplex::new(&0, &0, (0, 0));
        let x = Crossing::from_pd_code([0,1,2,3]);
        c.append(&x);

        assert_eq!(c.len(), 2);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 1);
    }

    #[test]
    fn single_x_resolved() { 
        let mut c = TngComplex::new(&0, &0, (0, 0));
        let x = Crossing::from_pd_code([0,1,2,3]).resolved(Resolution::Res0);
        c.append(&x);

        assert_eq!(c.len(), 1);
        assert_eq!(c.rank(0), 1);
    }

    #[test]
    fn two_x_disj() { 
        let mut c = TngComplex::new(&0, &0, (0, 0));
        let x0 = Crossing::from_pd_code([0,1,2,3]);
        let x1 = Crossing::from_pd_code([4,5,6,7]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.len(), 3);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);
    }

    #[test]
    fn two_x() { 
        let mut c = TngComplex::new(&0, &0, (0, 0));
        let x0 = Crossing::from_pd_code([0,4,1,5]);
        let x1 = Crossing::from_pd_code([3,1,4,2]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.len(), 3);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);
    }

    #[test]
    fn deloop_one() { 
        let mut c = TngComplex::new(&0, &0, (0, 0));
        let x0 = Crossing::from_pd_code([4,2,5,1]);
        let x1 = Crossing::from_pd_code([3,6,4,1]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.len(), 3);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);

        let e = c.find_loop(None);
        assert!(e.is_some());

        let Some((k, r, _)) = e else { panic!() };

        assert_eq!(k, KhEnhState::new(
            State::from_iter([1,0]), 
            KhLabel::from_iter([])
        ));
        assert_eq!(r, 2);

        c.deloop(&k, r, false);

        assert_eq!(c.len(), 3);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 3); // delooped here
        assert_eq!(c.rank(2), 1);
    }
}
