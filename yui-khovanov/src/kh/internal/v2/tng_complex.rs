use std::fmt::Display;
use std::ops::{Add, AddAssign, RangeInclusive};
use std::sync::RwLock;

use ahash::{AHashMap, AHashSet};
use auto_impl_ops::auto_ops;
use itertools::Itertools;
use num_traits::Zero;
use rayon::prelude::*;
use cartesian::cartesian;
use yui::{Ring, RingOps, Sign};
use yui_homology::{ChainComplex, Summand, Grid1};
use yui_link::{Crossing, Edge, State};
use yui::bitseq::Bit;

use crate::kh::{KhGen, KhAlg, KhChain, KhComplex, KhChainGen, KhTensor};
use super::cob::{Cob, Dot, Bottom, CobComp, LcCob, LcCobTrait};
use super::tng::{Tng, TngComp};

#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct TngKey { 
    pub state: State,
    pub label: KhTensor
}

impl TngKey { 
    pub(crate) fn init() -> Self { 
        Self { state: State::empty(), label: KhTensor::empty() }
    }

    pub fn weight(&self) -> usize { 
        self.state.weight()
    }

    fn append(&mut self, other: TngKey) { 
        self.state.append(other.state);
        self.label.append(other.label);
    }

    pub fn as_gen(&self, deg_shift: (isize, isize)) -> KhChainGen { 
        KhChainGen::new(self.state, self.label, deg_shift)
    }
}

#[auto_ops]
impl<'a> Add for &'a TngKey {
    type Output = TngKey;
    fn add(self, rhs: Self) -> Self::Output {
        let mut res = *self;
        res.append(*rhs);
        res
    }
}

#[auto_ops]
impl<'a> Add<KhGen> for &'a TngKey {
    type Output = TngKey;
    fn add(self, rhs: KhGen) -> Self::Output {
        let mut res = *self;
        res.label.push(rhs);
        res
    }
}

impl From<&KhChainGen> for TngKey {
    fn from(x: &KhChainGen) -> Self {
        TngKey { state: x.state, label: x.tensor }
    }
}

impl Display for TngKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.as_gen((0, 0)).fmt(f)
    }
}

#[derive(Clone, Debug)]
pub struct TngVertex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    key: TngKey,
    tng: Tng,
    in_edges: AHashSet<TngKey>,
    out_edges: AHashMap<TngKey, LcCob<R>>
}

impl<R> TngVertex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn init() -> Self { 
        let key = TngKey::init();
        let tng = Tng::empty();
        let in_edges = AHashSet::new();
        let out_edges = AHashMap::new();
        Self { key, tng, in_edges, out_edges }
    }

    pub fn tng(&self) -> &Tng { 
        &self.tng
    }

    pub fn in_edges(&self) -> impl Iterator<Item = &TngKey> {
        self.in_edges.iter()
    }

    pub fn out_edges(&self) -> impl Iterator<Item = &TngKey> {
        self.out_edges.keys()
    }

    fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge { 
        let key = self.key;
        let tng = self.tng.convert_edges(&f);
        let in_edges = self.in_edges.clone();
        let out_edges = self.out_edges.iter().map(|(k, cob)|
            (k.clone(), cob.convert_edges(&f))
        ).collect();
        TngVertex { key, tng, in_edges, out_edges }
    }
}

impl<R> Display for TngVertex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}; {})", self.key, self.tng)
    }
}

#[derive(Debug, Default)]
pub struct TngComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    ht: (R, R),
    deg_shift: (isize, isize),
    base_pt: Option<Edge>,
    vertices: AHashMap<TngKey, TngVertex<R>>,
    crossings: Vec<Crossing>,
}

impl<R> TngComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn new(h: &R, t: &R, deg_shift: (isize, isize), base_pt: Option<Edge>, vertices: AHashMap<TngKey, TngVertex<R>>, crossings: Vec<Crossing>) -> Self { 
        let ht = (h.clone(), t.clone());
        TngComplex{ ht, deg_shift, base_pt, vertices, crossings }
    }

    pub fn init(h: &R, t: &R, deg_shift: (isize, isize), base_pt: Option<Edge>) -> Self { 
        let mut vertices = AHashMap::new();
        let k0 = TngKey::init();
        let v0 = TngVertex::init();
        vertices.insert(k0, v0);

        TngComplex::new(h, t, deg_shift, base_pt, vertices, vec![])
    }

    pub fn ht(&self) -> &(R, R) { 
        &self.ht
    }

    pub fn deg_shift(&self) -> (isize, isize) { 
        self.deg_shift
    }

    pub fn set_deg_shift(&mut self, deg_shift: (isize, isize)) {
        self.deg_shift = deg_shift
    }

    pub fn base_pt(&self) -> Option<Edge> { 
        self.base_pt
    }

    pub fn contains_base_pt(&self, c: &TngComp) -> bool { 
        self.base_pt.map(|e| c.contains(e)).unwrap_or(false)
    }

    pub fn dim(&self) -> usize { 
        self.crossings.iter().filter(|x| !x.is_resolved()).count()
    }
    
    pub fn h_range(&self) -> RangeInclusive<isize> { 
        let i0 = self.deg_shift.0;
        let n = self.dim() as isize;
        i0 ..= i0 + n
    }

    pub fn rank(&self, i: isize) -> usize { 
        self.keys_of(i).count()
    }

    pub fn crossing(&self, i: usize) -> &Crossing {
        &self.crossings[i]
    }

    pub fn crossings(&self) -> &[Crossing] {
        &self.crossings
    }

    pub fn vertex(&self, v: &TngKey) -> &TngVertex<R> { 
        &self.vertices[v]
    }

    pub fn nverts(&self) -> usize { 
        self.vertices.len()
    }

    pub fn iter_verts(&self) -> impl Iterator<Item = (&TngKey, &TngVertex<R>)> {
        self.vertices.iter()
    }

    pub fn contains_key(&self, key: &TngKey) -> bool { 
        self.vertices.contains_key(key)
    }

    pub fn keys(&self) -> impl Iterator<Item = &TngKey> { 
        self.vertices.keys()
    }

    pub fn keys_of(&self, i: isize) -> impl Iterator<Item = &TngKey> { 
        let i0 = self.deg_shift.0;
        self.vertices.keys().filter(move |k| 
            (k.weight() as isize) + i0 == i
        )
    }

    pub fn keys_into(&self, k: &TngKey) -> impl Iterator<Item = &TngKey> { 
        self.vertex(k).in_edges()
    }

    pub fn keys_out_from(&self, k: &TngKey) -> impl Iterator<Item = &TngKey> { 
        self.vertex(k).out_edges()
    }

    pub fn add_vertex(&mut self, v: TngVertex<R>) { 
        assert!(!self.contains_key(&v.key));
        self.vertices.insert(v.key, v);
    }

    pub fn remove_vertex(&mut self, k: &TngKey) -> TngVertex<R> { 
        assert!(self.contains_key(k));

        let in_edges = self.keys_into(k).cloned().collect_vec();
        let out_edges = self.keys_out_from(k).cloned().collect_vec();

        let v = self.vertices.remove(k).unwrap();

        for j in in_edges { 
            self.vertices.get_mut(&j).unwrap().out_edges.remove(k);
        }
        
        for l in out_edges { 
            self.vertices.get_mut(&l).unwrap().in_edges.remove(k);
        }

        v
    }
    
    fn rename_vertex_key(&mut self, k_old: &TngKey, k_new: TngKey) { 
        assert_ne!(k_old, &k_new);

        let in_edges = self.keys_into(k_old).cloned().collect_vec(); 
        let in_removed = in_edges.into_iter().map(|j| {
            let f = self.remove_edge(&j, k_old);
            (j, f)
        }).collect_vec();

        let out_edges = self.keys_out_from(k_old).cloned().collect_vec();
        let out_removed = out_edges.into_iter().map(|l| {
            let f = self.remove_edge(k_old, &l);
            (l, f)
        }).collect_vec();

        let mut v = self.remove_vertex(k_old);
        v.key = k_new;
        self.add_vertex(v);

        for (j, f) in in_removed { 
            self.add_edge(&j, &k_new, f);
        }

        for (l, f) in out_removed { 
            self.add_edge(&k_new, &l, f);
        }
    }

    fn duplicate_vertex(&mut self, k: &TngKey, k_new: TngKey) { 
        assert_ne!(k, &k_new);

        let in_edges = self.vertex(k).in_edges.clone(); 
        let out_edges = self.vertex(k).out_edges.keys().cloned().collect_vec();

        let mut v_new = self.vertex(k).clone();
        v_new.key = k_new;
        v_new.in_edges.clear();
        v_new.out_edges.clear();

        self.add_vertex(v_new);

        for j in in_edges { 
            let f = self.edge(&j, k).clone();
            self.add_edge(&j, &k_new, f);
        }

        for l in out_edges { 
            let f = self.edge(k, &l).clone();
            self.add_edge(&k_new, &l, f);
        }
    }
    
    pub fn edge(&self, k: &TngKey, l: &TngKey) -> &LcCob<R> {
        &self.vertices[k].out_edges[l]
    }
    
    pub fn has_edge(&self, k: &TngKey, l: &TngKey) -> bool { 
        debug_assert_eq!(
            self.vertices[k].out_edges.contains_key(l), 
            self.vertices[l].in_edges.contains(k)
        );
        self.vertices[k].out_edges.contains_key(l)
    }

    fn add_edge(&mut self, k: &TngKey, l: &TngKey, f: LcCob<R>) { 
        assert!(!self.has_edge(k, l));
        assert!(!f.is_zero());

        let v = self.vertices.get_mut(k).unwrap();
        v.out_edges.insert(*l, f);

        let w = self.vertices.get_mut(l).unwrap();
        w.in_edges.insert(*k);
    }

    fn remove_edge(&mut self, k: &TngKey, l: &TngKey) -> LcCob<R> { 
        assert!(self.has_edge(k, l));
        let w = self.vertices.get_mut(l).unwrap();
        w.in_edges.remove(k);

        let v = self.vertices.get_mut(k).unwrap();
        v.out_edges.remove(l).unwrap()
    }

    fn modify_edge<F>(&mut self, k: &TngKey, l: &TngKey, map: F)
    where F: Fn(LcCob<R>) -> LcCob<R> {
        assert!(self.has_edge(k, l));

        let f = self.remove_edge(k, l);
        let map_f = map(f);

        if !map_f.is_zero() { 
            self.add_edge(k, l, map_f);
        }
    }

    pub fn append(&mut self, x: &Crossing) {
        let c = self.make_x(x);
        self.connect(c);
    }

    pub(crate) fn make_x(&self, x: &Crossing) -> Self { 
        let (h, t) = self.ht();
        let mut c = Self::new(h, t, (0, 0), None, AHashMap::new(), vec![]);

        if x.is_resolved() { 
            let mut v = TngVertex::init();
            v.tng = Tng::from_resolved(x);
            c.add_vertex(v);
        } else { 
            let mut v0 = TngVertex::init();
            v0.key.state.push_0();
            v0.tng = Tng::from_resolved(&x.resolved(Bit::Bit0));
            let k0 = v0.key;

            let mut v1 = TngVertex::init();
            v1.key.state.push_1();
            v1.tng = Tng::from_resolved(&x.resolved(Bit::Bit1));
            let k1 = v1.key;

            c.add_vertex(v0);
            c.add_vertex(v1);

            let sdl = LcCob::from(Cob::from(CobComp::sdl_from(&x)));
            c.add_edge(&k0, &k1, sdl);

            c.crossings.push(x.clone());
        }

        c
    }

    // See [Bar-Natan '05] Section 5.
    // https://arxiv.org/abs/math/0410495
    pub fn connect(&mut self, other: TngComplex<R>) { 
        let mut new = Self::connect_init(self, &other);
        for i in new.h_range() { 
            new.connect_vertices(&self, &other, i);
            new.connect_edges(&self, &other, i - 1);
        }
        *self = new
    }

    pub(crate) fn connect_init(&self, other: &TngComplex<R>) -> Self { 
        assert_eq!(self.ht(), other.ht());
        assert!(self.base_pt.is_none() || other.base_pt.is_none() || self.base_pt == other.base_pt);

        let (h, t) = self.ht();
        let base_pt = self.base_pt.or(other.base_pt);
        let deg_shift = (
            self.deg_shift.0 + other.deg_shift.0,
            self.deg_shift.1 + other.deg_shift.1
        );

        let mut new = TngComplex::init(h, t, deg_shift, base_pt);
        new.vertices.clear();
        new.crossings = Iterator::chain(self.crossings.iter(), other.crossings.iter()).cloned().collect();
        new
    }

    pub(crate) fn connect_vertices(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, i: isize) {
        let keys = self.collect_keys(left, right, i, false);
        keys.into_iter().for_each(|(k, l)| { 
            let v = left.vertex(k);
            let w = right.vertex(l);
            let kl = k + l;

            let mut vw = TngVertex::init();
            vw.key = kl;
            vw.tng = v.tng.connected(&w.tng); // D(v, w)

            self.add_vertex(vw);
        });
    }

    pub(crate) fn connect_edges(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, i: isize) {
        let (h, t) = self.ht().clone();
        let keys = self.collect_keys(left, right, i, true);

        let lock = RwLock::new(self);
        
        keys.into_par_iter().for_each(|(k0, l0)| { 
            let k0_l0 = k0 + l0;
            let v0 = left.vertex(&k0);
            let w0 = right.vertex(&l0);
            let i0 = (k0.state.weight() as isize) - left.deg_shift.0;

            let e1 = left.keys_out_from(&k0).map(|k1| { 
                let k1_l0 = k1 + l0;
                let f = left.edge(&k0, k1);
                let f_id = f.connected(&Cob::id(w0.tng())); // D(f, 1) 
                (k0_l0, k1_l0, f_id.part_eval(&h, &t))
            });

            let e2 = right.keys_out_from(&l0).map(|l1| { 
                let k0_l1 = k0 + l1;
                let f = right.edge(&l0, l1);
                let e = R::from_sign(Sign::from_parity(i0 as i64));
                let id_f = f.connected(&Cob::id(v0.tng())) * e; // (-1)^{deg(k0)} D(1, f) 
                (k0_l0, k0_l1, id_f.part_eval(&h, &t))
            });

            let mut this = lock.write().unwrap();

            Iterator::chain(e1, e2).for_each(|(k, l, f)| { 
                if this.contains_key(&l) && !f.is_zero() { 
                    this.add_edge(&k, &l, f);
                }
            });
        });
    }

    fn collect_keys<'a, 'b>(&mut self, left: &'a TngComplex<R>, right: &'b TngComplex<R>, i: isize, filter: bool) -> Vec<(&'a TngKey, &'b TngKey)> {
        let keys = left.h_range().flat_map(|i1| {
            let i2 = i - i1;
            left.keys_of(i1).flat_map(move |k| { 
                right.keys_of(i2).map(move |l|
                    (k, l)
                )
            })
        });
        
        if filter { 
            keys.filter(|(k, l)|
                self.contains_key(&(*k + *l))
            ).collect_vec()
        } else { 
            keys.collect_vec()
        }
    }

    pub fn deloop(&mut self, k: &TngKey, r: usize) -> Vec<TngKey> { 
        let c = self.vertex(k).tng.comp(r);
        assert!(c.is_circle());
        
        let based = self.contains_base_pt(c);

        #[allow(non_snake_case)]
        let updated_keys = if based { 
            let k_X = k + KhGen::X;

            self.rename_vertex_key(k, k_X);
            self.deloop_with(&k_X, r, Dot::X, Dot::None);

            vec![k_X]
        } else { 
            let k_X = k + KhGen::X;
            let k_1 = k + KhGen::I;

            self.rename_vertex_key(k, k_X);
            self.duplicate_vertex(&k_X, k_1);

            self.deloop_with(&k_X, r, Dot::X, Dot::None);
            self.deloop_with(&k_1, r, Dot::None, Dot::Y);

            vec![k_X, k_1]
        };

        updated_keys
    }

    fn deloop_with(&mut self, k: &TngKey, r: usize, birth_dot: Dot, death_dot: Dot) { 
        // remove circle
        let circ = self.vertices.get_mut(k).unwrap().tng.remove_at(r);

        let v_in = self.keys_into(k).cloned().collect_vec();
        let v_out = self.keys_out_from(k).cloned().collect_vec();

        // cap incoming cobs
        let (h, t) = self.ht.clone();
        for j in v_in.iter() { 
            self.modify_edge(j, k, |f|
                f.cap_off(Bottom::Tgt, &circ, death_dot).part_eval(&h, &t)
            );
        }
        
        // cup outgoing cobs
        for l in v_out.iter() { 
            self.modify_edge(k, l, |f|
                f.cap_off(Bottom::Src, &circ, birth_dot).part_eval(&h, &t)
            );
        }
    }

    //  Gaussian Elimination
    //
    //       a
    //  v0 - - -> v1         .             .
    //     \   / b
    //       /         ==>  
    //     /   \ c              d - ca⁻¹b
    //  w0 -----> w1         w0 ---------> w1
    //       d                

    pub fn eliminate(&mut self, k0: &TngKey, k1: &TngKey) {
        let a = self.edge(k0, k1);
        let Some(ainv) = a.inv() else { 
            panic!("{a} is not invertible.")
        };

        let (h, t) = self.ht();

        let keys = cartesian!(
            self.keys_into(k1).filter(|&l0| l0 != k0),
            self.keys_out_from(k0).filter(|&l1| l1 != k1)
        ).collect_vec();

        let values = keys.into_par_iter().map(|(l0, l1)|{
            let b = self.edge(l0, k1);
            let c = self.edge(k0, l1);
            let cab = (c * &ainv * b).part_eval(h, t);

            let s = if self.has_edge(l0, l1) { 
                let d = self.edge(l0, l1);
                d - cab
            } else { 
                -cab
            };

            (*l0, *l1, s)
        }).collect::<Vec<_>>();

        for (l0, l1, s) in values { 
            if self.has_edge(&l0, &l1) { 
                self.remove_edge(&l0, &l1);
            }
            if !s.is_zero() { 
                self.add_edge(&l0, &l1, s);
            }
        }

        self.remove_vertex(k0);
        self.remove_vertex(k1);
    }

    pub fn into_raw_complex(self) -> ChainComplex<KhChainGen, R> {
        assert!(self.is_completely_delooped());

        let summands = Grid1::generate(self.h_range(), |i| { 
            let gens = self.keys_of(i).map(|k|
                k.as_gen(self.deg_shift)
            ).sorted_by_key(|x|
                -x.q_deg()
            );
            Summand::from_raw_gens(gens)
        });

        let d = move |x: &KhChainGen| { 
            let (h, t) = self.ht();
            let k = TngKey::from(x);
            let v = self.vertex(&k);
            v.out_edges.iter().map(|(l, f)|
                (l.as_gen(x.deg_shift), f.eval(h, t))
            ).collect()
        };

        ChainComplex::new(summands, 1, move |_, z| { 
            z.apply(&d)
        })
    }

    pub fn into_kh_complex(self, canon_cycles: Vec<KhChain<R>>) -> KhComplex<R> { 
        assert!(self.is_completely_delooped());

        let (h, t) = self.ht();
        let str = KhAlg::new(h, t);
        let deg_shift = self.deg_shift;
        let reduced = self.base_pt.is_some();
        let inner = self.into_raw_complex();

        KhComplex::new_impl(inner, str, deg_shift, reduced, canon_cycles)
    }

    pub fn is_completely_delooped(&self) -> bool { 
        self.vertices.iter().all(|(_, v)|
            v.tng.is_empty()
        )
    }

    pub fn desc_d(&self) -> String { 
        let mut str = "".to_string();
        for i in self.h_range() { 
            str += &format!("C[{i}]: {}\n", self.rank(i));
            for (j, k) in self.keys_of(i).sorted().enumerate() { 
                let v = &self.vertices[k];
                str += &format!(" ({j}) {k}: {}", v.tng);
    
                for l in self.keys_out_from(k).sorted() { 
                    let f = self.edge(k, l);
                    str += &format!("\n  -> {l}: {f}");
                }
                str += "\n";
            }
        }
        str
    }

    pub fn print_d(&self) { 
        println!("{}", self.desc_d());
    }

    pub fn validate(&self) {
        for (k, v) in self.vertices.iter() { 
            // validate in_edges 
            for j in self.keys_into(k) {
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
            
            // validate out_edges 
            for l in self.keys_out_from(k) {
                assert!(
                    self.vertices.contains_key(l),
                    "no vertex for out-edge {k} -> {l}"
                );

                let w = self.vertex(l);

                assert!(
                    w.in_edges.contains(k),
                    "no in-edge {k} -> {l}"
                );
            }

            // validate cobordism
            for l in self.keys_out_from(k) {
                let w = self.vertex(l);
                let f = self.edge(k, l);

                assert!(!f.is_zero());

                f.iter().for_each(|(cob, _)| { 
                    assert_eq!(&cob.src(), v.tng(), "invalid source: {} for {cob}", v.tng());
                    assert_eq!(&cob.tgt(), w.tng(), "invalid target: {} for {cob}", w.tng());
                })
            }
        }
    }

    pub(crate) fn stat(&self) -> String { 
        format!("n: {}, v: {}", self.dim(), self.nverts())
    }

    pub fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge { 
        let (h, t) = self.ht();
        let base_pt = self.base_pt.map(|e| f(e));
        let crossings = self.crossings.iter().map(|x| x.convert_edges(&f)).collect();

        let vertices = self.iter_verts().map(|(k1, v1)| {
            let k2 = k1.clone();
            let v2 = v1.convert_edges(&f);
            (k2, v2)
        }).collect();

        TngComplex::new(h, t, self.deg_shift, base_pt, vertices, crossings)
    }
}

#[cfg(test)]
mod tests { 
    use super::*;
    use crate::kh::KhTensor;

    #[test]
    fn empty() { 
        let c = TngComplex::init(&0, &0, (0, 0), None);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 1);
    }

    #[test]
    fn single_x() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x = Crossing::from_pd_code([0,1,2,3]);
        c.append(&x);

        assert_eq!(c.dim(), 1);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 1);
    }

    #[test]
    fn single_x_resolved() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x = Crossing::from_pd_code([0,1,2,3]).resolved(Bit::Bit0);
        c.append(&x);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 1);
    }

    #[test]
    fn two_x_disj() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Crossing::from_pd_code([0,1,2,3]);
        let x1 = Crossing::from_pd_code([4,5,6,7]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.dim(), 2);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);
    }

    #[test]
    fn two_x() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Crossing::from_pd_code([0,4,1,5]);
        let x1 = Crossing::from_pd_code([3,1,4,2]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.dim(), 2);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);
    }

    #[test]
    fn deloop() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Crossing::from_pd_code([0, 1, 1, 0]).resolved(Bit::Bit0); // unknot
        c.append(&x0);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 1);

        let k = TngKey::init();
        let updated = c.deloop(&k, 0);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 2);

        assert_eq!(updated, vec![
            TngKey {
                state: State::empty(), 
                label: KhTensor::from_iter([KhGen::X])
            },
            TngKey {
                state: State::empty(), 
                label: KhTensor::from_iter([KhGen::I])
            }
        ]);
    }

    #[test]
    fn deloop_tangle() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Crossing::from_pd_code([4,2,5,1]);
        let x1 = Crossing::from_pd_code([3,6,4,1]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.dim(), 2);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);

        let k = TngKey {
            state: State::from([1,0]), 
            label: KhTensor::from_iter([])
        };
        let r = 2;

        assert!(c.vertex(&k).tng().comp(r).is_circle());

        let updated = c.deloop(&k, r);

        assert_eq!(c.dim(), 2);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 3); // delooped here
        assert_eq!(c.rank(2), 1);

        assert_eq!(updated, vec![
            TngKey {
                state: State::from([1,0]), 
                label: KhTensor::from_iter([KhGen::X])
            },
            TngKey {
                state: State::from([1,0]), 
                label: KhTensor::from_iter([KhGen::I])
            }
        ]);
    }

    #[test]
    fn deloop_based() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), Some(0)); // base point = 0
        let x0 = Crossing::from_pd_code([0, 1, 1, 0]).resolved(Bit::Bit0); // unknot
        c.append(&x0);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 1);

        let k = TngKey::init();
        let updated = c.deloop(&k, 0);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 1);

        assert_eq!(updated, vec![
            TngKey {
                state: State::empty(), 
                label: KhTensor::from_iter([KhGen::X])
            },
        ]);
    }

    #[test]
    fn connect() {
        let mut c0 = TngComplex::init(&0, &0, (0, 0), None);
        let mut c1 = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Crossing::from_pd_code([4,2,5,1]);
        let x1 = Crossing::from_pd_code([3,6,4,1]);

        c0.append(&x0);
        c1.append(&x1);

        c0.connect(c1);

        assert_eq!(c0.dim(), 2);
        assert_eq!(c0.rank(0), 1);
        assert_eq!(c0.rank(1), 2);
        assert_eq!(c0.rank(2), 1);

        c0.validate();
    }

    #[test]
    fn connect_trefoil() {
        let mut c0 = TngComplex::init(&0, &0, (0, 0), None);
        let mut c1 = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Crossing::from_pd_code([1,4,2,5]);
        let x1 = Crossing::from_pd_code([3,6,4,1]);
        let x2 = Crossing::from_pd_code([5,2,6,3]);

        c0.append(&x0);
        c0.append(&x1);
        c1.append(&x2);

        c0.connect(c1);

        assert_eq!(c0.dim(), 3);
        assert_eq!(c0.rank(0), 1);
        assert_eq!(c0.rank(1), 3);
        assert_eq!(c0.rank(2), 3);
        assert_eq!(c0.rank(3), 1);

        c0.validate();
    }
}
