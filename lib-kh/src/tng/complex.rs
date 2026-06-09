//! The tangle chain complex `[T]` over [`super::Cob`]: each vertex carries the
//! resolved tangle [`Tng`], each edge is a cobordism morphism, and the
//! differential squares to zero. Composing two complexes is the Bar-Natan
//! "tensor product" `D(Ω₁, Ω₂)` of BN05, §5, eq. (3). Local simplifications —
//! delooping and Gaussian elimination — come from BN07.
//!
//! References:
//! - BN05 — D. Bar-Natan, "Khovanov's homology for tangles and cobordisms",
//!   Geom. Topol. 9 (2005), 1443–1499.
//!   <https://doi.org/10.2140/gt.2005.9.1443>, <https://arxiv.org/abs/math/0410495>
//! - BN07 — D. Bar-Natan, "Fast Khovanov homology computations",
//!   J. Knot Theory Ramif. 16 (2007), 243–255.
//!   <https://doi.org/10.1142/S0218216507005294>, <https://arxiv.org/abs/math/0606318>

use std::fmt::Display;
use std::ops::{Add, AddAssign, RangeInclusive};

use ahash::{AHashMap, AHashSet};
use auto_impl_ops::auto_ops;
use itertools::Itertools;
use num_traits::Zero;
use yui_core::{CloneAnd, Ring, RingOps, Sign};
use yui_homology::{ChainComplex1, Summand, GrMod1};
use yui_link::{Edge, Node, Path, State};
use yui_core::bitseq::Bit;

use crate::kh::{KhAlgGen, KhGen, KhTensor};
use super::cob::{Cob, Dot, End, CobComp, LcCob, LcCobTrait};
use super::tng::{Tng, TngComp};

#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct TngComplexKey { 
    pub state: State,
    pub label: KhTensor
}

impl TngComplexKey { 
    pub(crate) fn init() -> Self { 
        Self { state: State::empty(), label: KhTensor::empty() }
    }

    pub fn weight(&self) -> usize { 
        self.state.weight()
    }

    fn append(&mut self, other: TngComplexKey) { 
        self.state.append(other.state);
        self.label.append(other.label);
    }

    pub fn as_gen(&self) -> KhGen {
        KhGen::new(self.state, self.label)
    }
}

#[auto_ops]
impl Add for &TngComplexKey {
    type Output = TngComplexKey;
    fn add(self, rhs: Self) -> Self::Output {
        let mut res = *self;
        res.append(*rhs);
        res
    }
}

#[auto_ops]
impl Add<KhAlgGen> for &TngComplexKey {
    type Output = TngComplexKey;
    fn add(self, rhs: KhAlgGen) -> Self::Output {
        let mut res = *self;
        res.label.push(rhs);
        res
    }
}

impl From<&KhGen> for TngComplexKey {
    fn from(x: &KhGen) -> Self {
        TngComplexKey { state: *x.state(), label: *x.tensor() }
    }
}

impl Display for TngComplexKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.as_gen().fmt(f)
    }
}

#[derive(Clone, Debug)]
pub struct TngComplexVertex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    tng: Tng,
    in_edges: AHashSet<TngComplexKey>,
    out_edges: AHashMap<TngComplexKey, LcCob<R>>
}

impl<R> TngComplexVertex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn init() -> Self { 
        let tng = Tng::empty();
        Self::from(tng)
    }

    pub fn tng(&self) -> &Tng { 
        &self.tng
    }

    pub fn in_edges(&self) -> impl Iterator<Item = &TngComplexKey> {
        self.in_edges.iter()
    }

    pub fn out_edges(&self) -> impl Iterator<Item = &TngComplexKey> {
        self.out_edges.keys()
    }

    pub(crate) fn c_weight(&self) -> usize {
        self.in_edges.len() * self.out_edges.len()
    }

    pub(crate) fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge {
        let tng = self.tng.convert_edges(&f);
        let in_edges = self.in_edges.clone();
        let out_edges = self.out_edges.iter().map(|(k, cob)|
            (*k, cob.map_ref(|c, r| (c.convert_edges(&f), r.clone())))
        ).collect();
        Self { tng, in_edges, out_edges }
    }
}

impl<R> From<Tng> for TngComplexVertex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn from(tng: Tng) -> Self {
        let in_edges = AHashSet::new();
        let out_edges = AHashMap::new();
        Self { tng, in_edges, out_edges }
    }
}

impl<R> Display for TngComplexVertex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.tng)
    }
}

#[derive(Debug, Default)]
pub struct TngComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    ht: (R, R),
    deg_shift: (isize, isize),
    base_pt: Option<Edge>,
    dim: usize, 
    vertices: AHashMap<TngComplexKey, TngComplexVertex<R>>,
}

impl<R> TngComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn new(h: &R, t: &R, deg_shift: (isize, isize), base_pt: Option<Edge>, dim: usize, vertices: AHashMap<TngComplexKey, TngComplexVertex<R>>) -> Self { 
        let ht = (h.clone(), t.clone());
        TngComplex{ ht, deg_shift, base_pt, dim, vertices }
    }

    pub fn init(h: &R, t: &R, deg_shift: (isize, isize), base_pt: Option<Edge>) -> Self { 
        let mut vertices = AHashMap::new();
        let k0 = TngComplexKey::init();
        let v0 = TngComplexVertex::init();
        vertices.insert(k0, v0);

        TngComplex::new(h, t, deg_shift, base_pt, 0, vertices)
    }

    pub fn from_node(h: &R, t: &R, x: &Node, base_pt: Option<Edge>) -> Self {
        if x.is_resolved() {
            let k = TngComplexKey::init();
            let tng = Tng::from_resolved(x, base_pt);
            let v = TngComplexVertex::from(tng);

            let mut c = Self::new(h, t, (0, 0), base_pt, 0, AHashMap::new());
            c.add_vertex(k, v);
            c
        } else {
            let k0 = TngComplexKey::init().clone_and(|k|
                k.state.push_0()
            );
            let t0 = Tng::from_resolved(&x.resolve(Bit::Bit0), base_pt);

            let k1 = TngComplexKey::init().clone_and(|k|
                k.state.push_1()
            );
            let t1 = Tng::from_resolved(&x.resolve(Bit::Bit1), base_pt);

            let sdl = LcCob::from(
                Cob::from(CobComp::plain(t0.clone(), t1.clone()))
            );

            let mut c = Self::new(h, t, (0, 0), base_pt, 1, AHashMap::new());
            c.add_vertex(k0, TngComplexVertex::from(t0));
            c.add_vertex(k1, TngComplexVertex::from(t1));
            c.add_edge(&k0, &k1, sdl);
            c
        }
    }

    pub fn from_loop(h: &R, t: &R, e: Edge, marked: bool) -> Self {
        let k = TngComplexKey::init();
        let tng = Tng::from(
            TngComp::from_path(Path::circ([e]), marked)
        );
        let v = TngComplexVertex::from(tng);

        let mut c = Self::new(h, t, (0, 0), None, 0, AHashMap::new());
        c.add_vertex(k, v);
        c
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

    pub fn dim(&self) -> usize {
        self.dim
    }
    
    pub fn h_range(&self) -> RangeInclusive<isize> { 
        let i0 = self.deg_shift.0;
        let n = self.dim() as isize;
        i0 ..= i0 + n
    }

    pub fn rank(&self, i: isize) -> usize { 
        self.keys_of_deg(i).count()
    }

    pub fn contains_key(&self, key: &TngComplexKey) -> bool { 
        self.vertices.contains_key(key)
    }

    pub fn keys(&self) -> impl Iterator<Item = &TngComplexKey> { 
        self.vertices.keys()
    }

    pub fn keys_of<F>(&self, pred: F) -> impl Iterator<Item = &TngComplexKey>
    where F: Fn(&TngComplexKey) -> bool {
        self.vertices.keys().filter(move |k| pred(k))
    }

    pub fn keys_of_deg(&self, i: isize) -> impl Iterator<Item = &TngComplexKey> {
        let i0 = self.deg_shift.0;
        self.keys_of(move |k| (k.weight() as isize) + i0 == i)
    }

    pub fn vertex(&self, v: &TngComplexKey) -> &TngComplexVertex<R> { 
        &self.vertices[v]
    }

    pub fn n_verts(&self) -> usize { 
        self.vertices.len()
    }

    pub fn iter_verts(&self) -> impl Iterator<Item = (&TngComplexKey, &TngComplexVertex<R>)> {
        self.vertices.iter()
    }

    /// Endpoints of the current tangle boundary. Same across all vertices
    /// (they share the partial diagram), so we read it off any one vertex.
    pub fn boundary_ends(&self) -> impl Iterator<Item = Edge> + '_ {
        self.iter_verts()
            .next()
            .into_iter()
            .flat_map(|(_, v)| v.tng().comps())
            .filter_map(|c| c.end_pts())
            .flat_map(|(e0, e1)| [e0, e1])
    }

    pub fn add_vertex(&mut self, k: TngComplexKey, v: TngComplexVertex<R>) {
        assert!(!self.contains_key(&k));
        self.vertices.insert(k, v);
    }

    pub fn remove_vertex(&mut self, k: &TngComplexKey) -> TngComplexVertex<R> { 
        assert!(self.contains_key(k));

        let in_edges = self.vertex(k).in_edges().cloned().collect_vec();
        let out_edges = self.vertex(k).out_edges().cloned().collect_vec();

        let v = self.vertices.remove(k).unwrap();

        for j in in_edges { 
            self.vertices.get_mut(&j).unwrap().out_edges.remove(k);
        }
        
        for l in out_edges { 
            self.vertices.get_mut(&l).unwrap().in_edges.remove(k);
        }

        v
    }

    pub(crate) fn clear_verts(&mut self) { 
        self.vertices.clear();
    }
    
    fn rename_vertex_key(&mut self, k_old: &TngComplexKey, k_new: TngComplexKey) { 
        assert_ne!(k_old, &k_new);

        let in_edges = self.vertex(k_old).in_edges().cloned().collect_vec(); 
        let in_removed = in_edges.into_iter().map(|j| {
            let f = self.remove_edge(&j, k_old);
            (j, f)
        }).collect_vec();

        let out_edges = self.vertex(k_old).out_edges().cloned().collect_vec();
        let out_removed = out_edges.into_iter().map(|l| {
            let f = self.remove_edge(k_old, &l);
            (l, f)
        }).collect_vec();

        let v = self.remove_vertex(k_old);
        self.add_vertex(k_new, v);

        for (j, f) in in_removed { 
            self.add_edge(&j, &k_new, f);
        }

        for (l, f) in out_removed { 
            self.add_edge(&k_new, &l, f);
        }
    }

    fn duplicate_vertex(&mut self, k: &TngComplexKey, k_new: TngComplexKey) { 
        assert_ne!(k, &k_new);

        let in_edges = self.vertex(k).in_edges.clone(); 
        let out_edges = self.vertex(k).out_edges.keys().cloned().collect_vec();

        let v_new = self.vertex(k).clone_and(|v| { 
            v.in_edges.clear();
            v.out_edges.clear();
        });

        self.add_vertex(k_new, v_new);

        for j in in_edges { 
            let f = self.edge(&j, k).clone();
            self.add_edge(&j, &k_new, f);
        }

        for l in out_edges { 
            let f = self.edge(k, &l).clone();
            self.add_edge(&k_new, &l, f);
        }
    }
    
    pub fn edge(&self, k: &TngComplexKey, l: &TngComplexKey) -> &LcCob<R> {
        &self.vertices[k].out_edges[l]
    }
    
    pub fn has_edge(&self, k: &TngComplexKey, l: &TngComplexKey) -> bool { 
        self.vertices[k].out_edges.contains_key(l) && 
        self.vertices[l].in_edges.contains(k) && 
        self.vertices[k].out_edges.contains_key(l)
    }

    pub fn add_edge(&mut self, k: &TngComplexKey, l: &TngComplexKey, f: LcCob<R>) { 
        assert!(!self.has_edge(k, l));
        assert!(!f.is_zero());

        let v = self.vertices.get_mut(k).unwrap();
        v.out_edges.insert(*l, f);

        let w = self.vertices.get_mut(l).unwrap();
        w.in_edges.insert(*k);
    }

    pub fn remove_edge(&mut self, k: &TngComplexKey, l: &TngComplexKey) -> LcCob<R> { 
        assert!(self.has_edge(k, l));
        
        let w = self.vertices.get_mut(l).unwrap();
        w.in_edges.remove(k);

        let v = self.vertices.get_mut(k).unwrap();
        v.out_edges.remove(l).unwrap()
    }

    pub fn replace_edge(&mut self, k: &TngComplexKey, l: &TngComplexKey, f: LcCob<R>) -> LcCob<R> { 
        assert!(self.has_edge(k, l));

        let v = self.vertices.get_mut(k).unwrap();
        v.out_edges.insert(*l, f).unwrap()
    }

    fn modify_edge<F>(&mut self, k: &TngComplexKey, l: &TngComplexKey, map: F)
    where F: Fn(LcCob<R>) -> LcCob<R> {
        assert!(self.has_edge(k, l));

        let v = self.vertices.get_mut(k).unwrap();
        let f = std::mem::take(v.out_edges.get_mut(l).unwrap());
        let map_f = map(f);

        if !map_f.is_zero() { 
            v.out_edges.insert(*l, map_f);
        } else { 
            self.remove_edge(k, l);
        }
    }

    pub fn append_node(&mut self, x: &Node) {
        let (h, t) = self.ht();
        let c = Self::from_node(h, t, x, self.base_pt);
        self.merge(c);
    }

    // See [Bar-Natan '05] Section 5.
    // https://arxiv.org/abs/math/0410495
    pub fn merge(&mut self, other: TngComplex<R>) { 
        let (left, right) = self.prepare_merge(other);

        for i in self.h_range() { 
            self.merge_vertices(&left, &right, i);
            self.merge_edges(&left, &right, i - 1);
        }
    }

    pub(crate) fn prepare_merge(&mut self, other: TngComplex<R>) -> (Self, Self) { 
        assert_eq!(self.ht(), other.ht());
        assert!(self.base_pt.is_none() || other.base_pt.is_none() || self.base_pt == other.base_pt);

        let (h, t) = self.ht();
        let base_pt = self.base_pt.or(other.base_pt);
        let deg_shift = (
            self.deg_shift.0 + other.deg_shift.0,
            self.deg_shift.1 + other.deg_shift.1
        );
        let dim = self.dim + other.dim;

        let mut new = TngComplex::init(h, t, deg_shift, base_pt);
        new.clear_verts();
        new.dim = dim;

        let left = std::mem::replace(self, new);
        (left, other)
    }

    pub(crate) fn merge_vertices(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, i: isize) {
        let keys = self.collect_keys(left, right, i).collect_vec();
        for (k, l) in keys { 
            let v = left.vertex(k);
            let w = right.vertex(l);
            let kl = k + l;

            let vw = TngComplexVertex::from(v.tng.connect(&w.tng));

            self.add_vertex(kl, vw);
        }
    }

    pub(crate) fn merge_edges(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, i: isize) {
        let (h, t) = self.ht().clone();
        let keys = self.collect_keys(left, right, i).filter(|(k, l)|
            self.contains_key(&(*k + *l))
        ).collect_vec();

        for (k_l, k_r) in keys {
            let k = k_l + k_r;
            
            let v_l = left.vertex(k_l);
            let v_r = right.vertex(k_r);

            let id_l = Cob::id(v_l.tng());
            let id_r = Cob::id(v_r.tng());

            let e1 = left.vertex(k_l).out_edges().map(|l_l| {
                let l = l_l + k_r;
                let f = left.edge(k_l, l_l).connect(&id_r); // D(f, 1)
                (l, f.reduce(&h, &t))
            });
            
            let i0 = (k_l.state.weight() as isize) - left.deg_shift.0;
            let sign = R::from_sign(Sign::from_parity(i0 as i64));

            let e2 = right.vertex(k_r).out_edges().map(|l_r| {
                let l = k_l + l_r;
                let id_f = right.edge(k_r, l_r).connect(&id_l) * &sign; // (-1)^{deg(k0)} D(1, f)
                (l, id_f.reduce(&h, &t))
            });

            for (l, f) in e1.chain(e2) {
                if !f.is_zero() {
                    self.add_edge(&k, &l, f);
                }
            }
        }
    }

    fn collect_keys<'a, 'b>(&self, left: &'a TngComplex<R>, right: &'b TngComplex<R>, i: isize) -> impl Iterator<Item = (&'a TngComplexKey, &'b TngComplexKey)> {
        left.h_range().filter_map(move |i1| {
            let i2 = i - i1;
            right.h_range().contains(&i2).then_some((i1, i2))
        }).flat_map(move |(i1, i2)|
            left.keys_of_deg(i1).flat_map(move |k1| 
                right.keys_of_deg(i2).map(move |k2|
                    (k1, k2)
            ))
        )
    }

    // Delooping (generalized to any `(h, t)` — [BN07, Lemma 4.1] handles only `h = 0`).
    //
    // A loop component `⚫` is isomorphic in `Cob_{/l}` to two empty objects,
    // one labeled with dot `X` and one with dot `1` (the dual of `Y = X − h`):
    //
    //     ⚫  ≅  ∅_X  ⊕  ∅_1     (with appropriate grading shifts)
    //
    // The two summands are inserted/projected by:
    //   - `∅_X`: include = cup with dot `X`, project = cap with no dot
    //   - `∅_1`: include = cup with no dot, project = cap with dot `Y`
    //
    // Orthogonality follows from `ε(X) = ε(Y) = 1`, `ε(1) = 0`, `XY = t·1`.
    // For the base-pointed (reduced) variant only the `X` summand is kept.
    pub fn deloop(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> {
        let c = self.vertex(k).tng.comp(r);
        assert!(c.is_circle());
        
        #[allow(non_snake_case)]
        let updated_keys = if c.is_marked() { 
            let k_X = k + KhAlgGen::X;

            self.rename_vertex_key(k, k_X);
            self.deloop_with(&k_X, r, Dot::X, Dot::None);

            vec![k_X]
        } else { 
            let k_X = k + KhAlgGen::X;
            let k_1 = k + KhAlgGen::I;

            self.rename_vertex_key(k, k_X);
            self.duplicate_vertex(&k_X, k_1);

            self.deloop_with(&k_X, r, Dot::X, Dot::None);
            self.deloop_with(&k_1, r, Dot::None, Dot::Y);

            vec![k_X, k_1]
        };

        updated_keys
    }

    fn deloop_with(&mut self, k: &TngComplexKey, r: usize, birth_dot: Dot, death_dot: Dot) { 
        // remove circle
        let circ = self.vertices.get_mut(k).unwrap().tng.remove_at(r);

        let v_in = self.vertex(k).in_edges().cloned().collect_vec();
        let v_out = self.vertex(k).out_edges().cloned().collect_vec();

        // cap incoming cobs
        let (h, t) = self.ht.clone();
        for j in v_in.iter() { 
            self.modify_edge(j, k, |f|
                f.cap_off(End::Tgt, &circ, death_dot).reduce(&h, &t)
            );
        }
        
        // cup outgoing cobs
        for l in v_out.iter() { 
            self.modify_edge(k, l, |f|
                f.cap_off(End::Src, &circ, birth_dot).reduce(&h, &t)
            );
        }
    }

    // Gaussian elimination — [BN07, Lemma 4.2].
    //
    // When `a: v0 → v1` is invertible (i.e. an iso in `Cob_{/l}`), the four-term
    // segment on the left is homotopy-equivalent to the simpler segment on the
    // right with `v0`, `v1` removed and the parallel edge `d` corrected by `−c·a⁻¹·b`:
    //
    //       a
    //  v0 - - -> v1         .             .
    //     \   / b
    //       /         ==>
    //     /   \ c              d - ca⁻¹b
    //  w0 -----> w1         w0 ---------> w1
    //       d                

    pub fn eliminate(&mut self, k0: &TngComplexKey, k1: &TngComplexKey) {
        let a = self.edge(k0, k1);
        let Some(ainv) = a.inv() else { 
            panic!("{a} is not invertible.")
        };

        let (h, t) = self.ht().clone();

        let in_keys: Vec<_> = self.vertex(k1).in_edges()
            .filter(|&l0| l0 != k0).copied().collect();
        let out_keys: Vec<_> = self.vertex(k0).out_edges()
            .filter(|&l1| l1 != k1).copied().collect();

        for l0 in &in_keys {
            let b = self.edge(l0, k1);
            let ainv_b = b.stack(&ainv);

            for l1 in &out_keys {
                let c = self.edge(k0, l1);
                let c_ainv_b = ainv_b.stack(c).reduce(&h, &t);
                if c_ainv_b.is_zero() {
                    continue
                }

                let s = if self.has_edge(l0, l1) {
                    let d = self.edge(l0, l1);
                    d - c_ainv_b
                } else { 
                    -c_ainv_b
                };

                match (self.has_edge(l0, l1), s.is_zero()) { 
                    (false, false) => self.add_edge(l0, l1, s),
                    (true,  false) => { self.replace_edge(l0, l1, s); },
                    (true,  true)  => { self.remove_edge(l0, l1); },
                    _              => ()
                }
            }
        }

        self.remove_vertex(k0);
        self.remove_vertex(k1);
    }

    pub fn into_raw_complex(self) -> ChainComplex1<KhGen, R> {
        assert!(self.is_completely_delooped());

        let summands = GrMod1::generate(self.h_range(), |i| {
            let gens = self.keys_of_deg(i).map(|k|
                k.as_gen()
            ).sorted_by_key(|x|
                -x.rel_q_deg()
            );
            Summand::from_raw_generators(gens)
        });

        let d = move |x: &KhGen| {
            let (h, t) = self.ht();
            let k = TngComplexKey::from(x);
            let v = self.vertex(&k);
            v.out_edges.iter().map(|(l, f)|
                (l.as_gen(), f.eval(h, t))
            ).collect()
        };

        ChainComplex1::new(summands, 1, move |_, z| { 
            z.apply(&d)
        })
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
            for (j, k) in self.keys_of_deg(i).sorted().enumerate() { 
                let v = &self.vertices[k];
                str += &format!(" ({j}) {k}: {}", v.tng);
    
                for l in self.vertex(k).out_edges().sorted() { 
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

    #[cfg(debug_assertions)]
    pub fn validate(&self) {
        for (k, v) in self.vertices.iter() { 
            // validate in_edges 
            for j in self.vertex(k).in_edges() {
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
            for l in self.vertex(k).out_edges() {
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
            for l in self.vertex(k).out_edges() {
                let w = self.vertex(l);
                let f = self.edge(k, l);

                assert!(!f.is_zero());

                f.iter().for_each(|(cob, _)| { 
                    assert_eq!(&cob.reconst_src(), v.tng(), "invalid source: {} for {cob}", v.tng());
                    assert_eq!(&cob.reconst_tgt(), w.tng(), "invalid target: {} for {cob}", w.tng());
                })
            }
        }
    }

    pub(crate) fn stat(&self) -> String {
        format!("(n: {}, v: {})", self.dim(), self.n_verts())
    }

    pub(crate) fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge {
        let (h, t) = self.ht();
        let base_pt = self.base_pt.map(&f);
        let vertices = self.iter_verts().map(|(k, v)|
            (*k, v.convert_edges(&f))
        ).collect();
        Self::new(h, t, self.deg_shift, base_pt, self.dim, vertices)
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
        let x = Node::from_pd_code([1,4,2,5]);
        c.append_node(&x);

        assert_eq!(c.dim(), 1);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 1);
    }

    #[test]
    fn single_x_resolved() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x = Node::from_pd_code([1,4,2,5]).resolve(Bit::Bit0);
        c.append_node(&x);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 1);
    }

    #[test]
    fn two_x_disj() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Node::from_pd_code([1,4,2,5]);
        let x1 = Node::from_pd_code([11,14,12,15]);

        c.append_node(&x0);
        c.append_node(&x1);

        assert_eq!(c.dim(), 2);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);
    }

    #[test]
    fn two_x() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Node::from_pd_code([4,2,5,1]);
        let x1 = Node::from_pd_code([3,6,4,1]);

        c.append_node(&x0);
        c.append_node(&x1);

        assert_eq!(c.dim(), 2);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);
    }

    #[test]
    fn deloop() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Node::from_pd_code([1,2,2,1]).resolve(Bit::Bit0); // unknot
        c.append_node(&x0);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 1);

        let k = TngComplexKey::init();
        let updated = c.deloop(&k, 0);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 2);

        assert_eq!(updated, vec![
            TngComplexKey {
                state: State::empty(), 
                label: KhTensor::from_iter([KhAlgGen::X])
            },
            TngComplexKey {
                state: State::empty(), 
                label: KhTensor::from_iter([KhAlgGen::I])
            }
        ]);
    }

    #[test]
    fn deloop_tangle() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Node::from_pd_code([4,2,5,1]);
        let x1 = Node::from_pd_code([3,6,4,1]);

        c.append_node(&x0);
        c.append_node(&x1);

        assert_eq!(c.dim(), 2);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);

        let k = TngComplexKey {
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
            TngComplexKey {
                state: State::from([1,0]), 
                label: KhTensor::from_iter([KhAlgGen::X])
            },
            TngComplexKey {
                state: State::from([1,0]), 
                label: KhTensor::from_iter([KhAlgGen::I])
            }
        ]);
    }

    #[test]
    fn deloop_based() { 
        let mut c = TngComplex::init(&0, &0, (0, 0), Some(1)); // base point = 1
        let x0 = Node::from_pd_code([1,2,2,1]).resolve(Bit::Bit0); // unknot
        c.append_node(&x0);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 1);

        let k = TngComplexKey::init();
        let updated = c.deloop(&k, 0);

        assert_eq!(c.dim(), 0);
        assert_eq!(c.rank(0), 1);

        assert_eq!(updated, vec![
            TngComplexKey {
                state: State::empty(), 
                label: KhTensor::from_iter([KhAlgGen::X])
            },
        ]);
    }

    #[test]
    fn merge() {
        let mut c0 = TngComplex::init(&0, &0, (0, 0), None);
        let mut c1 = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Node::from_pd_code([4,2,5,1]);
        let x1 = Node::from_pd_code([3,6,4,1]);

        c0.append_node(&x0);
        c1.append_node(&x1);

        c0.merge(c1);

        assert_eq!(c0.dim(), 2);
        assert_eq!(c0.rank(0), 1);
        assert_eq!(c0.rank(1), 2);
        assert_eq!(c0.rank(2), 1);

        c0.validate();
    }

    #[test]
    fn merge_trefoil() {
        let mut c0 = TngComplex::init(&0, &0, (0, 0), None);
        let mut c1 = TngComplex::init(&0, &0, (0, 0), None);
        let x0 = Node::from_pd_code([1,4,2,5]);
        let x1 = Node::from_pd_code([3,6,4,1]);
        let x2 = Node::from_pd_code([5,2,6,3]);

        c0.append_node(&x0);
        c0.append_node(&x1);
        c1.append_node(&x2);

        c0.merge(c1);

        assert_eq!(c0.dim(), 3);
        assert_eq!(c0.rank(0), 1);
        assert_eq!(c0.rank(1), 3);
        assert_eq!(c0.rank(2), 3);
        assert_eq!(c0.rank(3), 1);

        c0.validate();
    }
}
