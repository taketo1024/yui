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

use log::{debug, info};
use rustc_hash::{FxHashMap, FxHashSet};
use auto_impl_ops::auto_ops;
use itertools::{Itertools, iproduct};
use num_traits::Zero;
use rayon::prelude::*;
use yui_core::{CloneAnd, Ring, RingOps, Sign};
use yui_core::lc::LcKey;
use yui_homology::{ChainComplex1, Summand, GrMod1};
use yui_matrix::sparse::SpMat;
use yui_link::{Edge, Node, Path, State};
use yui_core::bitseq::Bit;

use crate::kh::{KhAlgGen, KhGen, KhTensor};
use crate::util::log_progress;
use super::cob::{Cob, Dot, End, CobComp, LcCob, LcCobTrait};
use super::tng::{Tng, TngComp};

// Pacing of the merge-loop progress lines — only bites on the huge central-degree slices
// (small merges stay silent), so it costs nothing on normal builds. Edge progress runs over
// key pairs, each expanding to many attempted edges — hence the smaller step.
const MERGE_LOG_STEP: usize = 1_000_000;
const MERGE_PAIR_LOG_STEP: usize = 50_000;

// Raw pointer made shareable across rayon tasks. SAFETY is the caller's: used only in
// `eliminate_par`'s parallel value-write, where each task holds a pointer to a *distinct* vertex.
struct SendPtr<T>(*mut T);
unsafe impl<T> Send for SendPtr<T> {}

/// Which invertible edges of a vertex count as its pivot candidates (`elim_cost`, `try_eliminate_at`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ElimDir {
    Incoming,
    Outgoing,
    Both,
}

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
    in_edges: FxHashSet<TngComplexKey>,
    out_edges: FxHashMap<TngComplexKey, LcCob<R>>
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
        let in_edges = FxHashSet::default();
        let out_edges = FxHashMap::default();
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
    nodes: Vec<Node>,
    vertices: FxHashMap<TngComplexKey, TngComplexVertex<R>>,
}

impl<R> TngComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub(crate) fn new(h: &R, t: &R, deg_shift: (isize, isize), base_pt: Option<Edge>, dim: usize, vertices: FxHashMap<TngComplexKey, TngComplexVertex<R>>) -> Self {
        let ht = (h.clone(), t.clone());
        TngComplex{ ht, deg_shift, base_pt, dim, nodes: vec![], vertices }
    }

    pub fn init(h: &R, t: &R, deg_shift: (isize, isize), base_pt: Option<Edge>) -> Self { 
        let mut vertices = FxHashMap::default();
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

            let mut c = Self::new(h, t, (0, 0), base_pt, 0, FxHashMap::default());
            c.add_vertex(k, v);
            c.nodes = vec![x.clone()];
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

            let mut c = Self::new(h, t, (0, 0), base_pt, 1, FxHashMap::default());
            c.add_vertex(k0, TngComplexVertex::from(t0));
            c.add_vertex(k1, TngComplexVertex::from(t1));
            c.add_edge(&k0, &k1, sdl);
            c.nodes = vec![x.clone()];
            c
        }
    }

    pub fn from_loop(h: &R, t: &R, e: Edge, marked: bool) -> Self {
        let k = TngComplexKey::init();
        let tng = Tng::from(
            TngComp::from_path(Path::circ([e]), marked)
        );
        let v = TngComplexVertex::from(tng);

        let mut c = Self::new(h, t, (0, 0), None, 0, FxHashMap::default());
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

    pub fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    pub fn h_range(&self) -> RangeInclusive<isize> {
        let i0 = self.deg_shift.0;
        let n = self.dim() as isize;
        i0 ..= i0 + n
    }

    pub fn rank(&self, i: isize) -> usize { 
        self.keys_of_deg(i).count()
    }

    /// Whether the diagram is closed — every component is a circle (no open arc-ends).
    pub fn is_closed(&self) -> bool {
        self.iter_verts().next().is_none_or(|(_, v)| v.tng().is_closed())
    }

    /// Endpoints of the current tangle boundary.
    pub fn boundary_ends(&self) -> impl Iterator<Item = Edge> + '_ {
        self.iter_verts()
            .next()
            .into_iter()
            .flat_map(|(_, v)| v.tng().comps())
            .filter_map(|c| c.end_pts())
            .flat_map(|(e0, e1)| [e0, e1])
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

    pub fn n_verts(&self) -> usize {
        self.vertices.len()
    }

    pub fn n_edges(&self) -> usize {
        self.vertices.values().map(|v| v.out_edges.len()).sum()
    }

    pub fn vertex(&self, v: &TngComplexKey) -> &TngComplexVertex<R> { 
        &self.vertices[v]
    }

    pub fn iter_verts(&self) -> impl Iterator<Item = (&TngComplexKey, &TngComplexVertex<R>)> {
        self.vertices.iter()
    }

    pub fn add_vertex(&mut self, k: TngComplexKey, v: TngComplexVertex<R>) {
        debug_assert!(!self.contains_key(&k));
        self.vertices.insert(k, v);
    }

    pub fn remove_vertex(&mut self, k: &TngComplexKey) -> TngComplexVertex<R> {
        debug_assert!(self.contains_key(k));

        let v = self.vertices.remove(k).unwrap();

        for j in v.in_edges.iter() {
            self.vertices.get_mut(&j).unwrap().out_edges.remove(k);
        }
        
        for l in v.out_edges.keys() {
            self.vertices.get_mut(&l).unwrap().in_edges.remove(k);
        }

        v
    }

    pub fn remove_vertices(&mut self, ks: &[TngComplexKey]) {
        for k in ks {
            self.remove_vertex(k);
        }
    }

    pub(crate) fn clear_verts(&mut self) {
        self.vertices.clear();
    }
    
    fn rename_vertex_key(&mut self, k_old: &TngComplexKey, k_new: TngComplexKey) { 
        debug_assert_ne!(k_old, &k_new);
        debug_assert!(self.contains_key(k_old));
        debug_assert!(!self.contains_key(&k_new));

        let v = self.vertices.remove(k_old).unwrap();

        for j in v.in_edges.iter() {
            let u = self.vertices.get_mut(j).unwrap();
            let f = u.out_edges.remove(k_old).unwrap();
            u.out_edges.insert(k_new, f);
        }

        for l in v.out_edges.keys() {
            let w = self.vertices.get_mut(l).unwrap();
            w.in_edges.remove(k_old);
            w.in_edges.insert(k_new);
        }

        self.add_vertex(k_new, v);
    }

    fn duplicate_vertex(&mut self, k: &TngComplexKey, k_new: TngComplexKey) { 
        debug_assert_ne!(k, &k_new);
        debug_assert!(self.contains_key(k));
        debug_assert!(!self.contains_key(&k_new));

        let v_new = self.vertex(k).clone();

        for j in v_new.in_edges.iter() { 
            let u = self.vertices.get_mut(j).unwrap();
            let f = u.out_edges[k].clone();
            u.out_edges.insert(k_new, f);
        }

        for l in v_new.out_edges.keys() {
            let w = self.vertices.get_mut(l).unwrap();
            w.in_edges.insert(k_new);
        }

        self.add_vertex(k_new, v_new);
    }
    
    pub fn edge(&self, k: &TngComplexKey, l: &TngComplexKey) -> &LcCob<R> {
        self.edge_opt(k, l).unwrap()
    }

    pub fn edge_opt(&self, k: &TngComplexKey, l: &TngComplexKey) -> Option<&LcCob<R>> {
        self.vertices.get(k)?.out_edges.get(l)
    }

    // private: raw `&mut` to a cob would let callers break the edge-set consistency.
    fn edge_mut(&mut self, k: &TngComplexKey, l: &TngComplexKey) -> Option<&mut LcCob<R>> {
        self.vertices.get_mut(k)?.out_edges.get_mut(l)
    }
    
    pub fn has_edge(&self, k: &TngComplexKey, l: &TngComplexKey) -> bool {
        self.vertices[k].out_edges.contains_key(l) &&
        self.vertices[l].in_edges.contains(k) &&
        self.vertices[k].out_edges.contains_key(l)
    }

    /// Fill cost of eliminating edge `k → l`: the Markowitz count `(out(k)-1)·(in(l)-1)` — the number
    /// of parallel edges the Schur complement corrects. O(1); the nterm-weighted variant priced the
    /// true fill more accurately but was too costly to compute over dense degrees (`collect_keys` wall).
    pub(crate) fn edge_weight(&self, k: &TngComplexKey, l: &TngComplexKey) -> usize {
        (self.vertices[k].out_edges.len() - 1) * (self.vertices[l].in_edges.len() - 1)
    }

    pub fn add_edge(&mut self, k: &TngComplexKey, l: &TngComplexKey, f: LcCob<R>) { 
        debug_assert!(!self.has_edge(k, l));
        debug_assert!(!f.is_zero());

        let v = self.vertices.get_mut(k).unwrap();
        v.out_edges.insert(*l, f);

        let w = self.vertices.get_mut(l).unwrap();
        w.in_edges.insert(*k);
    }

    pub fn remove_edge(&mut self, k: &TngComplexKey, l: &TngComplexKey) -> LcCob<R> { 
        debug_assert!(self.has_edge(k, l));
        
        let w = self.vertices.get_mut(l).unwrap();
        w.in_edges.remove(k);

        let v = self.vertices.get_mut(k).unwrap();
        v.out_edges.remove(l).unwrap()
    }

    pub fn append_node(&mut self, x: &Node) {
        let (h, t) = self.ht();
        let c = Self::from_node(h, t, x, self.base_pt);
        self.merge(c);
    }

    /// Bar-Natan tensor product with `other` — see [BN05, §5](https://arxiv.org/abs/math/0410495).
    /// Raw merge (no deloop / eliminate), degree by degree.
    pub fn merge(&mut self, other: TngComplex<R>) {
        let (left, right) = self.prepare_merge(other);
        for i in self.h_range() {
            debug!("build C[{i}]...");
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
        new.nodes = self.nodes.iter().chain(other.nodes.iter()).cloned().collect();

        let left = std::mem::replace(self, new);
        (left, other)
    }

    pub(crate) fn merge_vertices(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, i: isize) {
        let total = Self::count_keys(left, right, i);
        debug!("  C[{i}]: build {total} verts ...");

        let mut n = 0;
        for (k, l) in Self::collect_keys(left, right, i) {
            let v = left.vertex(k);
            let w = right.vertex(l);
            let kl = k + l;

            let vw = TngComplexVertex::from(v.tng.connect(&w.tng));

            self.add_vertex(kl, vw);
            n += 1;
            log_progress(n, n - 1, total, MERGE_LOG_STEP);
        }

        debug!("  C[{i}]: built {n} verts.");
    }

    pub(crate) fn merge_edges(&mut self, left: &TngComplex<R>, right: &TngComplex<R>, i: isize) {
        let (h, t) = self.ht().clone();

        // progress runs over the key pairs (created-edge counts alone have no meaningful total).
        let total = Self::count_keys(left, right, i);
        debug!("  C[{i}]: build edges over {total} pairs ...");

        let mut seen = 0;
        let mut n = 0;

        for (k_l, k_r) in Self::collect_keys(left, right, i) {
            seen += 1;
            log_progress(seen, seen - 1, total, MERGE_PAIR_LOG_STEP);

            let k = k_l + k_r;
            if !self.contains_key(&k) { continue }

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
                    n += 1;
                }
            }
        }
        
        debug!("  C[{i}]: built {n} edges");
    }

    // `collect_keys` pair count at degree `i` — one flat pass per side for the degree histograms,
    // then their convolution (no walk over the product).
    fn count_keys(left: &TngComplex<R>, right: &TngComplex<R>, i: isize) -> usize {
        let count_by_deg = |c: &TngComplex<R>| {
            c.vertices.keys().fold(FxHashMap::<isize, usize>::default(), |mut m, k| {
                *m.entry((k.weight() as isize) + c.deg_shift.0).or_insert(0) += 1;
                m
            })
        };
        let l = count_by_deg(left);
        let r = count_by_deg(right);
        l.iter().map(|(d, n)| n * r.get(&(i - d)).unwrap_or(&0)).sum()
    }

    pub(crate) fn collect_keys<'a, 'b>(left: &'a TngComplex<R>, right: &'b TngComplex<R>, i: isize) -> impl Iterator<Item = (&'a TngComplexKey, &'b TngComplexKey)> {
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

    /// Deloop the `r`-th circle of vertex `k` (generalized to any `(h, t)` —
    /// [BN07, Lemma 4.1] handles only `h = 0`).
    ///
    /// A loop component `◯` is isomorphic in `Cob_{/l}` to two empty objects,
    /// one labeled with dot `X` and one with dot `1` (the dual of `Y = X − h`):
    ///
    /// ```text
    ///     ◯  ≅  ∅_X  ⊕  ∅_1     (with appropriate grading shifts)
    /// ```
    ///
    /// The two summands are inserted/projected by:
    ///   - `∅_X`: include = cup with dot `X`, project = cap with no dot
    ///   - `∅_1`: include = cup with no dot, project = cap with dot `Y`
    ///
    /// Orthogonality follows from `ε(X) = ε(Y) = 1`, `ε(1) = 0`, `XY = t·1`.
    /// For the base-pointed (reduced) variant only the `X` summand is kept.
    pub fn deloop(&mut self, k: &TngComplexKey, c: &TngComp) -> Vec<TngComplexKey> {
        debug_assert!(c.is_circle());

        #[allow(non_snake_case)]
        let updated_keys = if c.is_marked() {
            let k_X = k + KhAlgGen::X;

            self.rename_vertex_key(k, k_X);
            self.deloop_with(&k_X, c, Some(Dot::X), None);

            vec![k_X]
        } else {
            let k_X = k + KhAlgGen::X;
            let k_1 = k + KhAlgGen::I;

            self.rename_vertex_key(k, k_X);
            self.duplicate_vertex(&k_X, k_1);

            self.deloop_with(&k_X, c, Some(Dot::X), None);
            self.deloop_with(&k_1, c, None, Some(Dot::Y));

            vec![k_X, k_1]
        };

        updated_keys
    }

    fn deloop_with(&mut self, k: &TngComplexKey, c: &TngComp, birth_dot: Option<Dot>, death_dot: Option<Dot>) {
        let (h, t) = self.ht.clone();

        // own the vertex so its in/out edge lists are walked without cloning, and the outgoing
        // cobs (which live in v) are capped in place; a cob that reduces to 0 drops its edge.
        let mut v = self.vertices.remove(k).unwrap();
        v.tng.remove(c);

        // cap incoming cobs (j → k); cobs live in the sources.
        v.in_edges.retain(|j| {
            let u = self.vertices.get_mut(j).unwrap();
            let f = std::mem::take(u.out_edges.get_mut(k).unwrap())
                .cap_off(End::Tgt, c, death_dot).reduce(&h, &t);
            let keep = !f.is_zero();
            if keep { u.out_edges.insert(*k, f); } else { u.out_edges.remove(k); }
            keep
        });

        // cup outgoing cobs (k → l); cobs owned in v, so cap in place and drop l's back-ref on 0.
        v.out_edges.retain(|l, f| {
            let g = std::mem::take(f).cap_off(End::Src, c, birth_dot).reduce(&h, &t);
            let keep = !g.is_zero();
            if keep { *f = g; } else { self.vertices.get_mut(l).unwrap().in_edges.remove(k); }
            keep
        });

        self.vertices.insert(*k, v);
    }

    // Least edge-weight (Markowitz cost) to eliminate `k`, over its invertible incident edges —
    // the pivot priority.
    pub(crate) fn elim_cost(&self, k: &TngComplexKey, dir: ElimDir) -> usize {
        let v = &self.vertices[k];
        let outs = || v.out_edges.keys().filter(|l| self.edge(k, l).is_invertible()).map(|l| self.edge_weight(k, l)).min();
        let ins = || v.in_edges.iter().filter(|j| self.edge(j, k).is_invertible()).map(|j| self.edge_weight(j, k)).min();
        let min = match dir {
            ElimDir::Outgoing => outs(),
            ElimDir::Incoming => ins(),
            ElimDir::Both => Iterator::chain(outs().into_iter(), ins()).min(),
        };
        min.unwrap_or(0)
    }

    /// Gaussian-eliminate the invertible pivot `a: k0 → k1` (an iso in `Cob_{/l}`) — [BN07, Lemma 4.2].
    /// Removes both endpoints and corrects each parallel edge `d` by `d − c·a⁻¹·b`, turning the
    /// four-term segment on the left into the right. Dispatches by Schur-block size: small → serial,
    /// large → parallel.
    ///
    /// ```text
    ///       a
    ///  v0 - - -> v1         .             .
    ///     \   / b
    ///       /         ==>
    ///     /   \ c              d - ca⁻¹b
    ///  w0 -----> w1         w0 ---------> w1
    ///       d
    /// ```
    pub fn eliminate(&mut self, k0: &TngComplexKey, k1: &TngComplexKey) {
        // the parallel path only pays off on large Schur blocks; below this its fixed overhead
        // (pre-fill + drop-zeros scans + dispatch) loses to the straightforward serial apply.
        const PARALLEL_THRESHOLD: usize = 1024;

        let v0 = self.vertices.remove(k0).unwrap();
        let v1 = self.vertices.remove(k1).unwrap();

        let a = &v0.out_edges[k1];
        let Some(ainv) = a.inv() else {
            panic!("{a} is not invertible.")
        };
        // fold the Schur minus sign into a⁻¹ so corrections are added, not subtracted.
        let neg_ainv = -ainv;

        // clear only the *non-Schur* dangling refs to the removed pair (k0's in-sources, k1's
        // out-targets); the Schur-side refs are deliberately left for ser/par — see the invariant below.
        for l0 in v0.in_edges.iter() {
            self.vertices.get_mut(l0).unwrap().out_edges.remove(k0);
        }
        for l1 in v1.out_edges.keys() {
            self.vertices.get_mut(l1).unwrap().in_edges.remove(k1);
        }

        // INVARIANT (fragile): the dispatched paths still depend on the Schur-side refs left dangling
        // here. Every in-source keeps its `b` edge (`out_edges[k1]`), which ser/par read as the
        // correction input; every out-target keeps `k0` in its `in_edges` (its `c` edge already left
        // with v0). ser/par consume and clear both — do NOT strip them here, it would drop the `b` cobs.

        // block size = #in-sources × #out-targets (both minus the pivot edge).
        if (v1.in_edges.len() - 1) * (v0.out_edges.len() - 1) < PARALLEL_THRESHOLD {
            self.eliminate_ser(k0, k1, v0, v1, &neg_ainv);
        } else {
            self.eliminate_par(k0, k1, v0, v1, &neg_ainv);
        }
    }

    // Serial elimination (small blocks): iterate the owned `v0`/`v1` directly and apply each
    // correction inline via `add_to_edge` — no `in_data`/`out_data` key-vecs.
    fn eliminate_ser(&mut self, k0: &TngComplexKey, k1: &TngComplexKey, v0: TngComplexVertex<R>, v1: TngComplexVertex<R>, neg_ainv: &LcCob<R>) {
        let (h, t) = self.ht().clone();

        // correct each parallel edge `l0 → l1` by `b·a⁻¹·c` (b: l0→k1, c: k0→l1), removing each `b`.
        for l0 in v1.in_edges.iter().filter(|&l0| l0 != k0) {
            let b = self.vertices.get_mut(l0).unwrap().out_edges.remove(k1).unwrap();
            let nainv_b = b.stack(neg_ainv);
            for (l1, c) in v0.out_edges.iter().filter(|(l1, _)| *l1 != k1) {
                let corr = nainv_b.stack(c).reduce(&h, &t);
                if !corr.is_zero() {
                    self.add_to_edge(l0, l1, corr);
                }
            }
        }

        // strip k0 from the Schur out-targets' back-refs (the dangling refs are cleared in `eliminate`).
        for l1 in v0.out_edges.keys().filter(|&l1| l1 != k1) {
            self.vertices.get_mut(l1).unwrap().in_edges.remove(k0);
        }
    }

    // Parallel elimination (large blocks): pull the Schur data into key-vecs, pre-fill a zero slot
    // for every target, then write corrections in parallel via one `*mut` per source vertex.
    fn eliminate_par(&mut self, k0: &TngComplexKey, k1: &TngComplexKey, v0: TngComplexVertex<R>, v1: TngComplexVertex<R>, neg_ainv: &LcCob<R>) {
        // out-edges k0 → l1 (cobs moved out of v0) except the pivot; strip k0 from each target's in_edges.
        let out_data: Vec<_> = v0.out_edges.into_iter().filter(|(l1, _)| l1 != k1).map(|(l1, c)| {
            self.vertices.get_mut(&l1).unwrap().in_edges.remove(k0);
            (l1, c)
        }).collect();

        // in-edges l0 → k1 (cobs pulled from the sources) except the pivot; the remove doubles as cleanup.
        let in_data: Vec<_> = v1.in_edges.iter().filter(|&l0| l0 != k0).map(|&l0| {
            let b = self.vertices.get_mut(&l0).unwrap().out_edges.remove(k1).unwrap();
            (l0, b)
        }).collect();

        let (h, t) = self.ht().clone();

        // Pre-fill a zero slot for every (l0, l1) target, so the parallel phase writes only values
        // (no structural mutation → no `in_edges` race). Per `l0`: back-refs, then `*mut`.
        let vert_ptrs: Vec<SendPtr<TngComplexVertex<R>>> = in_data.iter().map(|(l0, _)| {
            for (l1, _) in &out_data {
                self.vertices.get_mut(l1).unwrap().in_edges.insert(*l0);
            }
            let v = self.vertices.get_mut(l0).unwrap();
            for (l1, _) in &out_data {
                v.out_edges.entry(*l1).or_insert_with(LcCob::zero);
            }
            SendPtr(v as *mut _)
        }).collect();

        in_data.par_iter().zip(vert_ptrs).for_each(|((_, b), vptr)| {
            let nainv_b = b.stack(neg_ainv);
            // SAFETY: `in_data` keys are distinct, so each task holds a `*mut` to a different vertex;
            // the pre-fill created every target slot, so no map is structurally mutated here (values
            // only) and the pointers stay valid for the whole phase.
            let vert = unsafe { &mut *vptr.0 };
            for (l1, c) in &out_data {
                let corr = nainv_b.stack(c).reduce(&h, &t);
                if !corr.is_zero() {
                    *vert.out_edges.get_mut(l1).unwrap() += corr;
                }
            }
        });

        // Drop the slots that stayed/became zero (pre-filled-but-unwritten + cancellations).
        for ((l0, _), (l1, _)) in iproduct!(&in_data, &out_data) {
            if self.edge_opt(l0, l1).is_some_and(|e| e.is_zero()) {
                self.remove_edge(l0, l1);
            }
        }
    }

    // Add `val` to edge `k → l` (absent = 0), in place; drop the edge if the sum is 0.
    pub(crate) fn add_to_edge(&mut self, k: &TngComplexKey, l: &TngComplexKey, val: LcCob<R>) {
        if let Some(e) = self.edge_mut(k, l) {
            *e += val;
            if e.is_zero() {
                self.remove_edge(k, l);
            }
        } else if !val.is_zero() {
            self.add_edge(k, l, val);
        }
    }

    /// Convert to the raw Khovanov chain complex (`KhGen`-keyed). See `into_raw_complex_with`.
    pub fn into_raw_complex(self) -> ChainComplex1<KhGen, R> {
        self.into_raw_complex_with(|k, a| expanded_key(k, a).as_gen(), KhGen::rel_q_deg)
    }

    /// Convert to a raw chain complex, matrix-backed for `ChainReducer`. Any circles left on a vertex
    /// are delooped here, at the matrix level: the vertex expands into one generator per circle-label
    /// assignment (mapped to `X` by `into_gen`, ordered by `q_deg`), and each edge contributes the scalar
    /// `⟨b|f|a⟩` (the deloop pairing) as a closed cobordism. A delooped vertex gives one generator.
    pub(crate) fn into_raw_complex_with<X>(
        self,
        into_gen: impl Fn(&TngComplexKey, &KhTensor) -> X,
        q_deg: impl Fn(&X) -> isize,
    ) -> ChainComplex1<X, R>
    where X: LcKey {
        let mut c = self;
        assert!(c.is_closed(), "into_raw_complex requires a closed complex (only circles expand into generators)");
        let (h, t) = c.ht().clone();

        info!("build raw complex: {}", c.stat());

        // expanded generators (vertex, circle-labels) per degree, sorted by q-degree (descending) —
        // fixes both the summand generator order and the matrix row/column order.
        let keys: FxHashMap<isize, Vec<(TngComplexKey, KhTensor)>> = c.h_range().map(|i| {
            let expanded = c.keys_of_deg(i).flat_map(|k| {
                let circles = circles_of(c.vertex(k).tng());
                label_assignments(&circles).into_iter().map(|a| (*k, a)).collect_vec()
            }).sorted_by_key(|(k, a)|
                -q_deg(&into_gen(k, a))
            ).collect_vec();
            (i, expanded)
        }).collect();

        let summands = GrMod1::generate(c.h_range(), |i|
            Summand::from_raw_generators(keys[&i].iter().map(|(k, a)| into_gen(k, a)))
        );

        let matrices = c.h_range().map(|i| {
            let cols = &keys[&i];
            let rows: FxHashMap<TngComplexKey, usize> = keys.get(&(i + 1))
                .map(|ks| ks.iter().enumerate().map(|(r, (k, a))| (expanded_key(k, a), r)).collect())
                .unwrap_or_default();

            let mut entries = vec![];
            for (j, (k, a)) in cols.iter().enumerate() {
                let src_circles = circles_of(c.vertex(k).tng());
                for l in c.vertex(k).out_edges() {
                    let fa = cap_circles(c.edge(k, l).clone(), End::Src, &src_circles, a, &h, &t);
                    if fa.is_zero() {
                        continue;
                    }
                    let tgt_circles = circles_of(c.vertex(l).tng());
                    for b in label_assignments(&tgt_circles) {
                        let g = cap_circles(fa.clone(), End::Tgt, &tgt_circles, &b, &h, &t);
                        entries.push((rows[&expanded_key(l, &b)], j, g.eval(&h, &t)));
                    }
                }
            }

            debug!("  raw d[{i}]: {} -> {}, nnz: {}", cols.len(), rows.len(), entries.len());
            let m = SpMat::from_entries((rows.len(), cols.len()), entries);

            // degree i is now fully consumed (target of d[i-1], source of d[i]); raw-drop its
            // vertices + their cobordism edges so peak memory holds at most two adjacent degrees.
            for (k, _) in cols {
                c.vertices.remove(k);
            }

            (i, m)
        }).collect_vec();

        info!("raw complex done: {} gens", keys.values().map(|ks| ks.len()).sum::<usize>());
        ChainComplex1::new_with_d_matrices(summands, 1, matrices)
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

    #[cfg(any(test, debug_assertions))]
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
        format!("(n: {}, v: {}, e: {}, w: {})", self.dim(), self.n_verts(), self.n_edges(), self.cob_weight())
    }

    /// Total stored cobordism weight: number of cob monomials summed over all edges.
    fn cob_weight(&self) -> usize {
        self.vertices.values()
            .flat_map(|v| v.out_edges.values())
            .map(|f| f.nterms())
            .sum()
    }

    pub(crate) fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge {
        let (h, t) = self.ht();
        let base_pt = self.base_pt.map(&f);
        let vertices = self.iter_verts().map(|(k, v)|
            (*k, v.convert_edges(&f))
        ).collect();
        let mut c = Self::new(h, t, self.deg_shift, base_pt, self.dim, vertices);
        c.nodes = self.nodes.iter().map(|n| n.convert_edges(&f)).collect();
        c
    }
}

// ---- deloop-via-labels: shared by TngComplex::into_raw_complex and ConeBuilder::into_raw_complex ----
// A vertex's remaining circles expand into 2^circles generators at conversion time (delooping at the
// matrix level, never on the cobordisms). A marked (based) circle is fixed to the `X` label.

pub(crate) fn circles_of(tng: &Tng) -> Vec<TngComp> {
    debug_assert!(tng.comps().all(|c| c.is_circle()));
    let (marked, unmarked): (Vec<_>, Vec<_>) = tng.comps().cloned().partition(|c| c.is_marked());
    unmarked.into_iter().chain(marked).collect()
}

pub(crate) fn label_assignments(circles: &[TngComp]) -> Vec<KhTensor> {
    KhTensor::generate(circles.len()).filter(|a|
        circles.iter().enumerate().all(|(i, c)| !c.is_marked() || a[i].is_X())
    ).collect()
}

pub(crate) fn expanded_key(k: &TngComplexKey, a: &KhTensor) -> TngComplexKey {
    let mut kk = *k;
    kk.label.append(*a);
    kk
}

// Cap the circles by the deloop pairing: a source circle labeled `X` is cupped with `Dot::X`
// (`1` plain); a target circle labeled `X` is capped plain (`1` with `Dot::Y`).
pub(crate) fn cap_circles<R>(f: LcCob<R>, e: End, circles: &[TngComp], labels: &KhTensor, h: &R, t: &R) -> LcCob<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    circles.iter().enumerate().fold(f, |f, (i, c)| {
        let dot = match (e, labels[i].is_X()) {
            (End::Src, true)  => Some(Dot::X),
            (End::Src, false) => None,
            (End::Tgt, true)  => None,
            (End::Tgt, false) => Some(Dot::Y),
        };
        f.cap_off(e, c, dot).reduce(h, t)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kh::KhTensor;
    use yui_link::Link;

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
    fn nodes_track_append_order() {
        let l = Link::test_data("3_1");
        let mut c = TngComplex::init(&0, &0, (0, 0), None);
        for x in l.nodes() {
            c.append_node(x);
        }

        assert_eq!(c.nodes(), l.nodes().cloned().collect_vec());
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
        let circ = *c.vertex(&k).tng().comp(0);
        let updated = c.deloop(&k, &circ);

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
        let circ = *c.vertex(&k).tng().comp(2);
        assert!(circ.is_circle());

        let updated = c.deloop(&k, &circ);

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
        let circ = *c.vertex(&k).tng().comp(0);
        let updated = c.deloop(&k, &circ);

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
