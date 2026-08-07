//! Builds the involutive cone `Cone(1 + τ)` at the *cobordism* level — before delooping — for a
//! strongly invertible link. [`ConeBuilder`] owns a [`SymTngBuilder`], drives the symmetric build,
//! and converts the reduced cone at the boundary: `into_raw_complex` yields the matrix-backed
//! `KhIGen`-keyed chain complex, `eval_khi_elements` the canon classes. Char-2 only.
//!
//! The build always closes with the incremental `cone_merge`: the final chunk merge is fused with
//! cone construction + reduction per symmetric degree, so the full un-delooped product is never
//! materialized (a non-chunked build uses the degenerate two-piece plan — see `plan`).
//!
//! Reference:
//! - T. Sano, "Involutive Khovanov homology and equivariant knots",
//!   Algebr. Geom. Topol. 25 (2025), 5059–5111.

use std::ops::RangeInclusive;

use delegate::delegate;
use itertools::Itertools;
use log::{debug, info};
use num_traits::Zero;
use yui_core::bitseq::Bit;
use yui_core::abst::{Ring, RingOps};
use yui_link::{Edge, InvLink};

use rayon::prelude::*;
use rustc_hash::FxHashSet;
use yui_homology::ChainComplex1;

use crate::kh::{KhChain, KhGen};
use crate::khi::{KhIChain, KhIGen, KhIGenExt};
use crate::tng::{Cob, CobComp, End, LcCob, LcCobTrait, Tng, TngComplex, TngComplexElem, TngComplexKey, TngComplexVertex, circles_of, label_assignments, expanded_key, cap_circles};
use super::{reachable_range, SymTngBuilder, SymBuildConfig, TngComplexBuilder, BuildConfig};
use super::builder::PROGRESS_LOG_STEP;
use log::Level;
use yui_core::util::log::log_progress;

const CHUNK: usize = 4096;

// τ-orbit class of a symmetric-complex vertex: the representative (smaller key) and the τ-fixed
// survive into the direct cone; the non-representative column is never materialized.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum OrbitClass {
    Fixed, Rep, Drop
}

pub struct ConeBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    inner: SymTngBuilder<R>,
    cone: TngComplexBuilder<R>, // the reduced cone, filled in by the final `cone_merge`
    full_extend: bool, // build by the reference extension instead — see `cone_extend_full`
}

impl<R> ConeBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_inv_link(l: &InvLink, h: &R, t: &R, reduced: bool) -> Self {
        let inner = SymTngBuilder::from_inv_link(l, h, t, reduced);
        let cone = TngComplexBuilder::init(h, t, (0, 0), None); // replaced by `cone_merge`
        Self { inner, cone, full_extend: false }
    }

    pub fn with_config(mut self, config: SymBuildConfig) -> Self {
        self.inner = self.inner.with_config(config);
        self
    }

    /// Build the cone by the unreduced reference extension rather than the direct one.
    /// Slower, but both must give the same homology.
    pub fn with_full_extend(mut self, flag: bool) -> Self {
        self.full_extend = flag;
        self
    }

    // Build the symmetric complex via `SymTngBuilder::run` (chunked or incremental; its finalize
    // ends with the free-eliminate fixpoint), then build the cone over the finished complex.
    pub fn run(mut self) -> Self {
        self.inner = self.inner.run();
        self.build_cone();
        self.finalize();
        self
    }

    delegate! {
        to self.cone {
            pub fn into_tng_complex(self) -> TngComplex<R>;
            pub fn eval_elements(&self) -> Vec<KhChain<R>>;
        }
    }

    // Build the reduced cone over the FINISHED symmetric complex (Sano2026, Prop 4.6): per degree,
    // emit the symmetry-broken cone extension and collapse the invertible `1`-edges two degrees
    // behind. Nothing is merged here — the complex and its τ key map are already complete.
    fn build_cone(&mut self) {
        let range = reachable_range(self.inner.complex().h_range(), &self.inner.config().h_range, self.inner.n_nodes());
        let cone_config = cone_build_config(self.inner.config());
        self.cone = TngComplexBuilder::from_tng_complex(cone_shell(self.inner.complex()), cone_config);

        info!("build cone over {}, range: {range:?}", self.inner.complex().stat());

        self.seed_elements();

        let top = *range.end();
        for d in range {
            debug!("build cone C[{d}]...");
            if self.full_extend {
                self.cone_extend_full(d);
            } else {
                self.cone_extend_reduced(d);
                self.reduce_elements(d - 1); // degree d-1's out-edges are now complete
            }
            self.prune_consumed(d - 2); // free the consumed symmetric degree before the heavy eliminate
            self.cone.eliminate_in(d - 2); // stragglers: τ-fixed `1+τ` units, correction-created units
            debug!("  built cone C[{d}]: {}.", self.cone.complex().rank(d));
        }

        if !self.full_extend {
            self.reduce_elements(top); // the top degree has no further out-edges — Iτ pushes only
        }

        info!("cone eliminate top C[{}..={}]", top - 1, top);
        for d in (top - 1) ..= top {
            self.prune_consumed(d);
            self.cone.eliminate_in(d);
        }

        self.prune_isolated_top(top);

        info!("built cone: {}", self.cone.stat());
    }

    // Lift each completed symmetric canon cycle to its `B` (bit-0) and `Q` (bit-1) cone copies.
    fn seed_elements(&mut self) {
        let lifted = self.inner.elements().content().iter().flat_map(|e|
            [lift_elem(e, Bit::Bit0), lift_elem(e, Bit::Bit1)]
        ).collect_vec();
        self.cone.elements_mut().set(lifted);
    }

    // The full (non-reduced) cone extension — kept for debugging against `cone_extend_reduced`.
    // Adds the symmetric complex's degree-`d` slice to the cone: the two copies `k·0`, `k·1` of each
    // vertex, the horizontal edges *into* degree `d`, and the `1+τ` edges out of `k·0`. Processed
    // ascending, every edge lands exactly once (its target's degree). The free-orbit `1+τ` pivots
    // are collapsed two degrees behind by the loop's `eliminate_in`.
    fn cone_extend_full(&mut self, d: isize) {
        let keys = self.inner.complex().keys_of_deg(d).copied().collect_vec();

        debug!("  cone-extend C[{d}]: +{} verts", 2 * keys.len());

        for k in keys.iter() {
            let tng = self.inner.complex().vertex(k).tng().clone();
            self.cone.complex_mut().add_vertex(with_bit(k, Bit::Bit0), TngComplexVertex::from(tng.clone()));
            self.cone.complex_mut().add_vertex(with_bit(k, Bit::Bit1), TngComplexVertex::from(tng));
        }

        // horizontal edges into degree d (each layer copies the symmetric differential).
        for k in keys.iter() {
            for j in self.inner.complex().vertex(k).in_edges().copied().collect_vec() {
                let f = self.inner.complex().edge(&j, k).clone();
                self.cone.complex_mut().add_edge(&with_bit(&j, Bit::Bit0), &with_bit(k, Bit::Bit0), f.clone());
                self.cone.complex_mut().add_edge(&with_bit(&j, Bit::Bit1), &with_bit(k, Bit::Bit1), f);
            }
        }

        // connecting differential (1 + τ): k·0 → k·1 (id) and k·0 → τk·1 (τ-cyl).
        for k in keys.iter() {
            let tng = self.inner.complex().vertex(k).tng().clone();
            let id = LcCob::from(Cob::id(&tng));
            let tau = LcCob::from(tau_cob(&tng, |e| self.inner.inv_edge(e)));
            let tk = *self.inner.key_map().inv_key(k);
            let k0 = with_bit(k, Bit::Bit0);

            if &tk == k {
                let f = id + tau; // same target: sum over char 2
                if !f.is_zero() {
                    self.cone.complex_mut().add_edge(&k0, &with_bit(k, Bit::Bit1), f);
                }
            } else {
                self.cone.complex_mut().add_edge(&k0, &with_bit(k, Bit::Bit1), id);
                self.cone.complex_mut().add_edge(&k0, &with_bit(&tk, Bit::Bit1), tau);
            }
        }
    }

    // Symmetry-breaking reduction (Sano2026, Prop 4.6): per free τ-orbit only the
    // representative's two copies are created; the dropped column's contributions enter as
    // closed-form corrections. Corrections through two dropped corners vanish (an edge into a
    // dropped `l⁰` is an in-edge of a removed pivot source; one out of a dropped `j¹` is an
    // out-edge of a removed pivot target), so the rewrite is one elimination step deep:
    //   surv → surv : verbatim on both layers,
    //   surv → drop : (ii)  j¹ → (τk)¹ += Iτ(k)∘f,
    //   drop → surv : (iii) (τj)⁰ → k⁰ += f∘Iτ(τj),
    //   asymmetric elimination through dropped N at d−1 : (i) X¹ → k⁰ += (X→N)∘(N→k),
    // and only τ-fixed vertices keep a vertical `1 + τ` (a representative's identity cancels via τ²).
    fn cone_extend_reduced(&mut self, d: isize) {
        let keys = self.inner.complex().keys_of_deg(d).copied().collect_vec();
        self.add_cone_vertices(d, &keys);
        self.add_horizontal_edges(d, &keys);
        self.add_vertical_edges(d, &keys);
        self.add_diag_edges(d - 1);
    }

    // vertices: representatives and τ-fixed only.
    fn add_cone_vertices(&mut self, d: isize, keys: &[TngComplexKey]) {
        let dropped = keys.iter().filter(|k| self.tau_key(k).1 == OrbitClass::Drop).count();
        debug!("  C[{d}]: build {} verts ({dropped} dropped)", 2 * (keys.len() - dropped));

        for k in keys.iter() {
            if self.tau_key(k).1 == OrbitClass::Drop {
                continue;
            }
            let tng = self.inner.complex().vertex(k).tng().clone();
            self.cone.complex_mut().add_vertex(with_bit(k, Bit::Bit0), TngComplexVertex::from(tng.clone()));
            self.cone.complex_mut().add_vertex(with_bit(k, Bit::Bit1), TngComplexVertex::from(tng));
        }
    }

    // horizontal edges into degree d, rewritten per the tables in `cone_extend_reduced`. The Drop
    // arms each cost a cobordism composition (`Iτ∘f` / `f∘Iτ`), independent per edge — compute in
    // parallel (chunked over k to bound the transient), then apply the accumulating add_to_edge
    // serially (entries can collide: a Drop-arm correction may land on a genuine surviving edge).
    fn add_horizontal_edges(&mut self, d: isize, keys: &[TngComplexKey]) {
        debug!("  C[{d}]: add horizontal edges ({} keys)...", keys.len());

        let (h, t) = self.cone.complex().ht().clone();
        let inner = self.inner.complex();
        let mut done = 0;

        for keys_chunk in keys.chunks(CHUNK) {
            // reference shadows: the nested `move` closures can only capture `Copy` refs.
            let this = &*self;
            let (h, t) = (&h, &t);
            let entries: Vec<(TngComplexKey, TngComplexKey, LcCob<R>)> = keys_chunk.par_iter().flat_map_iter(|k| {
                let (tk, ck) = this.tau_key(k);
                inner.vertex(k).in_edges().flat_map(move |j| {
                    let f = inner.edge(j, k);
                    let (tj, cj) = this.tau_key(j);
                    // per the rewrite tables: at most one cone edge per layer (e0: bit-0, e1: bit-1).
                    let (e0, e1) = match (cj, ck) {
                        (OrbitClass::Drop, OrbitClass::Drop) => (None, None),
                        (OrbitClass::Drop, _) => {
                            let itau = LcCob::from(tau_cob(inner.vertex(&tj).tng(), |e| this.inner.inv_edge(e)));
                            let corr = itau.stack(f).reduce(h, t);
                            let e0 = (!corr.is_zero()).then(|| (with_bit(&tj, Bit::Bit0), with_bit(k, Bit::Bit0), corr));
                            (e0, None)
                        }
                        (_, OrbitClass::Drop) => {
                            let itau = LcCob::from(tau_cob(inner.vertex(k).tng(), |e| this.inner.inv_edge(e)));
                            let corr = f.stack(&itau).reduce(h, t);
                            let e1 = (!corr.is_zero()).then(|| (with_bit(j, Bit::Bit1), with_bit(&tk, Bit::Bit1), corr));
                            (None, e1)
                        }
                        _ => (
                            Some((with_bit(j, Bit::Bit0), with_bit(k, Bit::Bit0), f.clone())),
                            Some((with_bit(j, Bit::Bit1), with_bit(k, Bit::Bit1), f.clone())),
                        ),
                    };
                    [e0, e1].into_iter().flatten()
                })
            }).collect();

            for (src, dst, f) in entries {
                self.cone.complex_mut().add_to_edge(&src, &dst, f);
            }

            let prev = done;
            done += keys_chunk.len();
            log_progress(Level::Debug, done, prev, keys.len(), PROGRESS_LOG_STEP, 2);
        }
    }

    // verticals: τ-fixed only.
    fn add_vertical_edges(&mut self, d: isize, keys: &[TngComplexKey]) {
        let target = keys.iter().filter(|k| self.tau_key(k).1 == OrbitClass::Fixed).collect_vec();

        debug!("  C[{d}]: add vertical edges ({} keys)...", target.len());
        let mut done = 0;

        for k in target.iter() {
            let tng = self.inner.complex().vertex(k).tng().clone();
            let id = LcCob::from(Cob::id(&tng));
            let tau = LcCob::from(tau_cob(&tng, |e| self.inner.inv_edge(e)));
            let f = id + tau; // same target: sum over char 2
            if !f.is_zero() {
                self.cone.complex_mut().add_edge(&with_bit(k, Bit::Bit0), &with_bit(k, Bit::Bit1), f);
            }

            done += 1;
            log_progress(Level::Debug, done, done - 1, target.len(), PROGRESS_LOG_STEP, 2);
        }
    }

    // (i) asymmetric elimination: eliminating each dropped N's diagonal `1+τ` edge reconnects its
    // in-edges to its out-edges, X¹ → k⁰ += (X→N)∘(N→k) — one full cobordism composition per
    // (N, in-edge, out-edge) triple, O(dropped × in × out), the dominant cost. Independent per
    // triple, so compute in parallel (chunked over N to bound the transient), then apply the
    // accumulating add_to_edge serially (corrections collide with horizontal edges and each other).
    fn add_diag_edges(&mut self, d: isize) {
        let dropped = self.inner.complex().keys_of_deg(d)
            .filter(|n| self.tau_key(n).1 == OrbitClass::Drop)
            .copied().collect_vec();

        debug!("  C[{}]: add diagonal edges (dropped: {})", d, dropped.len());

        let (h, t) = self.cone.complex().ht().clone();
        let inner = self.inner.complex();
        let mut done = 0;

        for dropped_chunk in dropped.chunks(CHUNK) {
            let (h, t) = (&h, &t);
            let corrs: Vec<(TngComplexKey, TngComplexKey, LcCob<R>)> = dropped_chunk.par_iter().flat_map_iter(|n| {
                let ins = inner.vertex(n).in_edges().copied()
                    .filter(|x| self.tau_key(x).1 != OrbitClass::Drop).collect_vec();
                let outs = inner.vertex(n).out_edges().copied()
                    .filter(|k| self.tau_key(k).1 != OrbitClass::Drop).collect_vec();
                ins.into_iter().cartesian_product(outs).filter_map(move |(x, k)| {
                    let corr = inner.edge(&x, n).stack(inner.edge(n, &k)).reduce(h, t);
                    (!corr.is_zero()).then(||
                        (with_bit(&x, Bit::Bit1), with_bit(&k, Bit::Bit0), corr)
                    )
                })
            }).collect();

            for (src, dst, corr) in corrs {
                self.cone.complex_mut().add_to_edge(&src, &dst, corr);
            }

            let prev = done;
            done += dropped_chunk.len();
            log_progress(Level::Debug, done, prev, dropped.len(), PROGRESS_LOG_STEP, 2);
        }
    }

    // Retract canon-element components off the dropped column of degree `d`
    // (Sano2026, Prop 4.6 SDR): an entry at `N·0` is the pivot's dependent coordinate and drops;
    // an entry at `N·1` redirects to `(τN)·1` via Iτ and to `l·0` via each out-edge `N → l`
    // (mirroring `eliminate_from` with the identity pivot). Corrections landing on a dropped
    // `l·0` are removed by the next degree's rewrite, matching the sequential SDR composition.
    fn reduce_elements(&mut self, d: isize) {
        let n_elems = self.cone.elements().content().len();
        if n_elems == 0 {
            return;
        }

        // only keys the elements actually reference need a push table (the canon elements touch a
        // tiny fraction of the ~10⁵ dropped keys per degree).
        let referenced: FxHashSet<TngComplexKey> = self.cone.elements().content().iter()
            .flat_map(|e| e.out_cob().keys().copied())
            .map(|mut k| { k.state.remove(k.state.len() - 1); k })
            .collect();

        let dropped = self.inner.complex().keys_of_deg(d)
            .filter(|n| referenced.contains(n) && self.tau_key(n).1 == OrbitClass::Drop)
            .copied().collect_vec();

        if dropped.is_empty() {
            return;
        }

        debug!("    rewrite-elements {n_elems} in C[{d}]: {} referenced-dropped", dropped.len());

        let (h, t) = self.cone.complex().ht().clone();

        let pushes = dropped.iter().map(|n| {
            let tn = self.tau_key(n).0;
            let itau = LcCob::from(tau_cob(self.inner.complex().vertex(n).tng(), |e| self.inner.inv_edge(e)));
            let mut outs = vec![(with_bit(&tn, Bit::Bit1), itau)];
            for l in self.inner.complex().vertex(n).out_edges() {
                outs.push((with_bit(l, Bit::Bit0), self.inner.complex().edge(n, l).clone()));
            }
            (with_bit(n, Bit::Bit0), with_bit(n, Bit::Bit1), outs)
        }).collect_vec();

        let elems = self.cone.elements_mut().take().into_iter().map(|mut e| {
            e.modify_out_cob(|mut cob| {
                for (n0, n1, outs) in pushes.iter() {
                    cob.remove(n0);
                    let Some(f) = cob.remove(n1) else { continue };
                    for (target, edge) in outs.iter() {
                        let corr = f.stack(edge).reduce(&h, &t);
                        let s = if let Some(g) = cob.remove(target) { g - corr } else { -corr };
                        if !s.is_zero() {
                            cob.insert(*target, s);
                        }
                    }
                }
                cob
            });
            e
        }).collect_vec();

        self.cone.elements_mut().set(elems);
    }

    // τ-orbit classification via the inner key map; the representative is the smaller key.
    fn tau_key(&self, k: &TngComplexKey) -> (TngComplexKey, OrbitClass) {
        let tk = *self.inner.key_map().inv_key(k);
        let class = if tk == *k {
            OrbitClass::Fixed
        } else if *k < tk {
            OrbitClass::Rep
        } else {
            OrbitClass::Drop
        };
        (tk, class)
    }

    // Drop a consumed symmetric degree and its τ key-map entries — never needed again.
    fn prune_consumed(&mut self, d: isize) {
        let doomed = self.inner.complex().keys_of_deg(d).copied().collect_vec();
        self.inner.complex_mut().remove_vertices(&doomed);
        let i0 = self.inner.complex().deg_shift().0;
        self.inner.key_map_mut().drop(|k| k.weight() as isize + i0 == d);
    }

    // Restrict a truncated cone to what KhI within the window needs: drop the out-of-window layer
    // C[top+1] (the `Bit1` shift of C[top]), then prune no-in-edge vertices at C[top] — they fed
    // only the dropped layer. No-op for a full-range build (its C[top+1] is genuine).
    fn prune_isolated_top(&mut self, top: isize) {
        // prune only when the window truncates the KhI complex: if the requested top exceeds the
        // reachable Kh top, C[top+1] is the genuine KhI top degree (Bit1 of C[top]) and must stay.
        let Some(window_top) = self.inner.config().h_range.as_ref().map(|r| *r.end()) else { return };
        if window_top > top {
            return;
        }
        // keep vertices a canon cycle lands on: its `Q` copy lives at the boundary degree `top`, and
        // `eval_khi_elements` needs the vertex present (the homology window trims the class afterwards).
        let referenced: FxHashSet<TngComplexKey> = self.cone.elements().content().iter()
            .flat_map(|e| e.out_cob().keys().copied())
            .collect();

        let doomed = self.cone.complex().keys_of_deg(top + 1)
            .filter(|k| !referenced.contains(k))
            .copied().collect_vec();
        debug!("cone: drop out-of-window C[{}] ({} verts)", top + 1, doomed.len());
        self.cone.complex_mut().remove_vertices(&doomed);

        let total = self.cone.complex().keys_of_deg(top).count();
        let doomed = self.cone.complex().keys_of_deg(top)
            .filter(|k| self.cone.complex().vertex(k).in_edges().next().is_none() && !referenced.contains(k))
            .copied()
            .collect_vec();
        debug!("cone: window top C[{top}] {}/{total} no-in-edge verts pruned", doomed.len());
        self.cone.complex_mut().remove_vertices(&doomed);
    }

    /// The canon classes as KhI chains (cone bit → `B`/`Q`), expanded over each vertex's circles
    /// by the same pairing as `into_raw_complex`.
    pub fn eval_khi_elements(&self) -> Vec<KhIChain<R>> {
        let c = self.cone.complex();
        let (h, t) = c.ht().clone();

        // match `into_raw_complex_filtered`: skip out-of-window generators — dropped from the
        // complex, and the terms the full quotient mods out anyway.
        let q_range = self.cone.config().q_range.clone();
        let q_shift = c.deg_shift().1;
        let in_window = |g: &KhIGen|
            q_range.as_ref().is_none_or(|r| r.contains(&(q_shift + g.rel_q_deg())));

        self.cone.elements().content().iter().map(|e| {
            let init = LcCob::from(e.in_cob().clone());
            e.out_cob().iter().filter(|(k, _)| c.contains_key(k)).flat_map(|(k, retr)| {
                let circles = circles_of(c.vertex(k).tng());
                label_assignments(&circles).into_iter().filter_map(|b| {
                    let kg = into_khi_gen(&expanded_key(k, &b).as_gen());
                    if !in_window(&kg) {
                        return None;
                    }
                    let g = cap_circles(retr.clone(), End::Tgt, &circles, &b, &h, &t);
                    let x = init.stack(&g).eval(&h, &t);
                    Some((kg, x))
                }).collect_vec()
            }).collect()
        }).collect()
    }

    // Fully deloop the reduced cone so `eval_khi_elements` and `into_raw_complex` share the delooped
    // basis. Called from `run` after the last merge, matching the other builders.
    fn finalize(&mut self) {
        if self.skip_finalize() {
            info!("cone skip finalize (deloop deferred): {}", self.cone.stat());
            return;
        }
        info!("cone finalize: {}", self.cone.stat());
        self.deloop_all(false);
        self.deloop_all(true);
        info!("cone finalized: {}", self.cone.stat());
    }

    // Skip the finalize deloop when the build is complete and closed under `no_full_deloop`:
    // `into_raw_complex`/`eval_khi_elements` then expand the remaining circles at the matrix level.
    fn skip_finalize(&self) -> bool {
        self.inner.config().no_full_deloop
            && self.inner.n_nodes() == 0
            && self.cone.complex().is_closed()
    }

    // Deloop the whole cone one degree at a time, eliminating inline: a degree is eliminated once
    // its upper neighbor is delooped, so the delooped transient never spans more than the current
    // frontier. `eliminate_in` is capped by `max_elim_cost` — cheap pivots cascade here, heavy ones
    // survive to the matrix reduction. Called twice: non-based circles, then the based one.
    fn deloop_all(&mut self, based: bool) {
        debug!("cone deloop-all (based: {based})...");
        let auto_elim = self.cone.config().strategy.auto_elim();
        let range = self.cone.complex().h_range();
        let (start, end) = (*range.start(), *range.end());
        for d in start ..= end {
            self.cone.deloop_in_with(d, based);
            if auto_elim && d > start {
                self.cone.eliminate_in(d - 1);
            }
        }
        if auto_elim {
            self.cone.eliminate_in(end);
        }
    }

    /// Convert the cone to the KhI chain complex: each vertex expands into one generator per
    /// circle-label assignment, mapped to `KhIGen` via the cone bit (`into_khi_gen`) and ordered by
    /// `KhIGen` q-degree. `eval_khi_elements` reads the canon classes on this same basis.
    /// `h_range` restricts the conversion window (KhI degrees) — degrees outside are never expanded.
    pub fn into_raw_complex(self, h_range: Option<RangeInclusive<isize>>) -> ChainComplex1<KhIGen, R> {
        let q_range = self.cone.config().q_range.clone();
        self.cone.into_tng_complex().into_raw_complex_with(
            |k, a| into_khi_gen(&expanded_key(k, a).as_gen()),
            |g| g.rel_q_deg(),
            q_range,
            h_range,
        )
    }
}

// ---- cobordism-level cone construction (`1 + τ`, char-2) ----

// Config for the cone's own `TngComplexBuilder`, which drives deloop/eliminate on the coned complex:
// the simplify `strategy` and the elimination fill-cost cap carry over from the sym config.
fn cone_build_config(config: &SymBuildConfig) -> BuildConfig {
    BuildConfig { strategy: config.strategy, max_elim_cost: config.max_elim_cost, q_range: config.q_range.clone(), ..Default::default() }
}

// An empty cone shell: same `deg_shift`/base point, one extra h-degree for the cone bit.
fn cone_shell<R>(c: &TngComplex<R>) -> TngComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (h, t) = c.ht();
    TngComplex::new(h, t, c.deg_shift(), c.base_pt(), c.dim() + 1, Default::default())
}

// The cone copy of a symmetric-complex key: push the cone bit onto the state.
fn with_bit(k: &TngComplexKey, b: Bit) -> TngComplexKey {
    let mut key = *k;
    key.state.push(b);
    key
}

// The τ morphism on a closed tangle: one relabel cylinder `c → τc` per circle.
fn tau_cob<G>(tng: &Tng, inv_edge: G) -> Cob
where G: Fn(Edge) -> Edge {
    Cob::new(tng.comps().map(|c|
        CobComp::plain(Tng::from(c.clone()), Tng::from(c.convert_edges(&inv_edge)))
    ))
}

// Cone key → KhIGen: the cone bit is the last `state` bit (`0` → `B`/`Left`, `1` → `Q`/`Right`);
// stripping it gives the underlying `KhGen` (with the true q-grading).
fn into_khi_gen(x: &KhGen) -> KhIGen {
    let mut state = *x.state();
    let bit = state.iter().last().unwrap();
    state.remove(state.len() - 1);
    let under = KhGen::new(state, *x.tensor());
    match bit {
        Bit::Bit0 => KhIGen::from_left(under),
        Bit::Bit1 => KhIGen::from_right(under),
    }
}

// Lift a symmetric canon cycle to the cone layer `bit`: push the cone bit onto each `out_cob` key's
// state (`Bit0` = `B`, `Bit1` = `Q`); `in_cob`/`state` are unchanged.
fn lift_elem<R>(e: &TngComplexElem<R>, bit: Bit) -> TngComplexElem<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let mut e = e.clone();
    e.modify_out_cob(|out| out.into_iter().map(|(mut k, f)| {
        k.state.push(bit);
        (k, f)
    }).collect());
    e
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use yui_core::num::FF2;
    use yui_link::InvLink;
    use super::*;
    use super::super::{Strategy, CutOption};

    // Build the reduced cone, assert d² = 0, and return its nonzero homology ranks per degree.
    // (Full homology vs. the KhI reference is checked in `khi`.)
    fn cone_homology(l: &InvLink, reduced: bool, config: SymBuildConfig, full: bool) -> Vec<(isize, usize)> {
        let c = ConeBuilder::from_inv_link(l, &FF2::zero(), &FF2::zero(), reduced)
            .with_config(config).with_full_extend(full).run().into_raw_complex(None);
        c.check_d_all();
        let h = c.homology();
        h.support().map(|&i| (i, h[i].rank())).filter(|(_, r)| *r > 0).sorted().collect()
    }

    // ---- the two cone extensions agree ----

    // The symmetry-broken emission (`cone_extend_reduced`, Sano2026 Prop 4.6) is a deformation
    // retract of the doubled cone: its homology must agree with `cone_extend_full`'s.
    #[test]
    fn cone_reduced_matches_full() {
        for (name, l) in [
            ("3_1",    InvLink::test_data("3_1")),
            ("m3_1",   InvLink::test_data("3_1").mirror()),
            ("4_1",    InvLink::test_data("4_1")),
            ("5_1",    InvLink::test_data("5_1")),
            ("6_2a",   InvLink::test_data("6_2a")),
            ("m6_1a",  InvLink::test_data("6_1a").mirror()),
            ("6_3",    InvLink::test_data("6_3")),
            ("7_6a",   InvLink::test_data("7_6a")),
            ("8_21b",  InvLink::test_data("8_21b")),
            ("m9_46a", InvLink::test_data("9_46a").mirror()),
        ] {
            for reduced in [false, true] {
                let config = SymBuildConfig::default();
                let full = cone_homology(&l, reduced, config.clone(), true);
                let direct = cone_homology(&l, reduced, config, false);
                assert_eq!(full, direct, "{name}, reduced={reduced}");
            }
        }
    }

    // ---- the configuration does not change the homology ----

    // Chunking, deferring the final deloop, capping the elimination cost and the simplification
    // strategy are all performance knobs: every listed config must agree with the first.
    fn check_configs_agree(name: &str, l: &InvLink, configs: &[(&str, SymBuildConfig)]) {
        let (base_name, base) = &configs[0];
        for reduced in [false, true] {
            let expect = cone_homology(l, reduced, base.clone(), false);
            for (variant, config) in &configs[1..] {
                let got = cone_homology(l, reduced, config.clone(), false);
                assert_eq!(expect, got, "{name}: {variant} vs {base_name}, reduced={reduced}");
            }
        }
    }

    #[test]
    fn cone_config_independent() {
        for (name, l, chunks) in [
            ("3_1",  InvLink::test_data("3_1"), 2),
            ("4_1",  InvLink::test_data("4_1"), 2),
            ("5_2a", InvLink::test_data("5_2a"), 2),
            ("6_3",  InvLink::test_data("6_3"), 3),
            ("7_3a", InvLink::test_data("7_3a"), 3),
        ] {
            check_configs_agree(name, &l, &[
                ("default",            SymBuildConfig::default()),
                ("no-full-deloop",     SymBuildConfig { no_full_deloop: true, ..Default::default() }),
                ("chunked",            SymBuildConfig { cut: CutOption::Auto(chunks), ..Default::default() }),
                ("chunked, no-deloop", SymBuildConfig { cut: CutOption::Auto(chunks), no_full_deloop: true, ..Default::default() }),
                ("min-fill",           SymBuildConfig { strategy: Strategy::MinFill, ..Default::default() }),
            ]);
        }
    }

    // Capping the elimination cost pushes survivors into the matrix reduction, which is the
    // expensive direction — two diagrams are enough.
    #[test]
    fn cone_elim_cap_independent() {
        for (name, l, chunks) in [
            ("3_1", InvLink::test_data("3_1"), 2),
            ("6_3", InvLink::test_data("6_3"), 3),
        ] {
            check_configs_agree(name, &l, &[
                ("default",        SymBuildConfig::default()),
                ("cap 0",          SymBuildConfig { max_elim_cost: Some(0), ..Default::default() }),
                ("cap 4",          SymBuildConfig { max_elim_cost: Some(4), ..Default::default() }),
                ("chunked, cap 0", SymBuildConfig { cut: CutOption::Auto(chunks), max_elim_cost: Some(0), ..Default::default() }),
            ]);
        }
    }

    // `NoElim` / `None` skip the simplification entirely, so they blow up on anything but the
    // smallest diagram — `MinFill` is cheap and rides along with the other knobs above.
    #[test]
    fn cone_strategy_independent() {
        check_configs_agree("3_1", &InvLink::test_data("3_1"), &[
            ("greedy",      SymBuildConfig::default()),
            ("no-elim",     SymBuildConfig { strategy: Strategy::NoElim, ..Default::default() }),
            ("no-simplify", SymBuildConfig { strategy: Strategy::None, ..Default::default() }),
        ]);
    }

    // A windowed build's endpoint homology is wrong by construction, so only `d <= 0` compares.
    #[test]
    fn cone_windowed_chunked_9_46() {
        let l = InvLink::sym_pretzel(-3, 3, -3); // 9_46
        let window = Some(-64 ..= 1);
        let narrow = |h: Vec<(isize, usize)>| h.into_iter().filter(|&(d, _)| d <= 0).collect_vec();

        for reduced in [false, true] {
            let plain = cone_homology(&l, reduced, SymBuildConfig {
                h_range: window.clone(), ..Default::default()
            }, false);
            let chunked = cone_homology(&l, reduced, SymBuildConfig {
                cut: CutOption::Auto(2), strategy: Strategy::MinFill, h_range: window.clone(), ..Default::default()
            }, false);
            assert_eq!(narrow(plain), narrow(chunked), "reduced={reduced}");
        }
    }

    // ---- the canon classes ----

    // The cone's canon classes must give the same ssi as the matrix cone.
    #[test]
    fn cone_canon_ssi_matches_matrix() {
        use yui_core::poly::Poly;
        use crate::khi::KhIHomology;
        use crate::ss::{ssi_invariant_with, SsVersion};
        use crate::ss::div_vec;

        type P = Poly<'H', FF2>;
        let (c, t) = (P::variable(), P::zero());
        let knots = [
            ("3_1",  InvLink::test_data("3_1")),
            ("m3_1", InvLink::test_data("3_1").mirror()),
            ("4_1",  InvLink::test_data("4_1")),
            ("6_2a", InvLink::test_data("6_2a")),
            ("6_3",  InvLink::test_data("6_3")),
            ("7_6a", InvLink::test_data("7_6a")),
            // 9_46 and 8_21b have s̲ ≠ s̄ (ssi = (0, 2) and (2, 4)) — they exercise the canon-cycle
            // ordering, and are the diagrams here that reach `reduce_elements`' corrections.
            ("9_46",  InvLink::sym_pretzel(-3, 3, -3)),
            ("8_21b", InvLink::test_data("8_21b")),
        ];
        for (name, l) in knots {
            let matrix = ssi_invariant_with(&l, false, Default::default(), None, SsVersion::V1);

            let config = SymBuildConfig { h_range: Some(isize::MIN + 1 ..= 1), ..Default::default() };
            let kh = KhIHomology::new_with_config(&l, &c, &t, false, config);
            let zs = kh.canon_cycles();
            assert_eq!(zs.len(), 4);

            let ds = zs.iter().map(|z| {
                let h = kh.h_deg_of_chain(z);
                div_vec(&kh[h].vectorize_euc(z).subvec(0..2), &c).expect("invalid divisibility")
            }).collect_vec();

            let (w, r) = (l.writhe(), l.seifert_circles().len() as i32);
            let cone = (2 * ds[0] + w - r + 1, 2 * ds[2] + w - r + 1);
            assert_eq!(cone, matrix, "{name}");
        }
    }
}
