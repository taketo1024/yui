//! Equivariant variant of [`TngComplexBuilder`] for strongly invertible links:
//! processes axis-symmetric crossings singly and off-axis crossings in pairs
//! `(x, τx)`, then maintains a `key_map` matching each [`TngComplexKey`] with
//! its τ-image so that [`SymTngBuilder::tau_map`] can be used to assemble the
//! involutive Khovanov complex `CKhI = Cone(1 + τ)`.
//!
//! Reference:
//! - T. Sano, "Involutive Khovanov homology and equivariant knots",
//!   Algebr. Geom. Topol. 25 (2025), 5059–5111.
//!   <https://doi.org/10.2140/agt.2025.25.5059>, <https://arxiv.org/abs/2404.08568>

use std::collections::HashSet;
use std::ops::RangeInclusive;
use ahash::{AHashMap, AHashSet};
use cartesian::cartesian;
use itertools::Itertools;
use log::{debug, info};
use yui_core::algo::KeyedUnionFind;
use yui_core::bitseq::{Bit, BitSeq};
use yui_core::{RangeExt, Ring, RingOps};
use yui_link::{Node, Edge, InvLink};

use crate::kh::{KhGen, KhTensor};
use crate::tng::{LcCobTrait, TngComp, TngComplex, TngComplexElem, TngComplexKey};
use crate::tng::builder::{TngComplexBuilder, BuildConfig};

/// Toggles for the automatic simplification done while building (kept separate
/// from [`BuildConfig`] so the equivariant builder can gain its own flags).
#[derive(Clone, Debug)]
pub struct SymBuildConfig {
    pub auto_deloop: bool,
    pub auto_elim: bool,
    // build half the off-axis crossings and mirror via τ (see `preprocess`).
    pub preprocess: bool,
    // literal truncation: homology at the endpoints is wrong (build `(a-1)..=(b+1)` for correct `[a, b]`).
    pub h_range: Option<RangeInclusive<isize>>,
}

impl Default for SymBuildConfig {
    fn default() -> Self {
        Self { auto_deloop: true, auto_elim: true, preprocess: true, h_range: None }
    }
}

pub struct SymTngBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    inner: TngComplexBuilder<R>,
    x_map: AHashMap<Node, Node>,
    e_map: AHashMap<Edge, Edge>,
    key_map: AHashMap<TngComplexKey, TngComplexKey>,
    config: SymBuildConfig,
}

impl<R> SymTngBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_inv_link(l: &InvLink, h: &R, t: &R, reduced: bool) -> SymTngBuilder<R> {
        assert!(l.nodes().all(|x| x.is_crossing()));
        assert!(!reduced || l.base_pt().is_some());

        // the inner builder is driven by `self` — disable its own auto-simplify.
        let inner = TngComplexBuilder::from_link(l.inner(), h, t, reduced)
            .with_config(BuildConfig { auto_deloop: false, auto_elim: false, h_range: None });

        let x_map = l.nodes().map(|x|
            (x.clone(), l.inv_node(x).clone())
        ).collect();
        let e_map = l.edges().into_iter().map(|e| (e, l.inv_edge(e))).collect();
        let key_map = AHashMap::from_iter([(TngComplexKey::init(), TngComplexKey::init())]);

        SymTngBuilder { inner, x_map, e_map, key_map, config: SymBuildConfig::default() }
    }

    pub fn with_config(mut self, config: SymBuildConfig) -> Self {
        // propagate the window to the inner builder so the preprocess merges cap
        // to it; this also drops canon cycles when the window excludes h-degree 0.
        let inner_config = BuildConfig { auto_deloop: false, auto_elim: false, h_range: config.h_range.clone() };
        self.inner = self.inner.with_config(inner_config);
        self.config = config;
        self
    }

    pub fn run(mut self) -> Self {
        if self.config.preprocess {
            self.preprocess();
        }
        self.process_nodes();
        self.finalize();
        self
    }

    // Build one representative half of the off-axis crossings, mirror it via τ,
    // merge both, then leave the on-axis crossings for `process_nodes`.
    fn preprocess(&mut self) {
        assert_eq!(self.inner.complex().dim(), 0, "must start from init state.");

        let elements = self.inner.take_elements();
        let off_axis = self.off_axis_crossings(true).into_iter().cloned().collect_vec();

        info!("({}) preprocess off-axis: {}", self.inner.stat(), off_axis.len());

        let (c, tc, key_map, elements) = self.build_from_half(off_axis.iter(), elements);

        // merge left
        self.inner.drop_nodes(|x| off_axis.contains(x));
        self.inner.merge(c);

        // merge right
        let off_axis = off_axis.iter().map(|x| self.inv_node(x).clone()).collect_vec();
        self.inner.drop_nodes(|x| off_axis.contains(x));
        self.inner.merge(tc);

        self.key_map = key_map;
        self.inner.set_elements(elements);

        info!("({}) preprocess done.", self.inner.stat());
    }

    // Off-axis crossings; when `take_half`, one representative per τ-pair,
    // grouped by adjacency so a connected side is taken whole.
    fn off_axis_crossings(&self, take_half: bool) -> Vec<&Node> {
        let off_axis = self.inner.nodes().filter(|&x|
            self.inv_node(x) != x
        ).collect_vec();

        if !take_half {
            return off_axis
        }

        let is_adj = |x: &Node, y: &Node| {
            x.edges().iter().filter(|&&e|
                self.inv_edge(e) != e // no axis-crossing edge
            ).any(|e|
                y.edges().contains(e)
            )
        };

        let mut u = KeyedUnionFind::from_iter(off_axis.iter().cloned());

        for (i, x) in off_axis.iter().enumerate() {
            for j in 0 .. i {
                let y = &off_axis[j];
                if is_adj(x, y) {
                    u.union(x, y);
                }
            }
        }

        u.into_disjoint().into_iter().fold(vec![], |mut res, next| {
            if let Some(x) = next.first() {
                let tx = self.inv_node(x);
                if !res.contains(&tx) {
                    res.extend(next);
                }
            }
            res
        })
    }

    fn build_from_half<'a, I>(&self, crossings: I, elements: Vec<TngComplexElem<R>>) -> (TngComplex<R>, TngComplex<R>, AHashMap<TngComplexKey, TngComplexKey>, Vec<TngComplexElem<R>>)
    where I: IntoIterator<Item = &'a Node> {
        let (h, t) = self.inner.complex().ht();
        let base_pt = self.inner.complex().base_pt();
        let mut b = TngComplexBuilder::init(h, t, (0, 0), base_pt);

        // a half-vertex of weight > b - deg_shift.0 can never land in the window.
        if let Some(h_range) = &self.config.h_range {
            let s = self.inner.complex().deg_shift().0;
            b = b.with_config(BuildConfig { h_range: Some(0 ..= (*h_range.end() - s)), ..Default::default() });
        }

        b.set_nodes(crossings.into_iter().cloned());
        b.set_elements(elements);
        b.process_nodes();

        // take results
        let keys = b.complex().keys().cloned().collect_vec();
        let elements = b.take_elements();
        let c = b.into_tng_complex();
        let tc = c.convert_edges(|e| self.inv_edge(e));

        // build keys: k1 (from c) + k2 (from tc) ↦ τ(k) = k2 + k1
        let keys = cartesian!(keys.iter(), keys.iter());
        let key_map = keys.map(|(k1, k2)| {
            let k  = k1 + k2;
            let tk = k2 + k1;
            (k, tk)
        }).collect();

        // duplicate each element by connecting its τ-image
        let elements = elements.into_iter().map(|mut e| {
            e.modify(|k, c| {
                let kk = k + k;
                let cc = c.map(|c, r| {
                    let r = &r * &r;
                    let tc = c.convert_edges(|e| self.inv_edge(e));
                    (c.connect(&tc), r)
                });
                (kk, cc)
            });
            e
        }).collect();

        (c, tc, key_map, elements)
    }

    fn process_nodes(&mut self) {
        info!("process {} nodes", self.inner.nodes().count());

        while let Some(x) = self.choose_next_node_sym().cloned() {
            let tx = self.inv_node(&x).clone();
            if x == tx {
                self.append_on_axis(&x);
            } else {
                self.append_off_axis(&x, &tx);
            }

            info!("  appended: {x}, current size: {}", self.inner.stat());
        }
    }

    /// Pair-aware chooser: an off-axis node is scored as the pair `(x, τx)`
    /// it'll be appended as; on-axis is doubled to compare at the same scale.
    fn choose_next_node_sym(&self) -> Option<&Node> {
        let boundary_ends: AHashSet<Edge> = self.inner.complex().boundary_ends().collect();
        self.inner.nodes()
            .max_by_key(|x| {
                let tx = self.inv_node(x);
                let (l_x, w_x) = self.inner.score_node(x, &boundary_ends);
                if tx == *x {
                    (2 * l_x, 2 * w_x)
                } else {
                    let (l_tx, w_tx) = self.inner.score_node(tx, &boundary_ends);
                    (l_x + l_tx, w_x + w_tx)
                }
            })
    }

    fn append_on_axis(&mut self, x: &Node) { 
        info!("({}/{}) append on-axis: {x}", 
            self.inner.complex().dim() + 1, 
            self.inner.complex().dim() + self.inner.nodes().count() + 1, 
        );

        self.inner.prepare_append(x);

        let (h, t) = self.inner.complex().ht();
        let c = TngComplex::from_node(h, t, x, self.inner.complex().base_pt());
        let key_map = if x.is_crossing() {
            [Bit::Bit0, Bit::Bit1].map(|b| { 
                let k = TngComplexKey { state: BitSeq::from(b), label: KhTensor::empty() };
                (k, k)
            }).into_iter().collect()
        } else { 
            let k = TngComplexKey::init();
            [(k, k)].into_iter().collect()
        };

        self.merge(c, key_map);
    }

    fn append_off_axis(&mut self, x: &Node, tx: &Node) { 
        assert_eq!(self.inv_node(x), tx);

        info!("({}/{}) append off-axis: {x}, {tx}", 
            self.inner.complex().dim() + 1, 
            self.inner.complex().dim() + self.inner.nodes().count() + 1, 
        );

        self.inner.prepare_append(x);
        self.inner.prepare_append(tx);

        let c = {
            let (h, t) = self.inner.complex().ht();
            let mut c = TngComplex::from_node(h, t, x, self.inner.complex().base_pt());
            c.append_node(tx);
            c
        };
        let key_map = if x.is_crossing() { 
            [
                ([0, 0], [0, 0]),
                ([1, 0], [0, 1]),
                ([0, 1], [1, 0]),
                ([1, 1], [1, 1])
            ].map(|(b0, b1)| { 
                let k = TngComplexKey { state: BitSeq::from_iter(b0), label: KhTensor::empty() };
                let l = TngComplexKey { state: BitSeq::from_iter(b1), label: KhTensor::empty() };
                (k, l)
            }).into_iter().collect()
        } else {
            let k = TngComplexKey::init();
            [(k, k)].into_iter().collect()
        };

        self.merge(c, key_map);
    }

    fn merge(&mut self, c: TngComplex<R>, key_map: AHashMap<TngComplexKey, TngComplexKey>) { 
        self.key_map = cartesian!(
            self.key_map.iter(),
            key_map.iter()
        ).map(|((k1, l1), (k2, l2))|
            (k1 + k2, l1 + l2)
        ).collect();

        let (left, right) = self.inner.complex_mut().prepare_merge(c);

        // cap at the top of `h_range`: weight only grows, so higher is unreachable.
        let range = self.inner.complex().h_range().mv(0, 1);
        let top = match &self.config.h_range {
            Some(h_range) => *range.end().min(h_range.end()),
            None => *range.end(),
        };

        for i in *range.start() ..= top {
            self.inner.complex_mut().merge_vertices(&left, &right, i);
            self.inner.complex_mut().merge_edges(&left, &right, i - 1);

            if self.config.auto_elim {
                self.eliminate_in(i - 1);
            }
            if self.config.auto_deloop {
                self.deloop_in(i - 1, false);
            }
        }

        self.prune_h_range();

        // TODO merge elements
    }

    /// Prune doomed vertices *and* their `key_map` entries. τ preserves weight,
    /// so `deg(k) == deg(τk)` — pruning is symmetric and `key_map` stays a valid
    /// involution. The `key_map` is pruned independently of the complex because
    /// the cartesian merge over-generates it past the vertex cap.
    fn prune_h_range(&mut self) {
        let Some(h_range) = self.config.h_range.clone() else { return };
        let r = self.inner.nodes().count() as isize;
        let i0 = self.inner.complex().deg_shift().0;

        // a vertex of degree `d` reaches `[d, d + r]`, so it stays relevant iff
        // `d ∈ [a - r, b]` — current degrees that can still land in `h_range`.
        let live = (*h_range.start() - r) ..= *h_range.end();
        let doomed = |k: &TngComplexKey|
            !live.contains(&(k.weight() as isize + i0));

        let doomed_verts = self.inner.complex().keys_of(&doomed).copied().collect_vec();
        if !doomed_verts.is_empty() {
            debug!("prune {} verts outside h_range.", doomed_verts.len());
        }
        self.inner.prune_keys(&doomed_verts);

        let doomed_keys = self.key_map.keys().filter(|k| doomed(k)).copied().collect_vec();
        for k in &doomed_keys {
            self.key_map.remove(k);
        }
    }

    fn deloop_all(&mut self, allow_based: bool) {
        for i in self.inner.complex().h_range() {
            self.deloop_in(i, allow_based);
        }
    }

    fn deloop_in(&mut self, i: isize, allow_based: bool) {
        let mut keys = self.inner.pick_keys_in(i, |k|
            self.inner.is_deloopable(k, allow_based)
        );
        if keys.is_empty() { return }

        debug!("deloop in C[{i}], targets: {}.", keys.len());

        let before = self.inner.complex().rank(i) as isize;

        while !keys.is_empty() { 
            let k = keys.remove(0);
            if !self.inner.complex().contains_key(&k) { continue; } // already delooped or eliminated

            let mut list = vec![k];

            while !list.is_empty() { 
                let k = list.remove(0);
                if !self.inner.complex().contains_key(&k) { continue; }
                let Some(r) = self.inner.find_loop(&k, allow_based) else { continue };

                let added = self.deloop_equiv(&k, r);

                list.extend(added.into_iter().filter(|k|
                    self.inner.is_deloopable(k, allow_based)
                ));
            }
        }

        let after = self.inner.complex().rank(i) as isize;

        debug!("  delooped C[{i}]: {} (diff: {}).", after, after - before);
    }

    fn deloop_equiv(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> { 
        let added = if self.is_sym_key(k) { 
            let c = self.inner.complex().vertex(k).tng().comp(r);
            if self.is_sym_comp(c) {
                // symmetric loop on symmetric key
                self.deloop_on_axis_sym(k, r)
            } else {
                // asymmetric loop on symmetric key
                self.deloop_on_axis_asym(k, r)
            }
        } else { 
            // (symmetric or asymmetric) loop on asymmetric key
            self.deloop_off_axis(k, r)
        };

        if self.config.auto_elim { 
            added.into_iter().filter(|k| 
                self.inner.complex().contains_key(&k) && 
                !self.try_eliminate_equiv_at(&k)
            ).collect()
        } else { 
            added
        }
    }

    fn deloop_on_axis_sym(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> {
        let c = self.inner.complex().vertex(k).tng().comp(r);

        assert!(self.is_sym_key(k));
        assert!(self.is_sym_comp(c));

        let updated = self.inner.deloop(k, r);

        self.remove_key_pair(k);

        for &k_new in updated.iter() { 
            self.add_key_pair(k_new, k_new);
        }

        updated
    }

    #[allow(non_snake_case)]
    fn deloop_on_axis_asym(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> {
        let c = self.inner.complex().vertex(k).tng().comp(r);

        assert!(self.is_sym_key(k));
        assert!(!self.is_sym_comp(c));
        assert!(!c.is_marked());

        //          ⚪︎1 | ⚪︎1
        //  ⚪︎1 | ⚪︎X  <-->  ⚪︎X | ⚪︎1
        //          ⚪︎X | ⚪︎X

        let tc = c.convert_edges(|e| self.inv_edge(e));

        let ks = self.inner.deloop(k, r);

        let (k_X, k_1) = (ks[0], ks[1]);
        let (k_XX, k_X1) = { 
            let tr = self.inner.complex().vertex(&k_X).tng().index_of(&tc).unwrap();
            let tks = self.inner.deloop(&k_X, tr);
            (tks[0], tks[1])
        };
        let (k_1X, k_11) = { 
            let tr = self.inner.complex().vertex(&k_1).tng().index_of(&tc).unwrap();
            let tks = self.inner.deloop(&k_1, tr);
            (tks[0], tks[1])
        };

        self.remove_key_pair(k);

        self.add_key_pair(k_XX, k_XX);
        self.add_key_pair(k_X1, k_1X);
        self.add_key_pair(k_11, k_11);

        vec![k_XX, k_X1, k_1X, k_11]
    }

    #[allow(non_snake_case)]
    fn deloop_off_axis(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> {
        let c = self.inner.complex().vertex(k).tng().comp(r);

        assert!(!self.is_sym_key(k));

        //  ⚪︎1 | ..  <-->  .. | ⚪︎1
        //  ⚪︎X | ..  <-->  .. | ⚪︎X

        let tk = *self.inv_key(k);
        let tc = c.convert_edges(|e| self.inv_edge(e));
        let tr = self.inner.complex().vertex(&tk).tng().index_of(&tc).unwrap();

        let mut ks = self.inner.deloop(k, r);
        let mut tks = self.inner.deloop(&tk, tr);

        self.remove_key_pair(k);

        for (&k_new, &tk_new) in Iterator::zip(ks.iter(), tks.iter()) { 
            self.add_key_pair(k_new, tk_new);
        }

        ks.append(&mut tks);
        ks
    }

    fn eliminate_in(&mut self, i: isize) {
        let mut keys = self.inner.pick_keys_in(i, |k|
            self.inner.complex().vertex(k).out_edges()
                .any(|l| self.is_equiv_inv_edge(k, l))
        );
        if keys.is_empty() { return }

        debug!("eliminate in C[{i}], targets: {}", keys.len());
        
        let before = self.inner.complex().rank(i) as isize;

        while !keys.is_empty() { 
            let k = keys.remove(0);
            if !self.inner.complex().contains_key(&k) { continue; } // removed by other side
            self.try_eliminate_equiv_at(&k);
        }

        let after = self.inner.complex().rank(i) as isize;

        debug!("  eliminated C[{i}]: {} (diff: {}).", after, after - before);
    }

    fn choose_equiv_inv_edge_into(&self, k: &TngComplexKey) -> Option<&TngComplexKey> { 
        self.inner.complex().vertex(k).in_edges().filter_map(|j|
            self.is_equiv_inv_edge(j, k).then_some(j)
        )
        .min_by_key(|j| (self.inner.edge_weight(j, k), **j))
    }

    fn choose_equiv_inv_edge_from(&self, k: &TngComplexKey) -> Option<&TngComplexKey> { 
        self.inner.complex().vertex(k).out_edges().filter_map(|l|
            self.is_equiv_inv_edge(k, l).then_some(l)
        )
        .min_by_key(|l| (self.inner.edge_weight(k, l), **l))
    }

    fn try_eliminate_equiv_at(&mut self, k: &TngComplexKey) -> bool {
        if let Some(&j) = self.choose_equiv_inv_edge_into(&k) { 
            self.eliminate_equiv(&j, &k);
            true
        } else if let Some(&l) = self.choose_equiv_inv_edge_from(&k) { 
            self.eliminate_equiv(&k, &l);
            true
        } else { 
            false
        }
    }

    fn eliminate_equiv(&mut self, i: &TngComplexKey, j: &TngComplexKey) {
        assert_eq!(self.is_sym_key(i), self.is_sym_key(j));
        assert!(self.inner.complex().has_edge(i, j));

        if self.is_sym_key(i) { 
            self.inner.eliminate(i, j);
        } else { 
            let ti = *self.inv_key(i);
            let tj = *self.inv_key(j);

            assert!(self.inner.complex().has_edge(&ti, &tj));

            self.inner.eliminate(i, j);
            self.inner.eliminate(&ti, &tj);
        }

        self.remove_key_pair(i);
        self.remove_key_pair(j);
    }

    fn is_equiv_inv_edge(&self, i: &TngComplexKey, j: &TngComplexKey) -> bool { 
        let f = self.inner.complex().edge(i, j);
        f.is_invertible() && self.is_equiv_edge(i, j)
    }

    fn is_equiv_edge(&self, i: &TngComplexKey, j: &TngComplexKey) -> bool { 
        if self.is_sym_key(i) && self.is_sym_key(j) { 
            true
        } else if !self.is_sym_key(i) && !self.is_sym_key(j) { 
            //  i - - -> j 
            //    \   /   
            //      /     : not allowed
            //    /   \   
            // ti - - -> tj
            let ti = self.inv_key(i);
            let tj = self.inv_key(j);

            !self.inner.complex().vertex(j).in_edges().contains(ti) && 
            !self.inner.complex().vertex(tj).in_edges().contains(i)
        } else { 
            false
        }
    }

    fn finalize(&mut self) {
        if self.inner.complex().is_completely_delooped() { 
            return
        }

        info!("finalize: {}", self.inner.stat());

        self.deloop_all(false);
        self.deloop_all(true);

        info!("  finalized: {}", self.inner.stat());
    }

    pub fn tau_map(&self) -> impl Fn(&KhGen) -> KhGen + Send + Sync + 'static {
        let key_map = self.key_map.clone();

        move |x: &KhGen| -> KhGen {
            let k = TngComplexKey::from(x);
            let tk = key_map[&k];
            tk.as_gen()
        }
    }

    pub fn into_inner(self) -> TngComplexBuilder<R> { 
        self.inner
    }

    pub fn into_tng_complex(self) -> TngComplex<R> { 
        self.inner.into_tng_complex()
    }

    fn inv_node(&self, x: &Node) -> &Node { 
        &self.x_map[x]
    }

    fn inv_edge(&self, e: Edge) -> Edge { 
        self.e_map[&e]
    }

    fn inv_key(&self, k: &TngComplexKey) -> &TngComplexKey { 
        &self.key_map[k]
    }

    fn add_key_pair(&mut self, k: TngComplexKey, tk: TngComplexKey) { 
        if let Some(l) = self.key_map.get(&tk) {
            assert_eq!(k, *l);
            return;
        }

        self.key_map.insert(k, tk);
        if k != tk { 
            self.key_map.insert(tk, k);
        }
    }

    fn remove_key_pair(&mut self, k: &TngComplexKey) { 
        let tk = self.key_map.remove(k).unwrap();
        if k != &tk { 
            self.key_map.remove(&tk);
        }
    }

    fn is_sym_key(&self, k: &TngComplexKey) -> bool { 
        self.inv_key(k) == k
    }

    fn is_sym_comp(&self, c: &TngComp) -> bool { 
        &c.convert_edges(|e| self.inv_edge(e)) == c
    }

    #[allow(unused)]
    fn print_keys(&self) {
        let mut done = HashSet::new();
        for k in self.key_map.keys().sorted() { 
            if done.contains(&k) { continue }

            let tk = self.inv_key(k);
            if k == tk {
                println!("{}", self.inner.complex().vertex(k));
            } else { 
                println!("{} ↔ {}", self.inner.complex().vertex(k), self.inner.complex().vertex(tk));
            }

            done.insert(k);
            done.insert(tk);
        }
        println!();
    }

}

#[cfg(test)]
#[allow(unused)]
mod tests { 
    use crate::khi::{KhIGen, KhIHomology};

    use super::*;
    use num_traits::Zero;

    use yui_core::IteratorExt;
    use yui_core::lc::Lc;
    use yui_core::num::FF2;
    use yui_core::poly::Poly;
    use yui_core::RangeExt;
    use yui_homology::{ChainComplex1, ChainMap, ToSeqString, ToTableString};

    fn make_cone(b: SymTngBuilder<FF2>) -> ChainComplex1<KhIGen, FF2> { 
        let t = b.tau_map();
        let c = b.into_inner().into_tng_complex().into_raw_complex();
        let h_range = c.support().cloned().range().unwrap().mv(0, 1);
        let one_plus_tau = ChainMap::new(&c, &c, 0, move |_, z| {
            z.clone() + z.apply(|x| Lc::from(t(x)))
        });
        one_plus_tau.cone(h_range, false)
    }

    #[test]
    fn preprocess_matches() {
        let l = InvLink::test_data("6_3");
        let (h, t) = (FF2::zero(), FF2::zero());

        let build = |preprocess: bool, window: Option<(isize, isize)>| {
            let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
            b.config.preprocess = preprocess;
            b.config.h_range = window.map(|(a, b)| a..=b);
            b.run().into_tng_complex().into_raw_complex()
        };

        // full build: preprocess on/off agree on every degree.
        let (c_on, c_off) = (build(true, None), build(false, None));
        let range = c_off.support().cloned().range().unwrap();
        let (h_on, h_off) = (c_on.homology(), c_off.homology());
        for i in range.clone() {
            assert_eq!(h_on[i].rank(), h_off[i].rank(), "full rank at {i}");
        }

        // windowed build: agree on the interior of the window.
        let (a, b) = (*range.start() + 1, *range.end() - 1);
        let (c_on, c_off) = (build(true, Some((a, b))), build(false, Some((a, b))));
        let (h_on, h_off) = (c_on.homology(), c_off.homology());
        for i in (a + 1)..=(b - 1) {
            assert_eq!(h_on[i].rank(), h_off[i].rank(), "windowed rank at {i}");
        }
    }

    #[test]
    fn test_kh_3_1() {
        let l = InvLink::test_data("3_1");
        let (h, t) = (FF2::zero(), FF2::zero());

        let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
        b.config.auto_elim = false;
        let c = b.run().into_tng_complex().into_raw_complex();
        c.check_d_all();

        let h = c.homology();

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 2);
    }

    #[test]
    fn test_kh_3_1_h_range_full() {
        // A range covering the whole complex must reproduce the full homology.
        let l = InvLink::test_data("3_1");
        let (h, t) = (FF2::zero(), FF2::zero());

        let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
        b.config.auto_elim = false;
        b.config.h_range = Some(0..=3);
        let c = b.run().into_tng_complex().into_raw_complex();
        c.check_d_all();

        let h = c.homology();

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 2);
    }

    #[test]
    fn test_kh_6_3_h_range() {
        let l = InvLink::test_data("6_3");
        let (h, t) = (FF2::zero(), FF2::zero());

        // full build for reference.
        let full = SymTngBuilder::from_inv_link(&l, &h, &t, false).run()
            .into_tng_complex().into_raw_complex();
        let range = full.support().cloned().range().unwrap();
        let (lo, hi) = (*range.start(), *range.end());
        let full_h = full.homology();

        // restrict to a strict sub-window.
        let (a, b) = (lo + 1, hi - 1);
        let mut bld = SymTngBuilder::from_inv_link(&l, &h, &t, false);
        bld.config.h_range = Some(a..=b);
        let bld = bld.run();

        // key_map stays in sync with the truncated complex.
        let verts: HashSet<_> = bld.inner.complex().keys().copied().collect();
        let kmap: HashSet<_> = bld.key_map.keys().copied().collect();
        assert_eq!(verts, kmap);

        let c = bld.into_tng_complex().into_raw_complex();
        c.check_d_all();
        let trunc_h = c.homology();

        // chain groups vanish outside [a, b].
        assert_eq!(c[a - 1].rank(), 0);
        assert_eq!(c[b + 1].rank(), 0);

        // strict interior degrees of the window are correct.
        for i in (a + 1)..=(b - 1) {
            assert_eq!(trunc_h[i].rank(), full_h[i].rank(), "rank mismatch at h-deg {i}");
        }
    }

    #[test]
    fn test_khi_3_1() {
        let l = InvLink::test_data("3_1");
        let (h, t) = (FF2::zero(), FF2::zero());

        let b = SymTngBuilder::from_inv_link(&l, &h, &t, false).run();
        let c = make_cone(b);
        c.check_d_all();

        let h = c.homology();

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 4);
        assert_eq!(h[4].rank(), 2);
    }

    #[test]
    fn no_auto_deloop() {
        let l = InvLink::test_data("3_1");
        let (h, t) = (FF2::zero(), FF2::zero());

        let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
        b.config.auto_deloop = false;
        b.process_nodes();

        assert!(!b.inner.complex().is_completely_delooped());

        b.finalize();

        assert!(b.inner.complex().is_completely_delooped());

        let c = b.into_tng_complex().into_raw_complex();
        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);
        assert_eq!(c[2].rank(), 2);
        assert_eq!(c[3].rank(), 2);
        
        let h = c.homology();
        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 2);
    }

    #[test]
    fn no_auto_elim() { 
        let l = InvLink::test_data("3_1");
        let (h, t) = (FF2::zero(), FF2::zero());

        let mut b = SymTngBuilder::from_inv_link(&l, &h, &t, false);
        b.config.auto_elim = false;
        b.process_nodes();

        assert!(b.inner.complex().is_completely_delooped());

        let c = b.into_tng_complex().into_raw_complex();
        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 12);
        assert_eq!(c[3].rank(), 8);
        
        let h = c.homology();
        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 2);
    }
}