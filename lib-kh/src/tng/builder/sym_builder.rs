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
use ahash::{AHashMap, AHashSet};
use cartesian::cartesian;
use itertools::Itertools;
use log::{debug, info};
use yui_core::bitseq::{Bit, BitSeq};
use yui_core::{RangeExt, Ring, RingOps};
use yui_link::{Node, Edge, InvLink};

use crate::kh::{KhGen, KhTensor};
use crate::tng::{LcCobTrait, TngComp, TngComplex, TngComplexKey};
use crate::tng::builder::{TngComplexBuilder, BuildConfig};

/// Toggles for the automatic simplification done while building (kept separate
/// from [`BuildConfig`] so the equivariant builder can gain its own flags).
#[derive(Clone, Copy, Debug)]
pub struct SymBuildConfig {
    pub auto_deloop: bool,
    pub auto_elim: bool,
}

impl Default for SymBuildConfig {
    fn default() -> Self {
        Self { auto_deloop: true, auto_elim: true }
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
        self.config = config;
        self
    }

    pub fn run(mut self) -> Self {
        self.process_nodes();
        self.finalize();
        self
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

        for i in self.inner.complex().h_range().mv(0, 1) { 
            self.inner.complex_mut().merge_vertices(&left, &right, i);
            self.inner.complex_mut().merge_edges(&left, &right, i - 1);

            if self.config.auto_elim { 
                self.eliminate_in(i - 1);
            }
            if self.config.auto_deloop { 
                self.deloop_in(i - 1, false);
            }
        }

        // TODO merge elements
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