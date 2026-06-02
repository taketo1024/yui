//! Incremental builder for [`TngComplex`]: scan crossings one at a time,
//! tensor-merge with the new crossing's small complex, then deloop newborn
//! circles and gauss-eliminate invertible edges to keep the complex small.
//!
//! References:
//! - BN05 — D. Bar-Natan, "Khovanov's homology for tangles and cobordisms",
//!   Geom. Topol. 9 (2005), 1443–1499.
//!   <https://doi.org/10.2140/gt.2005.9.1443>, <https://arxiv.org/abs/math/0410495>
//! - BN07 — D. Bar-Natan, "Fast Khovanov homology computations",
//!   J. Knot Theory Ramif. 16 (2007), 243–255.
//!   <https://doi.org/10.1142/S0218216507005294>, <https://arxiv.org/abs/math/0606318>

use ahash::AHashSet;
use itertools::Itertools;
use log::{debug, info, trace};
use num_traits::Zero;
use yui_core::bitseq::Bit;
use yui_core::{RangeExt, Ring, RingOps};
use yui_link::{Node, Edge, Link};

use crate::kh::{KhChain, KhComplex};
use crate::tng::{TngComp, TngComplexElem, LcCobTrait, TngComplex, TngComplexKey};

pub struct TngComplexBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    complex: TngComplex<R>,
    nodes: Vec<Node>,
    loops: Vec<Edge>, 
    elements: Vec<TngComplexElem<R>>,
    pub auto_deloop: bool,
    pub auto_elim: bool
}

impl<R> TngComplexBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_link(l: &Link, h: &R, t: &R, reduced: bool) -> Self { 
        let base_pt = if reduced { l.base_pt() } else { None };
        let deg_shift = KhComplex::deg_shift_for(l, reduced);

        let mut b = Self::init(h, t, deg_shift, base_pt);
        b.set_nodes(l.nodes().cloned());
        b.set_loops(l.loops().iter().cloned());

        if t.is_zero() && l.is_knot() {
            let canon = TngComplexElem::canon_cycles(l, base_pt);
            b.set_elements(canon);
        }

        b
    }

    pub fn init(h: &R, t: &R, deg_shift: (isize, isize), base_pt: Option<Edge>) -> Self { 
        let complex = TngComplex::init(h, t, deg_shift, base_pt);
        Self { 
            complex, 
            nodes: vec![], 
            loops: vec![],
            elements: vec![], 
            auto_deloop: true, 
            auto_elim: true 
        }
    }

    pub fn complex(&self) -> &TngComplex<R> { 
        &self.complex
    }

    pub(crate) fn complex_mut(&mut self) -> &mut TngComplex<R> { 
        &mut self.complex
    }

    pub fn nodes(&self) -> impl Iterator<Item = &Node> { 
        self.nodes.iter()
    }

    pub fn set_nodes<I>(&mut self, nodes: I)
    where I: IntoIterator<Item = Node> {
        self.nodes = nodes.into_iter().collect_vec();
    }

    pub fn loops(&self) -> impl Iterator<Item = &Edge> { 
        self.loops.iter()
    }

    pub fn set_loops<I>(&mut self, loops: I)
    where I: IntoIterator<Item = Edge> {
        self.loops = loops.into_iter().collect_vec();
    }

    pub fn set_elements<I>(&mut self, elements: I)
    where I: IntoIterator<Item = TngComplexElem<R>> { 
        self.elements = elements.into_iter().collect_vec();
    }

    pub(crate) fn take_elements(&mut self) -> Vec<TngComplexElem<R>> {
        std::mem::take(&mut self.elements)
    }

    pub fn run(mut self) -> Self { 
        self.process_nodes();
        self.process_loops();
        self.finalize();
        self
    } 

    // See [BN07, §7] (scan-and-cancel algorithm).
    pub(crate) fn process_nodes(&mut self) {
        info!("process {} nodes", self.nodes.len());

        while let Some(x) = self.choose_next_node().cloned() {
            self.append_node(&x)
        }
    }

    /// Pick the next node to append by maximizing [`Self::score_node`].
    /// Removal happens in [`Self::prepare_append`] inside `append_*`.
    pub(crate) fn choose_next_node(&self) -> Option<&Node> {
        let boundary_ends: AHashSet<Edge> = self.complex.boundary_ends().collect();
        self.nodes.iter().max_by_key(|x| self.score_node(x, &boundary_ends))
    }

    /// Score `x` for the chooser. Higher is better.
    /// `(loop_bonus, width_score)` — loop closures first (they unlock
    /// delooping + elimination), width as tiebreaker.
    pub(crate) fn score_node(&self, x: &Node, boundary_ends: &AHashSet<Edge>) -> (usize, isize) {
        let arcs = self.node_arcs(x);

        let loops: usize = self.complex.iter_verts().map(|(_, v)| {
            v.tng().comps().map(|c|
                arcs.iter().filter(|a| c.is_connectable_bothends(a)).count()
            ).sum::<usize>()
        }).sum();

        let width_score: isize = arcs.iter().map(|a| {
            let Some((e0, e1)) = a.end_pts() else { return 0 };
            let m = boundary_ends.contains(&e0) as isize
                  + boundary_ends.contains(&e1) as isize;
            2 * m - 2 // matched − unmatched, range −2..2 per arc
        }).sum();

        (loops, width_score)
    }

    /// Arcs that will be added when appending `x` to the partial diagram,
    /// with the base-point edge filtered out.
    fn node_arcs(&self, x: &Node) -> Vec<TngComp> {
        let arcs = if x.is_resolved() {
            let (a0, a1) = x.arcs();
            vec![a0, a1]
        } else {
            let (a00, a01) = x.resolve(Bit::Bit0).arcs();
            let (a10, a11) = x.resolve(Bit::Bit1).arcs();
            vec![a00, a01, a10, a11]
        };
        let base_pt = self.complex.base_pt();
        arcs.into_iter()
            .filter(|a| base_pt.map(|e| !a.contains(e)).unwrap_or(true))
            .map(TngComp::from)
            .collect()
    }

    pub(crate) fn append_node(&mut self, x: &Node) { 
        info!("({}/{}) append: {x}", 
            self.complex.dim() + 1, 
            self.complex.dim() + self.nodes.len() + 1, 
        );

        self.prepare_append(x);

        let (h, t) = self.complex.ht();
        let cx = TngComplex::from_node(h, t, x);
        self.merge(cx);
    }

    pub(crate) fn prepare_append(&mut self, x: &Node) { 
        if let Some(i) = self.nodes.iter().find_position(|&e| e == x) { 
            self.nodes.remove(i.0);
        }

        for e in self.elements.iter_mut() { 
            e.append_node(x);
        }
    }

    pub(crate) fn merge(&mut self, other: TngComplex<R>) { 
        debug!("merge {} + {}", self.stat(), other.stat());

        let (left, right) = self.complex.prepare_merge(other); 

        for i in self.complex.h_range().mv(0, 1) { 
            self.complex.merge_vertices(&left, &right, i);
            self.complex.merge_edges(&left, &right, i - 1);
            
            if self.auto_elim { 
                self.eliminate_in(i - 1);
            }
            if self.auto_deloop {
                self.deloop_in(i - 1, false);
            }
        }

        debug!("  merged: {}", self.stat());
    }

    pub(crate) fn process_loops(&mut self) { 
        while !self.loops.is_empty() { 
            let c = self.loops.remove(0);

            for e in self.elements.iter_mut() { 
               e.insert_loop(c);
            }

            let (h, t) = self.complex.ht();
            let c = TngComplex::from_loop(h, t, c);
            self.merge(c);

            if self.auto_deloop { 
                self.deloop_all(false);
            }
        }
    }

    pub(crate) fn deloop_all(&mut self, allow_based: bool) {
        for i in self.complex.h_range() { 
            self.deloop_in(i, allow_based);
        }
    }

    fn deloop_in(&mut self, i: isize, allow_based: bool) {
        let mut keys = self.pick_keys_in(i, |k|
            self.is_deloopable(k, allow_based)
        );
        if keys.is_empty() { return }

        debug!("deloop in C[{i}], targets: {}.", keys.len());

        let before = self.complex.rank(i) as isize;

        while !keys.is_empty() { 
            let k = keys.remove(0);
            let mut list = vec![k];

            while !list.is_empty() { 
                let k = list.remove(0);
                if !self.complex.contains_key(&k) { continue; }
                let Some(r) = self.find_loop(&k, allow_based) else { continue };

                let added = self.deloop(&k, r);

                list.extend(added.into_iter().filter(|k|
                    self.is_deloopable(k, allow_based)
                ));
            }
        }

        let after = self.complex.rank(i) as isize;

       debug!("  delooped C[{i}]: {} (diff: {}).", after, after - before);
    }

    pub(crate) fn deloop(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> {
        let c = self.complex.vertex(k).tng().comp(r);

        trace!("{} deloop {c} in {}", self.stat(), self.complex.vertex(k));

        for e in self.elements.iter_mut() { 
            e.deloop(k, c);
        }

        let added = self.complex.deloop(k, r);

        if self.auto_elim { 
            // only retain keys are not eliminated
            added.into_iter().filter(|k|
                !self.try_eliminate_at(k)
            ).collect_vec()
        } else { 
            added
        }
    }

    pub(crate) fn is_deloopable(&self, k: &TngComplexKey, allow_based: bool) -> bool { 
        self.find_loop(k, allow_based).is_some()
    }

    pub(crate) fn find_loop(&self, k: &TngComplexKey, allow_based: bool) -> Option<usize> { 
        self.complex.vertex(k).tng().find_comp(|c|
            c.is_circle() && (allow_based || !self.complex.contains_base_pt(c))
        )
    }

    /// Keys at degree `i` matching `pred`, sorted ascending by the vertex's
    /// `c_weight` (so the cheapest vertices come first).
    pub(crate) fn pick_keys_in<F>(&self, i: isize, pred: F) -> Vec<TngComplexKey>
    where F: Fn(&TngComplexKey) -> bool {
        self.complex.keys_of_deg(i)
            .filter(|k| pred(k))
            .sorted_by_key(|k| self.complex.vertex(k).c_weight())
            .copied()
            .collect_vec()
    }

    pub(crate) fn eliminate_in(&mut self, i: isize) {
        let mut keys = self.pick_keys_in(i, |k|
            self.complex.vertex(k).out_edges().any(|l|
                self.complex.edge(k, l).is_invertible()
            )
        );
        if keys.is_empty() { return }

        debug!("eliminate in C[{i}], targets: {}", keys.len());

        let before = self.complex.rank(i);

        while !keys.is_empty() { 
            let k = keys.remove(0);
            self.try_eliminate_at(&k);
        }

        let after = self.complex.rank(i);

        debug!("  eliminated C[{i}]: {} (diff: {}).", after, after - before);
    }

    pub(crate) fn try_eliminate_at(&mut self, k: &TngComplexKey) -> bool {
        if let Some(&j) = self.choose_inv_edge_into(&k) { 
            self.eliminate(&j, &k);
            true
        } else if let Some(&l) = self.choose_inv_edge_from(&k) { 
            self.eliminate(&k, &l);
            true
        } else { 
            false
        }
    }

    fn choose_inv_edge_into(&self, k: &TngComplexKey) -> Option<&TngComplexKey> { 
        self.complex.vertex(k).in_edges().filter_map(|j|
            self.complex.edge(j, k).is_invertible().then_some(j)
        )
        .min_by_key(|j| self.edge_weight(j, k))
    }

    fn choose_inv_edge_from(&self, k: &TngComplexKey) -> Option<&TngComplexKey> { 
        self.complex.vertex(k).out_edges().filter_map(|l|
            self.complex.edge(k, l).is_invertible().then_some(l)
        )
        .min_by_key(|l| self.edge_weight(k, l))
    }

    pub(crate) fn edge_weight(&self, k: &TngComplexKey, l: &TngComplexKey) -> usize { 
        let nk = self.complex.vertex(k).out_edges().count(); // nnz in column k
        let nl = self.complex.vertex(l).in_edges().count();     // nnz in row l
        (nk - 1) * (nl - 1)
    }

    pub(crate) fn eliminate(&mut self, i: &TngComplexKey, j: &TngComplexKey) {
        trace!("{} eliminate {}: {} -> {}", self.stat(), self.complex.edge(i, j), self.complex.vertex(i), self.complex.vertex(j));
        
        self.eliminate_elements(i, j);
        self.complex.eliminate(i, j);
    }

    pub(crate) fn eliminate_elements(&mut self, i: &TngComplexKey, j: &TngComplexKey) {
        let mut elements = self.take_elements();
        for e in elements.iter_mut() { 
            self.eliminate_element(e, i, j);
        }
        self.elements = elements;
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
    
    fn eliminate_element(&self, e: &mut TngComplexElem<R>, i: &TngComplexKey, j: &TngComplexKey) {
        assert!(self.complex.has_edge(i, j));

        // mors into i can be simply dropped.
        e.remove_cob(i);

        // mors into j must be redirected by -ca^{-1}
        let Some(b) = e.remove_cob(j) else { return };

        let (h, t) = self.complex.ht();
        let a = self.complex.edge(i, j);
        let ainv = a.inv().unwrap();
        let ainv_b = ainv * b;

        for k in self.complex.vertex(i).out_edges() { 
            if k == j { continue }

            let c = self.complex.edge(i, k);
            let c_ainv_b = (c * &ainv_b).reduce(h, t);
            let s = if let Some(d) = e.remove_cob(k) {
                d - c_ainv_b
            } else {
                -c_ainv_b
            };

            if !s.is_zero() { 
                e.insert_cob(*k, s);
            }
        }
    }

    fn finalize(&mut self) { 
        if self.complex.is_completely_delooped() { 
            return;
        }

        info!("finalize");

        self.deloop_all(false);
        self.deloop_all(true);
    }

    pub fn into_tng_complex(self) -> TngComplex<R> { 
        self.complex
    }

    pub fn eval_elements(&self) -> Vec<KhChain<R>> {
        let (h, t) = self.complex.ht();
        self.elements.iter().map(|z|
            z.eval(h, t)
        ).collect()
    }

    pub(crate) fn stat(&self) -> String {
        self.complex.stat()
    }
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    
    use super::*;

    #[test]
    fn test_unknot() {
        let l = Link::unknot();
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).run();
        let c = b.into_tng_complex().into_raw_complex();

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_unknot_rm1() {
        let l = Link::test_data("unknot_l_twist");
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).run();
        let c = b.into_tng_complex().into_raw_complex();

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_unknot_rm1_neg() {
        let l = Link::test_data("unknot_r_twist");
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).run();
        let c = b.into_tng_complex().into_raw_complex();

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
    }

    #[test]
    fn test_unknot_rm2() {
        let l = Link::test_data("unknot_lr_twist");
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).run();
        let c = b.into_tng_complex().into_raw_complex();

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
        assert_eq!(c[ 1].rank(), 0);
    }

    #[test]
    fn test_unlink_2() {
        let l = Link::test_data("unlink2");
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).run();
        let c = b.into_tng_complex().into_raw_complex();

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 4);
        assert_eq!(c[ 1].rank(), 0);
    }

    #[test]
    fn test_tangle() { 
        let mut c = TngComplexBuilder::init(&0, &0, (0, 0), None);
        c.set_nodes([
            Node::from_pd_code([4,2,5,1]),
            Node::from_pd_code([3,6,4,1])
        ]);

        c.process_nodes();
        
        assert!(!c.complex.is_completely_delooped());
    }

    #[test]
    fn test_hopf_link() {
        let l = Link::test_data("L2a1");
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).run();
        let c = b.into_tng_complex().into_raw_complex();

        c.check_d_all();

        assert_eq!(c[-2].rank(), 2);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[ 0].rank(), 2);
    }

    #[test]
    fn test_8_19() {
        let l = Link::test_data("8_19");
        let b = TngComplexBuilder::from_link(&l, &0, &0, false).run();
        let c = b.into_tng_complex().into_raw_complex();

        c.check_d_all();

        let h = c.homology();

        for i in [1,6,7,8] {
            assert_eq!(h[i].rank(), 0);
            assert!(h[i].is_free());
        }

        for i in [0,4,5] {
            assert_eq!(h[i].rank(), 2);
            assert!(h[i].is_free());
        }

        assert_eq!(h[2].rank(), 1);
        assert!(h[2].is_free());

        assert_eq!(h[3].rank(), 1);
        assert_eq!(h[3].tors(), &vec![2]);
    }

    #[test]
    fn canon_cycle_trefoil() { 
        let l = Link::test_data("3_1");
        let b = TngComplexBuilder::from_link(&l, &1, &0, false).run();
        let zs = b.eval_elements();
        let c = b.into_tng_complex().into_raw_complex();

        assert_eq!(zs.len(), 2);
        assert_ne!(zs[0], zs[1]);
        
        for z in zs {
            assert!(c.d(0, &z).is_zero());
        }
    }
}