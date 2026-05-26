use std::collections::HashSet;

use itertools::Itertools;
use log::{debug, info};
use num_traits::Zero;
use yui_core::bitseq::Bit;
use yui_core::{Ring, RingOps};
use yui_link::{Node, Edge, Link};

use crate::kh::{KhChain, KhComplex};
use crate::tng::{TngComplexElem, LcCobTrait, TngComplex, TngComplexKey};

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

    pub(crate) fn remove_nodes<'a, I>(&mut self, nodes: I) 
    where I: IntoIterator<Item = &'a Node> { 
        let drop = nodes.into_iter().collect::<HashSet<_>>();
        self.nodes.retain(|x| !drop.contains(x));
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

    pub(crate) fn process_nodes(&mut self) { 
        while let Some(x) = self.choose_next_node() { 
            self.append_node(&x)
        }
    }

    pub(crate) fn choose_next_node(&mut self) -> Option<Node> { 
        let Some((i, _)) = self.nodes.iter().enumerate().max_by_key(|(_, x)|
            self.count_connections(x)
        ) else { 
            return None
        };

        let x = self.nodes.remove(i);
        Some(x)
    }

    fn count_connections(&self, x: &Node) -> usize { 
        let arcs = if x.is_resolved() { 
            let a = x.arcs();
            vec![a.0, a.1]
        } else { 
            let a0 = x.resolve(Bit::Bit0).arcs();
            let a1 = x.resolve(Bit::Bit1).arcs();
            vec![a0.0, a0.1, a1.0, a1.1]
        }.into_iter().filter(|a|
            self.complex.base_pt().map(|e| !a.contains(e)).unwrap_or(true)
        ).collect_vec();

        let count = self.complex.iter_verts().map(|(_, v)| {
            v.tng().comps().map(|c| 
                arcs.iter().map(|a| 
                    match c.path() {
                        p if p.is_connectable_bothends(a) => 2,
                        p if p.is_connectable(a)          => 1,
                        _                                 => 0
                    }
                ).sum::<usize>()
            ).sum::<usize>()
        }).sum::<usize>();

        count
    }

    pub(crate) fn append_node(&mut self, x: &Node) { 
        info!("({}) append: {x}", self.stat());

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
        info!("({}) merge <- ({})", self.stat(), other.stat());

        let (left, right) = self.complex.prepare_merge(other); 

        self.complex.merge_vertices(&left, &right);

        for i in self.complex.h_range() { 
            self.complex.merge_edges(&left, &right, i);
            if self.auto_deloop {
                self.deloop_in(i, false);
            }
        }
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
        let mut keys = self.complex.keys_of_deg(i).filter(|k| 
            self.is_deloopable(k, allow_based)
        ).sorted_by_key(|&k| self.complex.vertex(k).c_weight()).cloned().collect_vec();

        if keys.is_empty() { return }

        info!("({}) C[{i}]: {}, deloop targets: {}.", self.stat(), self.complex.rank(i), keys.len());

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

        info!("({}) -> C[{i}]: {} (diff: {}).", self.stat(), after, after - before);
    }

    pub(crate) fn deloop(&mut self, k: &TngComplexKey, r: usize) -> Vec<TngComplexKey> {
        let c = self.complex.vertex(k).tng().comp(r);

        debug!("({}) deloop {c} in {}", self.stat(), self.complex.vertex(k));

        for e in self.elements.iter_mut() { 
            e.deloop(k, c);
        }

        let added = self.complex.deloop(k, r);

        if self.auto_elim { 
            let mut res = vec![];
            for k in added.iter() { 
                if let Some(&j) = self.choose_inv_edge_into(&k) { 
                    self.eliminate(&j, &k);
                } else if let Some(&l) = self.choose_inv_edge_from(&k) { 
                    self.eliminate(&k, &l);
                } else { 
                    res.push(*k);
                }
            }
            res
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
        debug!("({}) eliminate {}: {} -> {}", self.stat(), self.complex.edge(i, j), self.complex.vertex(i), self.complex.vertex(j));
        
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

        let a = self.complex.edge(i, j);
        let ainv = a.inv().unwrap();
        let (h, t) = self.complex.ht();

        for k in self.complex.vertex(i).out_edges() { 
            if k == j { continue }

            let c = self.complex.edge(i, k);
            let cab = c * &ainv * &b;
            let s = if let Some(d) = e.remove_cob(k) {
                d - cab
            } else {
                -cab
            }.part_eval(h, t);

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

        assert!(self.complex.is_completely_delooped());
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