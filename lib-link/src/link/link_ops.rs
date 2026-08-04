//! Operations *on* a link — deriving a new diagram or a combinatorial object from the one at hand:
//! crossing changes, resolutions, and Seifert's algorithm. Contrast with [`crate::link::construct`],
//! which builds links out of patterns.

use petgraph::Graph;
use yui_core::{CloneAnd, Sign};
use yui_core::bitseq::Bit;

use super::{Edge, Link, LinkBuilder, Node, NodeOri, NodeType, Path, State};

impl Link {
    // Connected sum at the two base points (a PD-built link defaults to its minimal edge).
    pub fn conn_sum(&self, other: &Link) -> Link {
        let self_e = self.base_pt().expect("self needs a base point");
        let other_e = other.base_pt().expect("other needs a base point");
        self.conn_sum_at(other, self_e, other_e)
    }

    // Connected sum at edges `self_e`/`other_e`: import both links, then re-splice the two edges
    // tail-to-head (self out → other in, other out → self in), so the orientations are compatible.
    pub fn conn_sum_at(&self, other: &Link, self_e: Edge, other_e: Edge) -> Link {
        assert!(self.is_oriented(),  "conn_sum requires an oriented link (self)");
        assert!(other.is_oriented(), "conn_sum requires an oriented link (other)");

        let mut b = LinkBuilder::new();
        let v1 = b.add_link(self);
        let v2 = b.add_link(other);

        // the (tail, head) ends of each spliced edge, as builder ports.
        let port = |verts: &[_], (i, s): (usize, usize)| (verts[i], s);
        let (t1, h1) = self.edge_ends(self_e, true);
        let (t2, h2) = other.edge_ends(other_e, true);
        let (t1, h1) = (port(&v1, t1), port(&v1, h1));
        let (t2, h2) = (port(&v2, t2), port(&v2, h2));

        b.disconnect(h1);   // free self_e's two ports
        b.disconnect(h2);   // free other_e's two ports
        b.connect(t1, h2);  // self tail → other head
        b.connect(t2, h1);  // other tail → self head

        b.build().unwrap()
    }

    pub fn mirror(&self) -> Self {
        let l = Self::new(
            self.nodes().map(|x| x.mirror()),
            self.loops().iter().copied(),
        );
        // `new` defaults the base point to the minimal edge; mirroring preserves the edge set,
        // so re-imposing the original one is always valid.
        match self.base_pt() {
            Some(e) => l.with_base_pt(e),
            None => l,
        }
    }

    pub fn unoriented(&self) -> Self {
        if !self.is_oriented() {
            return self.clone();
        }
        let l = Self::new(
            self.nodes().map(|x| Node::new(x.node_type(), NodeOri::None, *x.edges())),
            self.loops().iter().copied(),
        );
        match self.base_pt() {
            Some(e) => l.with_base_pt(e),
            None => l,
        }
    }

    pub fn cc_at(&self, i: usize) -> Self {
        assert!(self.node(i).is_crossing());
        self.clone_and(|l|
            *l.node_mut(i) = l.node(i).mirror()
        )
    }

    pub fn resolve_at(&self, i: usize, r: Bit) -> Self {
        assert!(self.node(i).is_crossing());
        self.clone_and(|l|
            *l.node_mut(i) = l.node(i).resolve(r)
        )
    }

    pub fn resolve_by(&self, s: &State) -> Self {
        assert!(s.len() == self.n_crossings());

        let n = self.n_nodes();
        let itr = (0..n).filter(|&i| self.node(i).is_crossing());

        self.clone_and(|l| {
            for (i, r) in Iterator::zip(itr, s.iter()) {
                *l.node_mut(i) = self.node(i).resolve(r);
            }
        })
    }

    pub fn seifert_state(&self) -> State {
        assert!(self.is_oriented());

        let seq = self.crossings().map(|x|
            match x.sign() {
                Some(Sign::Pos) => 0,
                Some(Sign::Neg) => 1,
                None => panic!("Impossible.")
            }
        );
        State::from_iter(seq)
    }

    pub fn seifert_circles(&self) -> Vec<Path> {
        self.resolve_by(&self.seifert_state()).comps()
    }

    pub fn seifert_graph(&self) -> Graph<Path, usize> {
        assert!(self.is_oriented());

        type G = Graph<Path, usize>;

        let s0 = self.seifert_state();
        let l0 = self.resolve_by(&s0);
        let mut graph = Graph::new();

        // Vertices = Seifert circles (and free loops contribute their own circles).
        for c in l0.comps() {
            graph.add_node(c);
        }

        let find_node = |graph: &G, e| {
            graph.node_indices().find(|&i|
                graph[i].contains(e)
            )
        };

        // Edges = one per original crossing (now resolved into a V/H smoothing).
        // Free loops have no nodes, so they remain isolated vertices.
        for (i, x) in l0.nodes().enumerate() {
            let (e1, e2) = if x.node_type() == NodeType::V {
                (x.edge(0), x.edge(1))
            } else {
                (x.edge(0), x.edge(2))
            };
            let n1 = find_node(&graph, e1).unwrap();
            let n2 = find_node(&graph, e2).unwrap();
            graph.add_edge(n1, n2, i);
        }

        graph
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Braid, Node, NodeOri};
    use crate::NodeType::{XL, XR};
    use crate::misc::jones_polynomial;

    #[test]
    fn link_mirror() {
        let l = Link::test_data("unknot_l_twist");
        assert_eq!(l.node(0).node_type(), XL);

        let l = l.mirror();
        assert_eq!(l.node(0).node_type(), XR);
    }

    #[test]
    fn mirror_preserves_loops() {
        let l = Link::unlink(1).mirror();
        assert_eq!(l.n_loops(), 1);
        assert_eq!(l.loops(), &[1]);
    }

    #[test]
    fn mirror_preserves_base_pt() {
        let l = Link::test_data("3_1").with_base_pt(3).mirror();
        assert_eq!(l.base_pt(), Some(3));
    }

    #[test]
    fn crossing_change() {
        let l = Link::test_data("3_1");
        let l2 = l.cc_at(1);

        assert_eq!(l.node(1),  &Node::new(XL, NodeOri::Up, [3,6,4,1]));
        assert_eq!(l2.node(1), &Node::new(XR, NodeOri::Up, [3,6,4,1]));
    }

    #[test]
    fn link_resolve() {
        let s = State::from([0, 0, 0]);
        let l = Link::test_data("3_1").resolve_by(&s);

        let comps = l.comps();
        assert_eq!(comps.len(), 3);
        assert!(comps.iter().all(|c| c.is_circle()));

        let s = State::from([1, 1, 1]);
        let l = Link::test_data("3_1").resolve_by(&s);

        let comps = l.comps();
        assert_eq!(comps.len(), 2);
        assert!(comps.iter().all(|c| c.is_circle()));
    }

    #[test]
    fn seifert_graph_trefoil() {
        let l = Link::test_data("3_1");
        let g = l.seifert_graph();
        assert_eq!(g.node_count(), 2);
        assert_eq!(g.edge_count(), 3);
    }

    #[test]
    fn seifert_graph_unlink() {
        // Free loops contribute isolated vertices and no edges.
        let l = Link::unlink(3);
        let g = l.seifert_graph();
        assert_eq!(g.node_count(), 3);
        assert_eq!(g.edge_count(), 0);
    }

    #[test]
    fn conn_sum_is_jones_multiplicative() {
        // unreduced Jones: Ṽ(K1 # K2) · Ṽ(unknot) = Ṽ(K1) · Ṽ(K2)
        let k1 = Link::test_data("3_1");
        let k2 = Link::test_data("4_1");
        let cs = k1.conn_sum(&k2); // at the default base points (edge 1 each)
        assert_eq!(cs.n_comps(), 1);
        assert_eq!(cs.n_crossings(), k1.n_crossings() + k2.n_crossings());
        assert!(cs.is_oriented());

        let (vcs, vu) = (jones_polynomial(&cs), jones_polynomial(&Link::unknot()));
        let (v1, v2) = (jones_polynomial(&k1), jones_polynomial(&k2));
        assert_eq!(&vcs * &vu, &v1 * &v2);

        // the connected sum does not depend on the chosen edges.
        let cs2 = k1.conn_sum_at(&k2, 4, 6);
        assert_eq!(jones_polynomial(&cs2), jones_polynomial(&cs));
    }

    #[test]
    fn conn_sum_of_braid_closures() {
        // braid closures orient downward (`ori = Down`), exercising the non-PD arms of
        // `NodeOri::in_ports` in edge_dir.
        let k1 = Braid::from([1, 1, 1]).closure();      // 3_1
        let k2 = Braid::from([1, -2, 1, -2]).closure(); // 4_1
        let cs = k1.conn_sum(&k2);
        assert_eq!(cs.n_comps(), 1);
        assert!(cs.is_oriented());

        let (vcs, vu) = (jones_polynomial(&cs), jones_polynomial(&Link::unknot()));
        let (v1, v2) = (jones_polynomial(&k1), jones_polynomial(&k2));
        assert_eq!(&vcs * &vu, &v1 * &v2);
    }
}
