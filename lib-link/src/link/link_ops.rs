//! Operations *on* a link — deriving a new diagram or a combinatorial object from the one at hand:
//! crossing changes, resolutions, and Seifert's algorithm. Contrast with [`crate::link::construct`],
//! which builds links out of patterns.

use petgraph::Graph;
use yui_core::num::Sign;
use yui_core::ext::CloneAnd;
use yui_core::bitseq::Bit;

use super::{Edge, Link, LinkBuilder, NodeType, Path, Slot, State};

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
        assert!(
            !self.loops().contains(&self_e) && !other.loops().contains(&other_e),
            "connected sum on free loops is not supported yet"
        );

        let mut b = LinkBuilder::new();
        let v1 = b.add_link(self);
        let v2 = b.add_link(other);

        // the (tail, head) ends of each spliced edge, as builder ports.
        let port = |verts: &[_], (i, s): (usize, Slot)| (verts[i], s.index());
        let (t1, h1) = self.edge_ends(self_e, true);
        let (t2, h2) = other.edge_ends(other_e, true);
        let (t1, h1) = (port(&v1, t1), port(&v1, h1));
        let (t2, h2) = (port(&v2, t2), port(&v2, h2));

        b.disconnect(h1);   // free self_e's two ports
        b.disconnect(h2);   // free other_e's two ports
        b.connect(t1, h2);  // self tail → other head
        b.connect(t2, h1);  // other tail → self head

        // A base point elsewhere survives the splice; one sitting on the consumed `self_e` moves
        // onto the band, taking the `self tail → other head` side. Read the ids off before `build`
        // renumbers. Basing on the band matters for `K # K`, whose axis is the band itself.
        let base = self.base_pt().map(|e|
            if e != self_e {
                b.edge_at(port(&v1, self.edge_ends(e, false).0)).unwrap()
            } else {
                b.edge_at(t1).unwrap()
            }
        );

        // Orientation is inherited: the splice runs tail → head both ways, so every node keeps the
        // incoming slots it had in its own summand. Plain `build` only knows that *some* port is
        // incoming, and would be free to reverse the whole diagram.
        let n1 = self.n_nodes();
        let sum = b.build_with(|i, s| {
            let (l, i) = if i < n1 { (self, i) } else { (other, i - n1) };
            l.node(i).is_incoming(Slot::from(s))
        }).unwrap();

        match base {
            Some(e) => sum.with_base_pt(e),
            None => sum,
        }
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

    // Reverse the orientation. Both strands of every crossing turn around together, so the signs
    // — hence the writhe — are unchanged; only the direction of travel is.
    pub fn reversed(&self) -> Self {
        let l = Self::new(
            self.nodes().map(|x| x.reversed()),
            self.loops().iter().copied(),
        );
        // `new` defaults the base point to the minimal edge; reversing keeps the edge set, so the
        // original one is still valid.
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
        self.clone_and(|l| {
            *l.node_mut(i) = l.node(i).resolve(r);
            l.normalize_ori();
        })
    }

    pub fn resolve_by(&self, s: &State) -> Self {
        assert!(s.len() == self.n_crossings());

        let n = self.n_nodes();
        let itr = (0..n).filter(|&i| self.node(i).is_crossing());

        self.clone_and(|l| {
            for (i, r) in Iterator::zip(itr, s.iter()) {
                *l.node_mut(i) = self.node(i).resolve(r);
            }
            l.normalize_ori();
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
                (x.edge(Slot::SW), x.edge(Slot::SE))
            } else {
                (x.edge(Slot::SW), x.edge(Slot::NE))
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
    use crate::{Braid, Node};
    use crate::NodeType::{XL, XR};
    use crate::misc::jones_polynomial;

    #[test]
    fn link_reversed() {
        let l = Link::test_data("3_1");
        let r = l.reversed();

        assert!(r.is_oriented());
        assert_eq!(r.writhe(), l.writhe(), "reversing a knot keeps every crossing sign");
        assert_eq!(r.base_pt(), l.base_pt());

        // the underlying diagram is untouched — only the direction of travel changes.
        for (x, y) in Iterator::zip(l.nodes(), r.nodes()) {
            assert_eq!(y.node_type(), x.node_type());
            assert_eq!(y.edges(), x.edges());
        }
        assert_eq!(r.reversed(), l, "reversing twice is the identity");

        // the two incoming slots move to the far end of their own strands.
        for (x, y) in Iterator::zip(l.nodes(), r.nodes()) {
            let (p, q) = x.incoming().unwrap();
            assert_eq!(y.incoming(), Some((x.paired_slot(p).min(x.paired_slot(q)),
                                           x.paired_slot(p).max(x.paired_slot(q)))));
        }
    }

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

        assert_eq!(l.node(1),  &Node::new(XL, Some((Slot::SW, Slot::NW)), [3,1,4,6]));
        assert_eq!(l2.node(1), &Node::new(XR, Some((Slot::SW, Slot::NW)), [3,1,4,6]));
    }

    #[test]
    fn resolve_drops_the_orientation_wholesale() {
        // a smoothing that does not respect the orientation costs the whole diagram its own; only
        // the Seifert state keeps it.
        let l = Link::test_data("3_1");
        assert_eq!(l.seifert_state(), State::from([0, 0, 0]));

        assert!(l.resolve_by(&State::from([0, 0, 0])).is_oriented());
        for st in [[0, 0, 1], [0, 1, 1], [1, 1, 1]] {
            let r = l.resolve_by(&State::from(st));
            assert!(!r.is_oriented(), "{st:?} left a partial orientation");
            r.verify_ori();
        }
    }

    #[test]
    fn link_resolve() {
        let s = State::from([0, 0, 0]);
        let l = Link::test_data("3_1").resolve_by(&s);

        let comps = l.comps();
        assert_eq!(comps.len(), 2);
        assert!(comps.iter().all(|c| c.is_circle()));

        let s = State::from([1, 1, 1]);
        let l = Link::test_data("3_1").resolve_by(&s);

        let comps = l.comps();
        assert_eq!(comps.len(), 3);
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
        // braid closures orient downward, exercising the non-PD arms of `edge_ends`.
        let k1 = Braid::from([1, 1, 1]).closure();      // 3_1
        let k2 = Braid::from([1, -2, 1, -2]).closure(); // 4_1
        let cs = k1.conn_sum(&k2);
        assert_eq!(cs.n_comps(), 1);
        assert!(cs.is_oriented());

        let (vcs, vu) = (jones_polynomial(&cs), jones_polynomial(&Link::unknot()));
        let (v1, v2) = (jones_polynomial(&k1), jones_polynomial(&k2));
        assert_eq!(&vcs * &vu, &v1 * &v2);
    }

    // The four curl sums below come out on the builder's own labelling, so each PD code matches
    // with no reindexing: the band edges are 1 and 3, each curl keeps its own loop (2 and 4).
    // The splice consumes the base point, which lands on the band at edge 3.

    #[test]
    fn conn_sum_of_two_curls() {
        // The smallest non-trivial splice: the 1-crossing left curl with itself is the two-curl
        // unknot diagram. Pins the resulting *diagram*, not just an invariant of it.
        let k = Link::test_data("unknot_l_twist");
        let cs = k.conn_sum(&k);
        assert_eq!(cs.pd_code(), [[1,3,2,2],[3,1,4,4]]);
        assert_eq!(cs.base_pt(), Some(3));
    }

    #[test]
    fn conn_sum_of_a_curl_and_its_reverse() {
        // Same splice with the second curl traversed the other way: the summand enters the band
        // by its other end, so the two curls sit head-to-head.
        let k = Link::test_data("unknot_l_twist");
        let cs = k.conn_sum(&k.reversed());
        assert_eq!(cs.pd_code(), [[1,3,2,2],[4,4,1,3]]);
        assert_eq!(cs.base_pt(), Some(3));
    }

    #[test]
    fn conn_sum_of_a_curl_and_its_mirror() {
        // Mirroring flips the crossing sign but not the direction of travel, so the sum is the
        // R2-cancelling pair — a different diagram from the reversed case above.
        let k = Link::test_data("unknot_l_twist");
        let cs = k.conn_sum(&k.mirror());
        assert_eq!(cs.pd_code(), [[1,3,2,2],[4,3,1,4]]);
        assert_eq!(cs.base_pt(), Some(3));
    }

    #[test]
    fn conn_sum_of_a_curl_and_its_concordance_inverse() {
        // Mirror *and* reverse: the fourth of the four ways to glue the second curl on, and the
        // fourth distinct diagram.
        let k = Link::test_data("unknot_l_twist");
        let cs = k.conn_sum(&k.mirror().reversed());
        assert_eq!(cs.pd_code(), [[1,3,2,2],[3,4,4,1]]);
        assert_eq!(cs.base_pt(), Some(3));
    }

    #[test]
    fn conn_sum_at_of_two_curls_away_from_the_base_pt() {
        // Splicing at edge 2 instead — the curl's own loop — gives the two-curl diagram bundled as
        // `unknot_l_twist2`. Edge 1 carries the base point and the splice leaves it alone, so it
        // survives in place rather than moving onto the band.
        let k = Link::test_data("unknot_l_twist");
        let cs = k.conn_sum_at(&k, 2, 2);
        assert_eq!(cs.base_pt(), Some(1));

        // Unlike the four sums above, the raw code is [[1,1,4,2],[3,3,2,4]]: `LinkBuilder::build`
        // numbers edges in `connect` order, and here the two band edges are connected last, so
        // they take the highest ids instead of following the strand. Renumbering along the
        // orientation from the base point recovers the expected labelling.
        assert_eq!(cs.reindexed(1, 1).pd_code(), [[1,1,2,4],[3,3,4,2]]);
        assert_eq!(cs.reindexed(1, 1).pd_code(), Link::test_data("unknot_l_twist2").pd_code());
    }

    #[test]
    fn conn_sum_at_of_two_curls_on_different_edges() {
        // Asymmetric splice: self's loop (edge 2) to other's base edge (edge 1).
        let k = Link::test_data("unknot_l_twist");
        let cs = k.conn_sum_at(&k, 2, 1);
        assert_eq!(cs.base_pt(), Some(1));
        // raw is [[1,1,3,2],[3,2,4,4]], and the target below is not traversal-numbered either, so
        // both sides are renumbered from the base point before comparing.
        let expected = Link::from_pd_code([[1,1,2,3],[2,3,4,4]]);
        assert_eq!(cs.reindexed(1, 1).pd_code(), expected.reindexed(1, 1).pd_code());
    }

    #[test]
    fn conn_sum_keeps_each_summand_orientation() {
        // The splice joins tail to head both ways, so every node keeps the incoming slots it had
        // in its own summand. `add_link` appends self's nodes then other's, in order.
        // Reversing a summand is the case that catches it: a PD-loaded diagram happens to have
        // node 0's SW incoming, which is also what an unconstrained `build` picks.
        for (a, b) in [("3_1", "4_1"), ("4_1", "3_1"), ("unknot_l_twist", "3_1"), ("6_2", "5_1")] {
            for (rev1, rev2) in [(false, false), (true, false), (false, true), (true, true)] {
                let k1 = Link::test_data(a).clone_and(|l| if rev1 { *l = l.reversed() });
                let k2 = Link::test_data(b).clone_and(|l| if rev2 { *l = l.reversed() });
                let cs = k1.conn_sum(&k2);
                let n1 = k1.n_nodes();
                let case = format!("{a}{} # {b}{}", if rev1 { "*" } else { "" }, if rev2 { "*" } else { "" });

                for (i, x) in k1.nodes().enumerate() {
                    assert_eq!(cs.node(i).incoming(), x.incoming(), "{case}: node {i} of {a}");
                }
                for (i, x) in k2.nodes().enumerate() {
                    assert_eq!(cs.node(n1 + i).incoming(), x.incoming(), "{case}: node {i} of {b}");
                }
            }
        }
    }

    #[test]
    #[should_panic(expected = "free loops is not supported")]
    fn conn_sum_rejects_a_free_loop_base_pt() {
        // K # unknot = K mathematically, but the splice needs a base point sitting on a crossing.
        let _ = Link::unknot().conn_sum(&Link::test_data("3_1"));
    }
}
