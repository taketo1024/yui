// Constructions on `Link` kept out of the core type.

use itertools::Itertools;

use crate::{Link, Edge, LinkBuilder, Port};

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

    // Twist knot = numerator closure of the rational tangle [n, 2]: a horizontal |n|-twist region
    // summed with a vertical 2-twist clasp. twist_knot(0) = unknot, (1,2,3,4) = 3_1, 4_1, 5_2, 6_1;
    // n < 0 mirrors.
    pub fn twist_knot(n: i32) -> Link {
        use crate::NodeType::{XL, XR};
        let mut b = LinkBuilder::new();

        // the vertical 2-twist clasp; corners (v0, v1, v2, v3) = (SW, SE, NE, NW).
        let (v0, v1, v2, v3) = b.add_v_twist(XR, 2);

        if n == 0 {
            b.connect(v3, v2);
            b.connect(v0, v1);   // N([2] clasp) = unknot
            return b.build().unwrap();
        }

        // the horizontal |n|-twist region.
        let tt = if n >= 0 { XR } else { XL };
        let (h0, h1, h2, h3) = b.add_h_twist(tt, n.unsigned_abs() as usize);

        b.connect(v0, h3);
        b.connect(v3, h0);
        b.connect(v1, h2);
        b.connect(v2, h1);
        b.build().unwrap()
    }

    // Blackboard-framed 2-cable: each crossing → a 2×2 block of 4 sub-crossings of the same type,
    // each edge → 2 parallel edges. Every component doubles into its two parallel copies
    // (n components → 2n; framing = the diagram's writhe per component).
    pub fn cable2(&self) -> Link {
        let (b, _) = self.cable2_builder();
        b.build().unwrap()
    }

    // The 2-cable in an open builder, plus `cab[i][slot] = (copy-0 port, copy-1 port)` so callers can
    // re-splice the cable (e.g. the Whitehead clasp) before building.
    fn cable2_builder(&self) -> (LinkBuilder, Vec<[(Port, Port); 4]>) {
        let mut b = LinkBuilder::new();

        // per crossing: 4 sub-crossings + the 4 internal edges; record the two cable ports at each slot
        let cab: Vec<[(Port, Port); 4]> = self.nodes().map(|x| {
            let t = x.node_type();
            let [s0, s1, s2, s3] = [b.add_node(t), b.add_node(t), b.add_node(t), b.add_node(t)];
            b.connect((s0, 2), (s1, 0));
            b.connect((s0, 3), (s2, 1));
            b.connect((s1, 3), (s3, 1));
            b.connect((s2, 2), (s3, 0));
            // cab[slot] = the slot's two cable ports in a consistent (CCW) order, so edge-joins line up.
            [((s2, 0), (s0, 0)), ((s0, 1), (s1, 1)), ((s1, 2), (s3, 2)), ((s3, 3), (s2, 3))]
        }).collect();

        // join the two cables of each original edge across its two endpoints. the two ends list their
        // ports in CCW order, which reverses along the edge, so copy 0 pairs with the other's copy 1.
        for e in self.edges() {
            let ((i, s), (j, t)) = self.edge_ends(e, false);
            let ((a0, a1), (b0, b1)) = (cab[i][s], cab[j][t]);
            b.connect(a0, b1);
            b.connect(a1, b0);
        }
        for _ in self.loops() {   // each free loop doubles
            b.add_loop();
            b.add_loop();
        }
        (b, cab)
    }

    // The `tw`-twisted Whitehead double D±(K), `tw` from the Seifert (0) framing (tw = 0 = untwisted,
    // trivial Alexander); `positive` = clasp sign. Seifert sits at 2·writhe blackboard half-twists.
    pub fn whitehead_double(&self, positive: bool, tw: i32) -> Link {
        self.whitehead_double_bbf(positive, 2 * self.writhe() + tw)
    }

    // Whitehead double with the framing counted from the blackboard framing (tw = 0 = the diagram's
    // blackboard 2-cable): cut the cable to a 4-end tangle, add `tw` half-twists, close with the clasp.
    pub fn whitehead_double_bbf(&self, positive: bool, tw: i32) -> Link {
        // cut a clean edge (joining two distinct crossings): frees the 4 cable ends
        let e0 = self.edges().into_iter()
            .find(|&e| {
                let ((i, _), (j, _)) = self.edge_ends(e, false);
                i != j
            })
            .expect("the companion needs an edge joining two distinct crossings");
        self.whitehead_double_impl(positive, 0, tw, e0, None).0
    }

    // Whitehead double cutting the cable at edge `e0` (must join two distinct crossings), placing
    // `tw_a` framing half-twists on one side of the cut and `tw_b` on the other. Also returns the
    // result-edges of `base`'s two doubled strands (empty if `base` is `None`) so the caller can
    // place a base point there.
    pub(crate) fn whitehead_double_impl(&self, positive: bool, tw_a: i32, tw_b: i32, e0: Edge, base: Option<Edge>) -> (Link, Vec<Edge>) {
        use crate::NodeType::{XL, XR};
        assert!(self.is_knot(), "the Whitehead double requires a knot companion");

        let (mut b, cab) = self.cable2_builder();

        // the cut edge's two ends: the a-side (node ia, slot sa) and b-side (ib, sb), with their
        // cable ports (a0, a1) / (b0, b1) in CCW order.
        let ((ia, sa), (ib, sb)) = self.edge_ends(e0, false);
        let ((a0, a1), (b0, b1)) = (cab[ia][sa], cab[ib][sb]);
        b.disconnect(a0);   // the swapped join means a0–b1, a1–b0 are removed
        b.disconnect(a1);

        let twist_type = |tw: i32| if tw >= 0 { XR } else { XL };

        // |tw_a| half-twists growing left off the a-side ends (a1, a0): attach the row's right
        // corners (SE, NE), continue from its left corners (SW, NW) — the 180°-rotation image of
        // the b-side below.
        let (ta0, ta1) = if tw_a == 0 {
            (a1, a0)
        } else {
            let (h0, h1, h2, h3) = b.add_h_twist(twist_type(tw_a), tw_a.unsigned_abs() as usize);
            b.connect(a1, h1);
            b.connect(a0, h2);
            (h0, h3)
        };

        // |tw_b| half-twists growing right off the b-side ends (b1, b0): attach the row's left
        // corners (NW, SW), continue from its right corners (NE, SE).
        let (tb0, tb1) = if tw_b == 0 {
            (b1, b0)
        } else {
            let (h0, h1, h2, h3) = b.add_h_twist(twist_type(tw_b), tw_b.unsigned_abs() as usize);
            b.connect(b1, h3);
            b.connect(b0, h0);
            (h2, h1)
        };

        // vertical clasp glued as in twist_knot: the twisted a-side ends on its left, the twisted
        // b-side ends on its right. positive = a positive clasp = D⁺ (verified Kh-identical to the
        // reference Wh⁺(4_1)).
        let ct = if positive { XL } else { XR };
        let (c0, c1, c2, c3) = b.add_v_twist(ct, 2);
        b.connect(c3, ta0);
        b.connect(c0, ta1);
        b.connect(c2, tb1);
        b.connect(c1, tb0);

        // `base` is not at the cut, so its cable join survives: the two ports of cab at `base` each
        // carry one of its doubled strands; read off their result-edge ids before consuming the builder.
        let base_edges: Vec<Edge> = base.into_iter().flat_map(|base| {
            let ((i, s), _) = self.edge_ends(base, false);
            let (p0, p1) = cab[i][s];
            [b.edge_at(p0).unwrap(), b.edge_at(p1).unwrap()]
        }).collect();

        (b.build().unwrap(), base_edges)
    }

    // The two (node, slot) ends of edge `e`. When `directed`, they are ordered as (tail, head)
    // along the orientation — the strand exits at the tail and enters at the head (cf.
    // `NodeOri::in_ports`); otherwise the order carries no meaning.
    fn edge_ends(&self, e: Edge, directed: bool) -> ((usize, usize), (usize, usize)) {
        assert!(!directed || self.is_oriented(), "directed edge_ends requires an oriented link");

        let (x, y) = self.nodes().enumerate().flat_map(|(i, n)|
            (0..4).filter(move |&s| n.edge(s) == e).map(move |s| (i, s))
        ).collect_tuple().unwrap_or_else(||
            panic!("edge {e} must appear exactly twice")
        );
        if !directed {
            return (x, y);
        }

        let is_in = |(i, s): (usize, usize)| {
            let ports = self.node(i).ori().in_ports().expect("directed edge_ends requires an oriented link");
            ports.contains(&s)
        };
        debug_assert!(is_in(x) != is_in(y), "edge {e} must have one head and one tail");
        if is_in(x) { (y, x) } else { (x, y) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Braid;
    use crate::misc::{jones_polynomial, det};

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
    fn twist_knot_determinants() {
        // twist knot K_n is a knot of determinant 2n+1 (3_1, 4_1, 5_2, … → 3, 5, 7, …); K_0 = unknot.
        // negative n is the mirror (twist_knot(-1-n) = mirror K_n), with the same determinant.
        for n in 0..=6 {
            let k = Link::twist_knot(n);
            assert_eq!(k.n_comps(), 1);
            assert!(k.is_oriented());
            assert_eq!(det(&k), 2 * n + 1, "det twist_knot({n})");
            assert_eq!(det(&Link::twist_knot(-1 - n)), 2 * n + 1, "det twist_knot({})", -1 - n);
        }
    }

    #[test]
    fn whitehead_double_of_unknot() {
        // the untwisted double of the unknot is the unknot, from any companion diagram (framing 2·writhe).
        for word in [vec![1, -2], vec![1, 2], vec![-1, -2]] {
            let u = Braid::from_iter(word).closure();
            for positive in [true, false] {
                let d = u.whitehead_double(positive, 0);
                assert_eq!(d.n_comps(), 1);
                assert_eq!(det(&d), 1, "D±(unknot) must be the unknot (trivial Alexander)");
            }
        }
    }

    #[test]
    fn whitehead_double_is_a_knot() {
        // D±(K) of a nontrivial companion builds (exercising the planarity check) and is a knot with
        // 4·(crossings) + 2·|writhe| (framing) + 2 (clasp) crossings. det would confirm untwisted but
        // is exponential here — see _of_unknot.
        for name in ["3_1", "4_1", "5_2", "6_1"] {
            let k = Link::test_data(name);
            let expected = 4 * k.n_crossings() + 2 * k.writhe().unsigned_abs() as usize + 2;
            for positive in [true, false] {
                let d = k.whitehead_double(positive, 0);
                assert_eq!(d.n_comps(), 1, "D±({name}) is a knot");
                assert_eq!(d.n_crossings(), expected, "D±({name}) crossing count");
                assert!(d.is_oriented());
            }
        }
    }

    #[test]
    fn cable2_trefoil() {
        let k = Link::test_data("3_1");
        let c = k.cable2();
        assert_eq!(c.n_crossings(), 4 * k.n_crossings());
        assert_eq!(c.n_comps(), 2, "2-cable of a knot is a 2-component link");
        assert!(c.is_oriented());
        let _ = c.seifert_circles(); // exercises orientation consistency
    }

    #[test]
    fn cable2_is_hopf() {
        // 2-cable of a ±1-framed unknot is the Hopf link (linking ±1), not the 2-component unlink.
        let c = Link::test_data("unknot_l_twist").cable2();
        assert_eq!(c.n_comps(), 2);
        assert_ne!(jones_polynomial(&c), jones_polynomial(&Link::unlink(2)),
            "2-cable of a framed unknot must be linked (Hopf), not the unlink");
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
