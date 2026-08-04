//! Constructions producing new links from patterns — twist knots, cables and satellites.
//! Contrast with [`crate::link::link_ops`], which operates on links you already have.

use crate::{Link, Edge, LinkBuilder, Port};

impl Link {
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

        // The insert replaces the parallel cable strands (CCW at the cut ⇒ (a0, a1) = (lower, upper),
        // (b0, b1) = (upper, lower)):
        //
        //   a1 ────[       ]────[       ]────[       ]──── b0   (upper strand)
        //          [ row_a ]    [ clasp ]    [ row_b ]
        //   a0 ────[       ]────[       ]────[       ]──── b1   (lower strand)
        //
        // The clasp turns each side's pair back on itself and the turn-backs hook (winding 0).
        // Corners are destructured by picture position (ul/ur/ll/lr = upper/lower × left/right);
        // the wiring is identical for every piece: enter at (ul, ll), continue from (ur, lr).
        let twist_type = |tw: i32| if tw >= 0 { XR } else { XL };
        let (mut upper, mut lower) = (a1, a0);

        if tw_a != 0 {
            let (ll, lr, ur, ul) = b.add_h_twist(twist_type(tw_a), tw_a.unsigned_abs() as usize);
            b.connect(upper, ul);
            b.connect(lower, ll);
            (upper, lower) = (ur, lr);
        }

        // positive = a positive clasp = D⁺ (pinned by the clasp-sign determinant test).
        let ct = if positive { XL } else { XR };
        let (ll, lr, ur, ul) = b.add_v_twist(ct, 2);
        b.connect(upper, ul);
        b.connect(lower, ll);
        (upper, lower) = (ur, lr);

        if tw_b != 0 {
            let (ll, lr, ur, ul) = b.add_h_twist(twist_type(tw_b), tw_b.unsigned_abs() as usize);
            b.connect(upper, ul);
            b.connect(lower, ll);
            (upper, lower) = (ur, lr);
        }

        b.connect(upper, b0);
        b.connect(lower, b1);

        // `base` is not at the cut, so its cable join survives: the two ports of cab at `base` each
        // carry one of its doubled strands; read off their result-edge ids before consuming the builder.
        let base_edges: Vec<Edge> = base.into_iter().flat_map(|base| {
            let ((i, s), _) = self.edge_ends(base, false);
            let (p0, p1) = cab[i][s];
            [b.edge_at(p0).unwrap(), b.edge_at(p1).unwrap()]
        }).collect();

        (b.build().unwrap(), base_edges)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Braid;
    use crate::misc::{jones_polynomial, det};

    fn same_knot(a: &Link, b: &Link) -> bool {
        let ja = jones_polynomial(a);
        ja == jones_polynomial(b) || ja == jones_polynomial(&b.mirror())
    }

    #[test]
    fn braid_closures_are_jones_faithful() {
        assert!(same_knot(&Braid::from([1, 1, 1]).closure(), &Link::test_data("3_1")));
        assert!(same_knot(&Braid::from([1, -2, 1, -2]).closure(), &Link::test_data("4_1")));
    }

    #[test]
    fn twist_knot_matches_references() {
        // twist_knot(n) vs the independent PD-code reference, by Jones (mirror allowed);
        // twist_knot(-1-n) is the mirror.
        let table = ["3_1", "4_1", "5_2", "6_1", "7_2"];
        for (i, name) in table.iter().enumerate() {
            let n = i as i32 + 1;
            assert!(same_knot(&Link::twist_knot(n), &Link::test_data(name)), "twist_knot({n}) ≠ {name}");
            assert!(same_knot(&Link::twist_knot(-1 - n), &Link::test_data(name)), "twist_knot({}) ≠ mirror {name}", -1 - n);
        }
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
    fn whitehead_double_clasp_sign() {
        // twisted doubles of the unknot are twist knots, whose determinants separate the clasp
        // signs; values pinned against the Kh-verified reference implementation.
        let u = Braid::from([1, -2]).closure(); // writhe-0 unknot diagram
        for (tw, pos_det, neg_det) in [(2, 3, 5), (4, 7, 9)] {
            assert_eq!(det(&u.whitehead_double(true, tw)), pos_det, "D+(U, tw={tw})");
            assert_eq!(det(&u.whitehead_double(false, tw)), neg_det, "D-(U, tw={tw})");
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

}
