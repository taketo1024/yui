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

    // The 3-pretzel P(a, b, c): three vertical twist regions of |a|, |b|, |c| half-twists, chained by
    // their inner arcs and closed by one spanning arc on the top and one on the bottom. A positive
    // parameter is a right-handed (XR) region; all three must be nonzero. Edge 1 is the top spanning
    // arc.
    pub fn pretzel(a: i32, b: i32, c: i32) -> Link {
        use crate::NodeType::{XL, XR};
        assert!(a != 0 && b != 0 && c != 0, "pretzel parameters must be nonzero");

        let mut bld = LinkBuilder::new();
        let ends: Vec<_> = [a, b, c].iter().map(|&v| {
            let ty = if v > 0 { XR } else { XL };
            bld.add_v_twist(ty, v.unsigned_abs() as usize) // (sw, se, ne, nw)
        }).collect();
        let (sw1, se1, ne1, nw1) = ends[0];
        let (sw2, se2, ne2, nw2) = ends[1];
        let (sw3, se3, ne3, nw3) = ends[2];

        bld.connect(ne1, nw2);
        bld.connect(ne2, nw3);
        bld.connect(se1, sw2);
        bld.connect(se2, sw3);
        bld.connect(nw1, ne3); // spanning top
        bld.connect(sw1, se3); // spanning bottom

        let start = bld.edge_at(nw1).unwrap();
        bld.build().expect("pretzel must be planar").reindexed(start, 1)
    }

    // Blackboard-framed 2-cable: each crossing → a 2×2 block of 4 sub-crossings of the same type,
    // each edge → 2 parallel edges. Every component doubles into its two parallel copies
    // (n components → 2n; framing = the diagram's writhe per component).
    pub fn cable2(l: &Link) -> Link {
        let (b, _) = Self::cable2_builder(l);
        b.build().unwrap()
    }

    // The 2-cable in an open builder, plus `cab[i][slot.index()] = (copy-0 port, copy-1 port)` so callers can
    // re-splice the cable (e.g. the Whitehead clasp) before building.
    fn cable2_builder(l: &Link) -> (LinkBuilder, Vec<[(Port, Port); 4]>) {
        let mut b = LinkBuilder::new();

        // per crossing: 4 sub-crossings + the 4 internal edges; record the two cable ports at each slot
        let cab: Vec<[(Port, Port); 4]> = l.nodes().map(|x| {
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
        for e in l.edges() {
            let ((i, s), (j, t)) = l.edge_ends(e, false);
            let ((a0, a1), (b0, b1)) = (cab[i][s.index()], cab[j][t.index()]);
            b.connect(a0, b1);
            b.connect(a1, b0);
        }
        for _ in l.loops() {   // each free loop doubles
            b.add_loop();
            b.add_loop();
        }
        (b, cab)
    }

    // The `tw`-twisted Whitehead double D±(K), `tw` from the Seifert (0) framing (tw = 0 = untwisted,
    // trivial Alexander); `positive` = clasp sign. Seifert sits at 2·writhe blackboard half-twists.
    pub fn whitehead_double(l: &Link, positive: bool, tw: i32) -> Link {
        Self::whitehead_double_bbf(l, positive, 2 * l.writhe() + tw)
    }

    // Whitehead double with the framing counted from the blackboard framing (tw = 0 = the diagram's
    // blackboard 2-cable): cut the cable to a 4-end tangle, add `tw` half-twists, close with the clasp.
    pub fn whitehead_double_bbf(l: &Link, positive: bool, tw: i32) -> Link {
        // cut a clean edge (joining two distinct crossings): frees the 4 cable ends
        let cut = l.edges().into_iter()
            .find(|&e| {
                let ((i, _), (j, _)) = l.edge_ends(e, false);
                i != j
            })
            .expect("the companion needs an edge joining two distinct crossings");
        Self::whitehead_double_impl(l, positive, 0, tw, cut, None)
    }

    // Whitehead double cutting the cable at edge `cut` (must join two distinct crossings), placing
    // `tw_a` framing half-twists on one side of the cut and `tw_b` on the other. The result is based
    // at one of `base`'s two doubled strands.
    pub(crate) fn whitehead_double_impl(l: &Link, positive: bool, tw_a: i32, tw_b: i32, cut: Edge, base: Option<Edge>) -> Link {
        use crate::NodeType::{XL, XR};
        assert!(l.is_knot(), "the Whitehead double requires a knot companion");
        assert!(base.is_none_or(|b| b != cut), "the base point must be away from the cut");

        let (mut b, cab) = Self::cable2_builder(l);

        // the cut edge's two ends: the a-side (node ia, slot sa) and b-side (ib, sb), with their
        // cable ports (a0, a1) / (b0, b1) in CCW order.
        let ((ia, sa), (ib, sb)) = l.edge_ends(cut, false);
        let ((a0, a1), (b0, b1)) = (cab[ia][sa.index()], cab[ib][sb.index()]);
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

        // `base` is away from the cut, so its cable join survives: `cab` holds the two ports its
        // doubled strands run through. Read the id off before `build` consumes the builder.
        let base_pt = base.map(|e| {
            let ((i, s), _) = l.edge_ends(e, false);
            b.edge_at(cab[i][s.index()].0).unwrap()
        });

        let double = b.build().unwrap();
        match base_pt {
            Some(e) => double.with_base_pt(e),
            None => double,
        }
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
    fn pretzel_determinants() {
        // det P(a, b, c) = |ab + bc + ca| — includes the (-2, 3, 7)-pretzel (det 1). `det` sums over
        // all 2^n resolutions, so the 15-crossing cases are left to `pretzel_band_symmetries`.
        for (a, b, c) in [(1, 1, 1), (-1, -1, -1), (1, 3, 5), (-2, 3, 7)] {
            let l = Link::pretzel(a, b, c);
            let n = (a.unsigned_abs() + b.unsigned_abs() + c.unsigned_abs()) as usize;
            assert_eq!(l.n_crossings(), n, "P({a},{b},{c}) crossing count");
            assert_eq!(det(&l), (a * b + b * c + c * a).abs(), "det P({a},{b},{c})");
        }
    }

    #[test]
    fn pretzel_band_symmetries() {
        // the three bands sit in a cycle, so rotating them — or reversing their order — leaves the
        // diagram itself unchanged, not merely the knot type.
        // Canonical form of the *unoriented* diagram: least PD code over all start edges and both
        // strand directions. Rotating or reversing the bands may reverse the direction (it does
        // whenever a band is even), so the orientation must be quotiented out here.
        let canon = |l: &Link| {
            // reverse the strand: the under-strand enters at the far end, CCW order unchanged.
            let rev = Link::from_pd_code(l.pd_code().into_iter().map(|[a, b, c, d]| [c, d, a, b]));
            l.reindexed_canon().pd_code().min(rev.reindexed_canon().pd_code())
        };

        for (a, b, c) in [(1, 3, 5), (3, 5, 7), (-3, 3, -3), (-2, 3, 7), (-5, 5, -5)] {
            let p = canon(&Link::pretzel(a, b, c));
            assert_eq!(canon(&Link::pretzel(b, c, a)), p, "cyclic P({a},{b},{c})");
            assert_eq!(canon(&Link::pretzel(c, a, b)), p, "cyclic P({a},{b},{c})");
            assert_eq!(canon(&Link::pretzel(c, b, a)), p, "reversal P({a},{b},{c})");
        }
    }

    #[test]
    fn mirror_identities() {
        // mirroring flips every band, every cable crossing and the clasp.
        for (a, b, c) in [(1, 3, 5), (3, 5, 7), (-2, 3, 7), (-5, 5, -5)] {
            assert_eq!(Link::pretzel(a, b, c).mirror().pd_code(), Link::pretzel(-a, -b, -c).pd_code(),
                "mirror P({a},{b},{c})");
        }
        for name in ["3_1", "4_1", "5_2"] {
            let k = Link::test_data(name);
            assert_eq!(Link::cable2(&k).mirror().pd_code(), Link::cable2(&k.mirror()).pd_code(),
                "mirror cable2({name})");
            assert_eq!(Link::whitehead_double(&k, true, 0).mirror().pd_code(),
                       Link::whitehead_double(&k.mirror(), false, 0).pd_code(),
                "mirror D+({name}) vs D-(mirror {name})");
        }
    }

    #[test]
    fn pretzel_trefoil() {
        assert!(same_knot(&Link::pretzel(1, 1, 1), &Link::test_data("3_1")));
    }

    #[test]
    #[ignore = "slow: det sums over all 2^n resolutions, and the framing test needs 14-crossing doubles"]
    fn whitehead_double_of_unknot() {
        // the untwisted double of the unknot is the unknot, from any companion diagram (framing 2·writhe).
        for word in [vec![1, -2], vec![1, 2], vec![-1, -2]] {
            let u = Braid::from_iter(word).closure();
            for positive in [true, false] {
                let d = Link::whitehead_double(&u, positive, 0);
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
                let d = Link::whitehead_double(&k, positive, 0);
                assert_eq!(d.n_comps(), 1, "D±({name}) is a knot");
                assert_eq!(d.n_crossings(), expected, "D±({name}) crossing count");
                assert!(d.is_oriented());
            }
        }
    }

    #[test]
    fn whitehead_double_clasp_sign() {
        // twisted doubles of the unknot are twist knots, whose determinants separate the clasp
        // signs; values pinned against the Kh-verified reference implementation. `tw` counts
        // half-twists, and det is 2*tw ∓ 1 — two values are enough to pin that line, and the
        // smallest two keep the 2^n determinant cheap.
        let u = Braid::from([1, -2]).closure(); // writhe-0 unknot diagram
        for (tw, pos_det, neg_det) in [(1, 1, 3), (2, 3, 5)] {
            assert_eq!(det(&Link::whitehead_double(&u, true, tw)), pos_det, "D+(U, tw={tw})");
            assert_eq!(det(&Link::whitehead_double(&u, false, tw)), neg_det, "D-(U, tw={tw})");
        }
    }

    #[test]
    fn cable2_trefoil() {
        let k = Link::test_data("3_1");
        let c = Link::cable2(&k);
        assert_eq!(c.n_crossings(), 4 * k.n_crossings());
        assert_eq!(c.n_comps(), 2, "2-cable of a knot is a 2-component link");
        assert!(c.is_oriented());
        let _ = c.seifert_circles(); // exercises orientation consistency
    }

    #[test]
    fn cable2_is_hopf() {
        // 2-cable of a ±1-framed unknot is the Hopf link (linking ±1), not the 2-component unlink.
        let c = Link::cable2(&Link::test_data("unknot_l_twist"));
        assert_eq!(c.n_comps(), 2);
        assert_ne!(jones_polynomial(&c), jones_polynomial(&Link::unlink(2)),
            "2-cable of a framed unknot must be linked (Hopf), not the unlink");
    }

}
