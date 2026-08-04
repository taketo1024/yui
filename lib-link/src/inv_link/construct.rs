//! Constructions specific to [`InvLink`] — the equivariant counterparts of the operations in
//! [`crate::link::construct`]. Each one produces a diagram carrying a strong inversion, recovered by
//! reindexing to the standard involution `e ↦ (n+1-e)%n+1`.

use num_integer::Integer;

use crate::{Edge, InvLink, Link};

impl InvLink {
    // Equivariant connected sum at the two (on-axis) base points.
    pub fn conn_sum(&self, other: &InvLink) -> InvLink {
        let self_e = self.base_pt().expect("self needs a base point");
        let other_e = other.base_pt().expect("other needs a base point");
        self.conn_sum_at(other, self_e, other_e)
    }

    // Equivariant connected sum: splice along on-axis edges (`inv_edge(e) == e`) of each summand,
    // then recover the combined strong inversion by reindexing to the standard involution.
    pub fn conn_sum_at(&self, other: &InvLink, self_e: Edge, other_e: Edge) -> InvLink {
        assert_eq!(self.inv_edge(self_e), self_e, "self_e {self_e} must be on-axis");
        assert_eq!(other.inv_edge(other_e), other_e, "other_e {other_e} must be on-axis");

        let inner = self.inner().conn_sum_at(other.inner(), self_e, other_e);
        let starts = inner.edges();
        Self::from_standard_reindex(inner, starts).expect("connected sum is not τ-symmetric")
    }

    /// The 3-pretzel `P(a, b, a)` with its strong inversion (the π-rotation through the middle
    /// band), reindexed to the standard involution. Standard convention: all-odd 3-pretzel bands
    /// are anti-parallel, so a positive (right-handed) half-twist is a NEGATIVE crossing —
    /// `writhe(P(a, b, a)) = -(2a + b)`.
    pub fn sym_pretzel(a: i32, b: i32, c: i32) -> InvLink {
        assert_eq!(a, c, "the strong inversion needs P(a, b, a)");
        assert!(a % 2 != 0 && b % 2 != 0, "all-odd parameters required");

        // the two spanning arcs are τ-fixed; `Link::pretzel` numbers from the top one, which is what
        // the standard involution expects.
        let link = Link::pretzel(a, b, c);
        InvLink::from_symmetric_pd_code(link.pd_code())
    }

    // Strongly-invertible Whitehead double of a symmetric companion. The 2-cable inherits τ; the
    // clasp and `tw` framing twists go in at the *other* on-axis edge (`inv_edge(e) == e`,
    // e ≠ base_pt), split evenly across the axis so the diagram stays τ-invariant. `tw` counts from
    // the Seifert framing and must be even. The base point lands on the doubled on-axis strand.
    pub fn whitehead_double(&self, positive: bool, tw: i32) -> InvLink {
        assert!(tw.is_even(), "tw must be even for a τ-symmetric diagram");

        let base = self.base_pt().expect("companion needs a base point");
        assert_eq!(self.inv_edge(base), base, "base point must be on-axis");
        let cut = self.edges().into_iter()
            .find(|&e| e != base && self.inv_edge(e) == e)
            .expect("need a second on-axis edge for the clasp");

        let half = self.writhe() + tw / 2;   // (2·writhe + tw) / 2 = half the blackboard framing
        let (inner, base_edges) = Link::whitehead_double_impl(self.inner(), positive, half, half, cut, Some(base));

        // base point on the on-axis doubled base_pt strand: reindex from the copy that realizes the
        // standard τ, so edge 1 lands there rather than at the clasp.
        Self::from_standard_reindex(inner, base_edges)
            .expect("whitehead double diagram is not τ-symmetric at the base point")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::misc::det;

    #[test]
    fn conn_sum_is_equivariant() {
        // construction succeeding ⟺ from_standard_reindex found a valid τ on the sum.
        let k1 = InvLink::test_data("3_1");
        let k2 = InvLink::test_data("4_1");
        let cs = k1.conn_sum(&k2);
        assert!(cs.is_knot());
        assert_eq!(det(cs.inner()), 3 * 5, "det is multiplicative under conn sum");
    }

    #[test]
    fn sym_pretzel_is_symmetric() {
        // construction succeeding ⟺ the standard involution is realized by the pretzel numbering.
        let k = InvLink::sym_pretzel(-3, 3, -3);
        assert!(k.is_knot());
        assert_eq!(k.n_crossings(), 9);
        assert_eq!(det(k.inner()), 9);  // |ab + bc + ca|
        assert_eq!(k.writhe(), 3);      // = -(2a + b) for anti-parallel bands
    }

    #[test]
    fn whitehead_double_is_symmetric() {
        // building succeeding ⟺ try_new found an on-axis reindex start ⟺ the diagram is τ-symmetric.
        for name in ["3_1", "4_1"] {
            let k = InvLink::test_data(name);
            let w = k.whitehead_double(true, 0);
            assert!(w.is_knot());
            assert_eq!(w.base_pt(), Some(1), "base point on the doubled base_pt strand");
            assert_eq!(w.inv_edge(1), 1, "base point is on-axis");
            // 4·n_crossings (cable) + |blackboard framing| (split each side) + 2 (clasp)
            let bl = (2 * k.writhe()).unsigned_abs() as usize;
            assert_eq!(w.n_crossings(), 4 * k.n_crossings() + bl + 2);
        }
    }
}
