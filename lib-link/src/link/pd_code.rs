//! PD (planar diagram) codes — the KnotAtlas `X[i,j,k,l]` presentation of a link, and loading a
//! link from the data directory by name.
//!
//! A PD code is a sequence of crossings of the form:
//!
//! ```text
//!     d   c
//!      \ /
//!       \     = [a, b, c, d]
//!      / \
//!     a   b
//! ```
//!
//! The lower edge is always oriented a -> c.
//! see: <http://katlas.math.toronto.edu/wiki/Planar_Diagrams>

use itertools::Itertools;

use super::{Edge, Link, Node, NodeType};

/// One crossing of a PD code: `X[i,j,k,l]`, listed counter-clockwise from the incoming under-strand.
pub type PDCodeX = [Edge; 4];

impl Link {
    pub fn from_pd_code<I>(pd_code: I) -> Self
    where I: IntoIterator<Item = PDCodeX> {
        let nodes = pd_code.into_iter().map(Node::from_pd_code).collect_vec();
        let mut l = Self::from_nodes(nodes); // unoriented
        l.reorient(|_, s| s.index() == 0); // PD convention: the under-strand enters at slot 0.
        l
    }

    /// The KnotAtlas PD code — one `X[i,j,k,l]` per crossing, CCW from the incoming under-strand `i`.
    /// Built by traversing the components and emitting each crossing at its under-pass, so `i` is the
    /// under-strand's incoming edge; round-trips through `from_pd_code`. Resolved (V/H) nodes are skipped.
    pub fn pd_code(&self) -> Vec<PDCodeX> {
        let mut pd = Vec::with_capacity(self.n_nodes());
        self.traverse_comps(|_, i, j| {
            let x = self.node(i);
            // under-pass: `XL`'s under-strand enters at an even port (0/2), `XR`'s at an odd port (1/3);
            // the other pass is the over-strand and is skipped, so each crossing is emitted exactly once.
            let under = match x.node_type() {
                NodeType::XL => j.index() % 2 == 0,
                NodeType::XR => j.index() % 2 == 1,
                _ => return,
            };
            if under {
                pd.push([x.edge(j), x.edge(j.shift(1)), x.edge(j.shift(2)), x.edge(j.shift(3))]);
            }
        });
        pd
    }

    pub fn load(name: &str) -> Result<Link, Box<dyn std::error::Error>> {
        // the data dir is external and empty on a fresh checkout; tests must not depend on it.
        assert!(!cfg!(feature = "test-utils"), "`load` reads the data directory — use `test_data` in tests");
        let json = yui_core::util::data_dir::load_json("links", name)?;
        let data: Vec<PDCodeX> = serde_json::from_str(&json)?;
        Ok(Link::from_pd_code(data))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::NodeType::XL;

    #[test]
    fn pd_code_roundtrip() {
        use crate::misc::jones_polynomial;

        let l = Link::test_data("3_1"); // chiral trefoil (oriented)

        // XL round-trip: pd_code -> from_pd_code recovers the same link.
        let l2 = Link::from_pd_code(l.pd_code());
        assert_eq!(jones_polynomial(&l), jones_polynomial(&l2), "XL round-trip");

        // XR round-trip: `mirror` flips XL<->XR, exercising the rotated PD emission.
        let m = l.mirror();
        let m2 = Link::from_pd_code(m.pd_code());
        assert_eq!(jones_polynomial(&m), jones_polynomial(&m2), "XR round-trip");

        // chirality guard: 3_1 differs from its mirror, so the XR case isn't vacuous.
        assert_ne!(jones_polynomial(&l), jones_polynomial(&m));
    }

    #[test]
    fn pd_code_depends_only_on_the_diagram() {
        // Rebuilding changes the node order, so this holds only if traversal follows the orientation.
        for name in ["3_1", "4_1", "5_2", "6_1", "L2a1", "L4a1"] {
            let l = Link::test_data(name);
            assert_eq!(Link::from_pd_code(l.pd_code()).pd_code(), l.pd_code(), "{name}");
        }
        let l = Link::pretzel(1, 3, 5);
        assert_eq!(Link::from_pd_code(l.pd_code()).pd_code(), l.pd_code(), "pretzel(1,3,5)");
    }

    #[test]
    fn link_from_pd_code() {
        let l = Link::test_data("unknot_l_twist");
        assert_eq!(l.n_nodes(), 1);
        assert_eq!(l.node(0).node_type(), XL);
    }
}
