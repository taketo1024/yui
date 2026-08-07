//! [`Node`]: one vertex of a diagram — a crossing (`XL`/`XR`) or a smoothing
//! (`V`/`H`) — its four incident edges, and the orientation given by the pair of
//! [`Slot`]s the two strands enter by.
//!
//! The four [`NodeType`]s, with the incident edges numbered counter-clockwise
//! from the lower left (`0 = SW`, `1 = SE`, `2 = NE`, `3 = NW`):
//!
//! ```text
//!     3   2         3   2         3   2         3   2
//!      \ /           \ /           \ /           \_/
//!       \    = XL,    /    = XR,   | |   = V,     _    = H
//!      / \           / \           / \           / \
//!     0   1         0   1         0   1         0   1
//! ```

use std::fmt::Display;

use yui_core::bitseq::Bit;
use yui_core::num::Sign;
use yui_core::ext::CloneAnd;

use crate::Path;
use super::Edge;

use NodeType::{XL, XR, V, H};

/// One of a node's four ends, counter-clockwise from the lower left
/// (see the module diagram).
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, derive_more::Display, Debug)]
pub enum Slot {
    SW, SE, NE, NW
}

impl Slot {
    pub const ALL: [Slot; 4] = [Slot::SW, Slot::SE, Slot::NE, Slot::NW];

    pub fn index(&self) -> usize {
        *self as usize
    }

    // the slot `k` steps counter-clockwise from this one.
    pub fn shift(&self, k: usize) -> Self {
        Self::from((self.index() + k) % 4)
    }
}

impl From<usize> for Slot {
    fn from(i: usize) -> Self {
        Self::ALL.get(i).copied().unwrap_or_else(|| panic!("slot {i} out of range"))
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, derive_more::Display, Debug)]
pub enum NodeType {
    XL, XR, V, H
}

impl NodeType {
    // The slot at the other end of the strand passing through `slot` — the strand pairing.
    pub fn paired_slot(&self, slot: Slot) -> Slot {
        let i = slot.index();
        Slot::from(match self {
            XL | XR => (i + 2) % 4,   // SW<->NE, SE<->NW
            V => 3 - i,               // SW<->NW, SE<->NE
            H => i ^ 1,               // SW<->SE, NE<->NW
        })
    }

    pub fn mirror(&self) -> Self {
        match self {
            XL => XR,
            XR => XL,
            _  => *self
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Node {
    node_type: NodeType,
    // the two slots where the strands enter, sorted; `None` when the node is not coherently
    // oriented. The two must lie on different strands — see `Node::orientable`.
    incoming: Option<(Slot, Slot)>,
    edges: [Edge; 4],
}

impl Node {
    pub fn new(node_type: NodeType, incoming: Option<(Slot, Slot)>, edges: [Edge; 4]) -> Self {
        let incoming = incoming.map(|(p, q)| (p.min(q), p.max(q)));
        assert!(
            incoming.is_none_or(|(p, q)| Self::orientable(node_type, p, q)),
            "the two incoming slots must lie on different strands of {node_type}: {incoming:?}"
        );
        Node { node_type, edges, incoming }
    }

    pub fn from_pd_code(edges: [Edge; 4]) -> Self {
        Node::new(NodeType::XL, None, edges)
    }

    // A node is coherently oriented only when its two incoming slots sit on different strands;
    // `paired_slot` is the strand pairing.
    pub fn orientable(node_type: NodeType, p: Slot, q: Slot) -> bool {
        p != q && node_type.paired_slot(p) != q
    }

    pub fn node_type(&self) -> NodeType {
        self.node_type
    }

    pub fn edge(&self, s: Slot) -> Edge {
        self.edges[s.index()]
    }

    pub fn edges(&self) -> &[Edge; 4] {
        &self.edges
    }

    pub fn min_edge(&self) -> Edge {
        *self.edges.iter().min().unwrap()
    }

    pub fn is_crossing(&self) -> bool {
        matches!(self.node_type, XL | XR)
    }

    pub fn is_resolved(&self) -> bool {
        matches!(self.node_type, V | H)
    }

    pub fn resolve(&self, r: Bit) -> Self {
        use Bit::{Bit0, Bit1};

        self.clone_and(|x| {
            x.node_type = match (x.node_type, r) {
                (XL, Bit0) | (XR, Bit1) => H,
                (XL, Bit1) | (XR, Bit0) => V,
                _ => panic!("cannot resolve node-type: {}", x.node_type)
            };

            // the smoothing re-pairs the strands, so an orientation survives only one of the two:
            // incoming {0,1}/{2,3} survives V, incoming {1,2}/{0,3} survives H.
            if x.incoming.is_some_and(|(p, q)| !Self::orientable(x.node_type, p, q)) {
                x.incoming = None
            }
        })
    }

    pub fn is_oriented(&self) -> bool {
        self.incoming.is_some()
    }

    // The two slots where the strands enter, sorted.
    pub fn incoming(&self) -> Option<(Slot, Slot)> {
        self.incoming
    }

    // False on an unoriented node, where no slot is claimed either way.
    pub fn is_incoming(&self, s: Slot) -> bool {
        self.incoming.is_some_and(|(p, q)| p == s || q == s)
    }

    // Goes through `new`, so the pair is sorted and validated however the caller passes it.
    pub(crate) fn set_incoming(&mut self, incoming: Option<(Slot, Slot)>) {
        *self = Self::new(self.node_type, incoming, self.edges);
    }

    pub fn is_pos(&self) -> bool {
        self.sign().map(|x| x.is_positive()).unwrap_or(false)
    }

    pub fn is_neg(&self) -> bool {
        self.sign().map(|x| x.is_negative()).unwrap_or(false)
    }

    // Only a crossing has a sign, read off the type and which pair of slots the strands enter by.
    pub fn sign(&self) -> Option<Sign> {
        use Slot::{SW, SE, NE, NW};
        match (self.node_type, self.incoming?) {
            (XL, (SE, NE)) | (XL, (SW, NW)) | (XR, (SW, SE)) | (XR, (NE, NW)) => Some(Sign::Pos),
            (XL, (SW, SE)) | (XL, (NE, NW)) | (XR, (SE, NE)) | (XR, (SW, NW)) => Some(Sign::Neg),
            _ => None,
        }
    }

    pub fn mirror(&self) -> Self {
        self.clone_and(|x|
            x.node_type = self.node_type.mirror()
        )
    }

    // Reverse both strands: each now enters by the slot it used to leave. The two incoming slots
    // stay on different strands, so orientability — and the sign — survive.
    pub fn reversed(&self) -> Self {
        self.clone_and(|x|
            x.incoming = self.incoming.map(|(p, q)| {
                let (p, q) = (self.paired_slot(p), self.paired_slot(q));
                (p.min(q), p.max(q))
            })
        )
    }

    pub fn is_adj_to(&self, x: &Node) -> bool {
        self.edges.iter().any(|e| x.edges.contains(e))
    }

    pub fn arcs(&self) -> (Path, Path) {
        let comp = |i: usize, j: usize| {
            let (ei, ej) = (self.edges[i], self.edges[j]);
            if ei == ej {
                Path::circ([ei])
            } else {
                Path::arc([ei, ej])
            }
        };
        match self.node_type {
            XL |
            XR => (comp(0, 2), comp(1, 3)),
            V  => (comp(0, 3), comp(1, 2)),
            H  => (comp(0, 1), comp(2, 3))
        }
    }

    pub fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge {
        Self {
            node_type: self.node_type,
            incoming: self.incoming,
            edges: self.edges.map(f)
        }
    }

    pub fn paired_slot(&self, s: Slot) -> Slot {
        self.node_type.paired_slot(s)
    }
}

impl Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{:?}", self.node_type, self.edges)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // the four orientations a crossing can carry, by the direction the strands run.
    const UP:    Option<(Slot, Slot)> = Some((Slot::SW, Slot::SE));
    const LEFT:  Option<(Slot, Slot)> = Some((Slot::SE, Slot::NE));
    const DOWN:  Option<(Slot, Slot)> = Some((Slot::NE, Slot::NW));
    const RIGHT: Option<(Slot, Slot)> = Some((Slot::SW, Slot::NW));

    fn node(ntype: NodeType, incoming: Option<(Slot, Slot)>) -> Node {
        Node::new(ntype, incoming, [0, 1, 2, 3])
    }

    #[test]
    fn test_is_resolved() {
        for ntype in [XL, XR] {
            let c = node(ntype, None);
            assert!(c.is_crossing());
            assert!(!c.is_resolved());
        }
        for ntype in [V, H] {
            let c = node(ntype, None);
            assert!(!c.is_crossing());
            assert!(c.is_resolved());
        }
    }

    #[test]
    fn test_resolve() {
        use Bit::{Bit0, Bit1};

        // (ntype, incoming, bit, expected_ntype, expected_incoming)
        let cases = [
            (XL, None, Bit0, H, None),
            (XL, None, Bit1, V, None),
            (XR, None, Bit0, V, None),
            (XR, None, Bit1, H, None),
            // an orientation survives exactly one of the two smoothings: the one that leaves the
            // two incoming slots on different strands.
            (XL, UP,   Bit1, V, UP),
            (XL, LEFT, Bit0, H, LEFT),
            (XL, UP,   Bit0, H, None),
            (XL, LEFT, Bit1, V, None),
            (XR, LEFT, Bit0, V, None),
            (XR, UP,   Bit1, H, None),
        ];

        for (ntype, incoming, bit, expected_ntype, expected) in cases {
            let c = node(ntype, incoming).resolve(bit);
            assert!(c.is_resolved());
            assert_eq!(c.node_type(), expected_ntype);
            assert_eq!(c.incoming(), expected);
        }
    }

    #[test]
    fn test_anti_parallel() {
        // A resolved node may have its two strands running opposite ways — the diagonal slot pairs,
        // which no crossing can carry. `resolve` never produces them, but a V/H diagram oriented in
        // its own right does (e.g. orienting the circles of a resolution).
        for ntype in [V, H] {
            for (a, b) in [(Slot::SW, Slot::NE), (Slot::SE, Slot::NW)] {
                let pair = (a, b);
                assert!(Node::orientable(ntype, pair.0, pair.1));
                assert_eq!(node(ntype, Some(pair)).incoming(), Some(pair));
            }
        }
        for ntype in [XL, XR] {
            assert!(!Node::orientable(ntype, Slot::SW, Slot::NE), "a crossing's diagonal is one strand");
            assert!(!Node::orientable(ntype, Slot::SE, Slot::NW));
        }
    }

    #[test]
    fn test_mirror() {
        // mirror flips XL <-> XR and preserves the orientation.
        let cases = [
            (XL, None,  XR, None),
            (XR, None,  XL, None),
            (H,  None,  H,  None),
            (V,  None,  V,  None),
            (XL, UP,    XR, UP),
            (XR, LEFT,  XL, LEFT),
            (V,  DOWN,  V,  DOWN),
            (H,  RIGHT, H,  RIGHT),
        ];

        for (ntype, incoming, expected_ntype, expected) in cases {
            let c = node(ntype, incoming).mirror();
            assert_eq!(c.node_type(), expected_ntype);
            assert_eq!(c.incoming(), expected);
        }
    }

    #[test]
    fn test_sign() {
        let cases = [
            (XL, LEFT,  Some(Sign::Pos)),
            (XL, RIGHT, Some(Sign::Pos)),
            (XR, UP,    Some(Sign::Pos)),
            (XR, DOWN,  Some(Sign::Pos)),
            (XL, UP,    Some(Sign::Neg)),
            (XL, DOWN,  Some(Sign::Neg)),
            (XR, LEFT,  Some(Sign::Neg)),
            (XR, RIGHT, Some(Sign::Neg)),
            // unsigned: no orientation, or a resolved node.
            (XL, None,  Option::None),
            (XR, None,  Option::None),
            (V,  UP,    Option::None),
            (V,  DOWN,  Option::None),
            (H,  LEFT,  Option::None),
            (H,  RIGHT, Option::None),
            (V,  None,  Option::None),
            (H,  None,  Option::None),
        ];

        for (ntype, incoming, expected) in cases {
            assert_eq!(node(ntype, incoming).sign(), expected);
        }
    }

    #[test]
    fn test_traverse() {
        let cases = [
            (XL, [2, 3, 0, 1]),
            (XR, [2, 3, 0, 1]),
            (V,  [3, 2, 1, 0]),
            (H,  [1, 0, 3, 2]),
        ];
        for (ntype, expected) in cases {
            let c = node(ntype, None);
            for s in Slot::ALL {
                assert_eq!(c.paired_slot(s), Slot::from(expected[s.index()]));
            }
        }
    }
}
