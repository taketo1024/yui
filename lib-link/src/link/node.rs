use std::fmt::Display;

use yui_core::bitseq::Bit;
use yui_core::{CloneAnd, Sign};

use crate::Path;
use super::Edge;

use NodeType::{XL, XR, V, H};

// NodeType:
//
//     3   2         3   2         3   2         3   2        
//      \ /           \ /           \ /           \_/         
//       \    = XL,    /    = XR,   | |   = V,     _    = H, 
//      / \           / \           / \           / \         
//     0   1         0   1         0   1         0   1        
//

#[derive(Clone, Copy, PartialEq, Eq, Hash, derive_more::Display, Debug)]
pub enum NodeType { 
    XL, XR, V, H 
}

impl NodeType { 
    pub fn mirror(&self) -> Self {
        match self { 
            XL => XR,
            XR => XL,
            _  => *self
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, Default, derive_more::Display, Debug)]
pub enum NodeOri {
    #[default]
    #[display("-")] None,
    #[display("↑")] Up,
    #[display("↓")] Down,
    #[display("←")] Left,
    #[display("→")] Right,
}

impl NodeOri { 
    pub fn rev(&self) -> NodeOri {
        use NodeOri::*;
        match self {
            Up    => Down,
            Down  => Up,
            Left  => Right,
            Right => Left,
            None  => None,
        }
    }

    fn is_compatible(&self, ntype: NodeType) -> bool {
        use NodeType::*;
        use NodeOri::*;

        match (ntype, self) {
            (V, Right) | (V, Left) |
            (H, Up)    | (H, Down) => false,
            _ => true
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Node { 
    ntype: NodeType,
    pub(crate) ori: NodeOri,
    edges: [Edge; 4],
}

impl Node {
    pub fn new(ntype: NodeType, ori: NodeOri, edges: [Edge; 4]) -> Self { 
        assert!(ori.is_compatible(ntype), "Invalid (node-type, ori) combination: ({ntype}, {ori})");
        Node { ntype, edges, ori }
    }

    pub fn from_pd_code(edges: [Edge; 4]) -> Self { 
        Node::new(NodeType::XL, NodeOri::None, edges)
    }

    pub fn ntype(&self) -> NodeType { 
        self.ntype
    }

    pub fn edge(&self, i: usize) -> Edge { 
        assert!(i < 4);
        self.edges[i]
    }

    pub fn edges(&self) -> &[Edge; 4] { 
        &self.edges
    }

    pub fn min_edge(&self) -> Edge { 
        self.edges.iter().min().unwrap().clone()
    }

    pub fn is_crossing(&self) -> bool { 
        matches!(self.ntype, XL | XR)
    }

    pub fn is_resolved(&self) -> bool { 
        matches!(self.ntype, V | H)
    }

    pub fn resolve(&self, r: Bit) -> Self {
        use Bit::{Bit0, Bit1};

        self.clone_and(|x| {
            x.ntype = match (x.ntype, r) {
                (XL, Bit0) | (XR, Bit1) => H,
                (XL, Bit1) | (XR, Bit0) => V,
                _ => panic!("cannot resolve node-type: {}", x.ntype)
            };

            if !x.ori.is_compatible(x.ntype) {
                x.ori = NodeOri::None
            }
        })
    }

    pub fn is_oriented(&self) -> bool { 
        self.ori != NodeOri::None
    }

    pub fn ori(&self) -> NodeOri { 
        self.ori
    }

    pub fn is_pos(&self) -> bool { 
        self.sign().map(|x| x.is_positive()).unwrap_or(false)
    }

    pub fn is_neg(&self) -> bool { 
        self.sign().map(|x| x.is_negative()).unwrap_or(false)
    }

    pub fn sign(&self) -> Option<Sign> { 
        use NodeType::*;
        use NodeOri::*;

        match (self.ntype, self.ori) { 
            (XL, Left) | (XL, Right) | (XR, Up) | (XR, Down)  => Some(Sign::Pos),
            (XL, Up) | (XL, Down) |(XR, Left) | (XR, Right)   => Some(Sign::Neg),
            (_, None) | (V, Up) | (V, Down) | (H, Left) | (H, Right) => Option::None, 
            _ => panic!("Invalid (node-type, ori) combination: ({}, {})", self.ntype, self.ori)
        }
    }

    pub fn mirror(&self) -> Self { 
        self.clone_and(|x| 
            x.ntype = self.ntype.mirror()
        )
    }

    pub fn is_adj_to(&self, x: &Node) -> bool { 
        self.edges.iter().any(|e| x.edges.contains(e))
    }

    pub fn arcs(&self) -> (Path, Path) {
        let comp = |i: usize, j: usize| {
            let (ei, ej) = (self.edges[i], self.edges[j]);
            if ei == ej { 
                Path::new(vec![ei], true)
            } else { 
                Path::new(vec![ei, ej], false)
            }
        };
        match self.ntype { 
            XL | 
            XR => (comp(0, 2), comp(1, 3)),
            V  => (comp(0, 3), comp(1, 2)),
            H  => (comp(0, 1), comp(2, 3))
        }
    }

    pub fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge { 
        Self { 
            ntype: self.ntype, 
            ori:   self.ori,
            edges: self.edges.map(|e| f(e))
        }
    }

    pub(crate) fn traverse_inner(&self, index:usize) -> usize { 
        debug_assert!((0..4).contains(&index));

        match self.ntype {
            XL | XR => (index + 2) % 4,
            V => 3 - index,
            H => (5 - index) % 4
        }
    }
}

impl Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{:?}", self.ntype, self.edges)
    }
}

#[cfg(test)]
mod tests { 
    use super::*;
    
    fn node(ntype: NodeType, ori: NodeOri) -> Node {
        Node::new(ntype, ori, [0, 1, 2, 3])
    }

    #[test]
    fn test_is_resolved() {
        let c = node(XL, NodeOri::None);
        assert!(c.is_crossing());
        assert!(!c.is_resolved());

        let c = node(XR, NodeOri::None);
        assert!(c.is_crossing());
        assert!(!c.is_resolved());

        let c = node(H, NodeOri::None);
        assert!(!c.is_crossing());
        assert!(c.is_resolved());

        let c = node(V, NodeOri::None);
        assert!(!c.is_crossing());
        assert!(c.is_resolved());
    }

    #[test]
    fn test_resolve() {
        use Bit::{Bit0, Bit1};
        use NodeOri::*;

        // (ntype, ori, bit, expected_ntype, expected_ori)
        let cases = [
            (XL, None, Bit0, H, None),
            (XL, None, Bit1, V, None),
            (XR, None, Bit0, V, None),
            (XR, None, Bit1, H, None),
            // ori preserved when compatible with the resolved ntype.
            (XL, Up,   Bit1, V, Up),
            (XL, Left, Bit0, H, Left),
            // ori reset to None when incompatible with the resolved ntype.
            (XL, Up,   Bit0, H, None),
            (XL, Left, Bit1, V, None),
            (XR, Left, Bit0, V, None),
            (XR, Up,   Bit1, H, None),
        ];

        for (ntype, ori, bit, expected_ntype, expected_ori) in cases {
            let c = node(ntype, ori).resolve(bit);
            assert!(c.is_resolved());
            assert_eq!(c.ntype(), expected_ntype);
            assert_eq!(c.ori(), expected_ori);
        }
    }

    #[test]
    fn test_mirror() {
        use NodeOri::*;

        // (ntype, ori, expected_ntype, expected_ori) — mirror flips XL <-> XR and preserves ori.
        let cases = [
            (XL, None,  XR, None),
            (XR, None,  XL, None),
            (H,  None,  H,  None),
            (V,  None,  V,  None),
            (XL, Up,    XR, Up),
            (XR, Left,  XL, Left),
            (V,  Down,  V,  Down),
            (H,  Right, H,  Right),
        ];

        for (ntype, ori, expected_ntype, expected_ori) in cases {
            let c = node(ntype, ori).mirror();
            assert_eq!(c.ntype(), expected_ntype);
            assert_eq!(c.ori(), expected_ori);
        }
    }

    #[test]
    fn test_sign() {
        use NodeOri::*;

        // (ntype, ori, expected_sign)
        let cases = [
            (XL, Left,  Some(Sign::Pos)),
            (XL, Right, Some(Sign::Pos)),
            (XR, Up,    Some(Sign::Pos)),
            (XR, Down,  Some(Sign::Pos)),
            (XL, Up,    Some(Sign::Neg)),
            (XL, Down,  Some(Sign::Neg)),
            (XR, Left,  Some(Sign::Neg)),
            (XR, Right, Some(Sign::Neg)),
            // unsigned: no orientation, or resolved node.
            (XL, None,  Option::None),
            (XR, None,  Option::None),
            (V,  Up,    Option::None),
            (V,  Down,  Option::None),
            (H,  Left,  Option::None),
            (H,  Right, Option::None),
            (V,  None,  Option::None),
            (H,  None,  Option::None),
        ];

        for (ntype, ori, expected) in cases {
            assert_eq!(node(ntype, ori).sign(), expected);
        }
    }

    #[test]
    fn test_traverse() {
        let c = node(XL, NodeOri::None);
        assert_eq!(c.traverse_inner(0), 2);
        assert_eq!(c.traverse_inner(1), 3);
        assert_eq!(c.traverse_inner(2), 0);
        assert_eq!(c.traverse_inner(3), 1);

        let c = node(XR, NodeOri::None);
        assert_eq!(c.traverse_inner(0), 2);
        assert_eq!(c.traverse_inner(1), 3);
        assert_eq!(c.traverse_inner(2), 0);
        assert_eq!(c.traverse_inner(3), 1);

        let c = node(V, NodeOri::None);
        assert_eq!(c.traverse_inner(0), 3);
        assert_eq!(c.traverse_inner(1), 2);
        assert_eq!(c.traverse_inner(2), 1);
        assert_eq!(c.traverse_inner(3), 0);

        let c = node(H, NodeOri::None);
        assert_eq!(c.traverse_inner(0), 1);
        assert_eq!(c.traverse_inner(1), 0);
        assert_eq!(c.traverse_inner(2), 3);
        assert_eq!(c.traverse_inner(3), 2);
    }
}