use std::fmt::Display;

use yui_core::bitseq::Bit;
use yui_core::{CloneAnd, Sign};

use crate::Path;
use super::Edge;

use NodeType::{XL, XR, V, H};

#[derive(Clone, Copy, PartialEq, Eq, Hash, derive_more::Display, Debug)]
pub enum NodeType { 
    XL, XR, V, H 
}

impl NodeType { 
    pub fn mirror(&self) -> Self {
        match self { 
            XR => XL,
            XL => XR,
            _  => *self
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Node { 
    ntype: NodeType,
    edges: [Edge; 4]
}

impl Node {
    pub fn new(ntype: NodeType, edges: [Edge; 4]) -> Self { 
        Node { ntype, edges }
    }

    pub fn from_pd_code(edges: [Edge; 4]) -> Self { 
        Node::new(NodeType::XL, edges)
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

        self.clone_and(|x| 
            x.ntype = match (x.ntype, r) {
                (XL, Bit0) | (XR, Bit1) => H,
                (XL, Bit1) | (XR, Bit0) => V,
                _ => panic!()
            }
        )
    }

    pub fn sign(&self, j: usize) -> Option<Sign> {
        match (self.ntype, j) { 
            (XR, 1) | (XL, 3) => Some(Sign::Pos),
            (XR, 3) | (XL, 1) => Some(Sign::Neg),
            _ => None
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

impl From<[Edge; 4]> for Node {
    fn from(edges: [Edge; 4]) -> Self {
        Self::new(NodeType::XL, edges)
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
    
    fn a_crossing(ntype:NodeType) -> Node {
        Node{
            ntype, 
            edges: [0,1,2,3]
        }
    }

    #[test]
    fn crossing_is_resolved() {
        let c = a_crossing(XL);
        assert!(c.is_crossing());

        let c = a_crossing(XR);
        assert!(c.is_crossing());

        let c = a_crossing(H);
        assert!(c.is_resolved());

        let c = a_crossing(V);
        assert!(c.is_resolved());
    }

    #[test]
    fn crossing_resolve() {
        use Bit::{Bit0, Bit1};

        let c = a_crossing(XL).resolve(Bit0);
        assert!(c.is_resolved());
        assert_eq!(c.ntype(), H);

        let c = a_crossing(XL).resolve(Bit1);
        assert!(c.is_resolved());
        assert_eq!(c.ntype(), V);

        let c = a_crossing(XR).resolve(Bit0);
        assert!(c.is_resolved());
        assert_eq!(c.ntype(), V);

        let c = a_crossing(XR).resolve(Bit1);
        assert!(c.is_resolved());
        assert_eq!(c.ntype(), H);
    }

    #[test]
    fn crossing_mirror() {
        let c = a_crossing(XL).mirror();
        assert_eq!(c.ntype(), XR);

        let c = a_crossing(XR).mirror();
        assert_eq!(c.ntype(), XL);

        let c = a_crossing(H).mirror();
        assert_eq!(c.ntype(), H);

        let c = a_crossing(V).mirror();
        assert_eq!(c.ntype(), V);
    }

    #[test]
    fn crossing_pass() {
        let c = a_crossing(XL);
        assert_eq!(c.traverse_inner(0), 2);
        assert_eq!(c.traverse_inner(1), 3);
        assert_eq!(c.traverse_inner(2), 0);
        assert_eq!(c.traverse_inner(3), 1);

        let c = a_crossing(XR);
        assert_eq!(c.traverse_inner(0), 2);
        assert_eq!(c.traverse_inner(1), 3);
        assert_eq!(c.traverse_inner(2), 0);
        assert_eq!(c.traverse_inner(3), 1);

        let c = a_crossing(V);
        assert_eq!(c.traverse_inner(0), 3);
        assert_eq!(c.traverse_inner(1), 2);
        assert_eq!(c.traverse_inner(2), 1);
        assert_eq!(c.traverse_inner(3), 0);

        let c = a_crossing(H);
        assert_eq!(c.traverse_inner(0), 1);
        assert_eq!(c.traverse_inner(1), 0);
        assert_eq!(c.traverse_inner(2), 3);
        assert_eq!(c.traverse_inner(3), 2);
    }
}