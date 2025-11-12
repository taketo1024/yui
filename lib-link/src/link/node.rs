use std::fmt::Display;

use yui::bitseq::Bit;
use yui::CloneAnd;

use crate::Path;
use super::Edge;

use NodeType::{X, Xm, V, H};

#[derive(Clone, Copy, PartialEq, Eq, Hash, derive_more::Display, Debug)]
pub enum NodeType { 
    X, Xm, V, H 
}

impl NodeType { 
    pub fn cc(&mut self) {
        match self { 
            Xm => *self = X,
            X  => *self = Xm,
            _ => ()
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
        Node::new(NodeType::X, edges)
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
        matches!(self.ntype, X | Xm)
    }

    pub fn is_resolved(&self) -> bool { 
        matches!(self.ntype, V | H)
    }

    pub fn resolve(&mut self, r: Bit) {
        use Bit::{Bit0, Bit1};

        match (self.ntype, r) {
            (X, Bit0) | (Xm, Bit1) => self.ntype = H,
            (X, Bit1) | (Xm, Bit0) => self.ntype = V,
            _ => panic!()
        }
    }

    pub fn resolved(&self, r: Bit) -> Self { 
        self.clone_and(|x|
            x.resolve(r)
        )
    }

    pub fn cc(&mut self) { 
        self.ntype.cc()
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
            X | 
            Xm => (comp(0, 2), comp(1, 3)),
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
            X | Xm => (index + 2) % 4,
            V => 3 - index,
            H => (5 - index) % 4
        }
    }
}

impl From<[Edge; 4]> for Node {
    fn from(edges: [Edge; 4]) -> Self {
        Self::new(NodeType::X, edges)
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
        let c = a_crossing(X);
        assert!(c.is_crossing());

        let c = a_crossing(Xm);
        assert!(c.is_crossing());

        let c = a_crossing(H);
        assert!(c.is_resolved());

        let c = a_crossing(V);
        assert!(c.is_resolved());
    }

    #[test]
    fn crossing_resolve() {
        use Bit::{Bit0, Bit1};

        let mut c = a_crossing(X);
        
        c.resolve(Bit0);
        assert!(c.is_resolved());
        assert_eq!(c.ntype(), H);

        let mut c = a_crossing(X);
        
        c.resolve(Bit1);
        assert!(c.is_resolved());
        assert_eq!(c.ntype(), V);

        let mut c = a_crossing(Xm);
        
        c.resolve(Bit0);
        assert!(c.is_resolved());
        assert_eq!(c.ntype(), V);

        let mut c = a_crossing(Xm);
        
        c.resolve(Bit1);
        assert!(c.is_resolved());
        assert_eq!(c.ntype(), H);
    }

    #[test]
    fn crossing_mirror() {
        let c = a_crossing(X).clone_and(|c| c.cc());
        assert_eq!(c.ntype(), Xm);

        let c = a_crossing(Xm).clone_and(|c| c.cc());
        assert_eq!(c.ntype(), X);

        let c = a_crossing(H).clone_and(|c| c.cc());
        assert_eq!(c.ntype(), H);

        let c = a_crossing(V).clone_and(|c| c.cc());
        assert_eq!(c.ntype(), V);
    }

    #[test]
    fn crossing_pass() {
        let c = a_crossing(X);
        assert_eq!(c.traverse_inner(0), 2);
        assert_eq!(c.traverse_inner(1), 3);
        assert_eq!(c.traverse_inner(2), 0);
        assert_eq!(c.traverse_inner(3), 1);

        let c = a_crossing(Xm);
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