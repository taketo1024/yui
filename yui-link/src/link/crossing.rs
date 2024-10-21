use derive_more::Display;
use yui::bitseq::Bit;

use crate::Path;
use super::Edge;

use CrossingType::{X, Xm, V, H};

#[derive(Display, Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CrossingType { 
    X, Xm, V, H 
}

impl CrossingType { 
    pub fn mirror(self) -> CrossingType {
        match self { 
            Xm => X,
            X => Xm,
            other => other
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Crossing { 
    ctype: CrossingType,
    edges: [Edge; 4]
}

impl Crossing {
    pub fn new(ctype: CrossingType, edges: [Edge; 4]) -> Self { 
        Crossing { ctype, edges }
    }

    pub fn from_pd_code(edges: [Edge; 4]) -> Self { 
        Crossing::new(CrossingType::X, edges)
    }

    pub fn ctype(&self) -> CrossingType { 
        self.ctype
    }

    pub fn edge(&self, i: usize) -> Edge { 
        assert!(i < 4);
        self.edges[i]
    }

    pub fn edges(&self) -> &[Edge; 4] { 
        &self.edges
    }

    pub fn is_resolved(&self) -> bool { 
        matches!(self.ctype, V | H)
    }

    pub fn resolve(&mut self, r: Bit) {
        use Bit::{Bit0, Bit1};

        match (self.ctype, r) {
            (X, Bit0) | (Xm, Bit1) => self.ctype = H,
            (X, Bit1) | (Xm, Bit0) => self.ctype = V,
            _ => panic!()
        }
    }

    pub fn resolved(&self, r: Bit) -> Self { 
        let mut x = self.clone();
        x.resolve(r);
        x
    }

    pub fn mirror(&self) -> Self { 
        Self { 
            ctype: self.ctype.mirror(),
            edges: self.edges.clone()
        }
    }

    pub fn is_adj_to(&self, x: &Crossing) -> bool { 
        self.edges.iter().any(|e| x.edges.contains(e))
    }

    pub fn pass(&self, index:usize) -> usize { 
        debug_assert!((0..4).contains(&index));

        match self.ctype {
            X | Xm => (index + 2) % 4,
            V => 3 - index,
            H => (5 - index) % 4
        }
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
        match self.ctype { 
            X | 
            Xm => (comp(0, 2), comp(1, 3)),
            V  => (comp(0, 3), comp(1, 2)),
            H  => (comp(0, 1), comp(2, 3))
        }
    }

    pub fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge { 
        Self { 
            ctype: self.ctype, 
            edges: self.edges.map(|e| f(e)) 
        }
    }
}

impl Display for Crossing {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{:?}", self.ctype, self.edges)
    }
}

#[cfg(test)]
mod tests { 
    use super::*;
    
    fn a_crossing(ctype:CrossingType) -> Crossing {
        Crossing{
            ctype, 
            edges: [0,1,2,3]
        }
    }

    #[test]
    fn crossing_is_resolved() {
        let c = a_crossing(X);
        assert!(!c.is_resolved());

        let c = a_crossing(Xm);
        assert!(!c.is_resolved());

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
        assert_eq!(c.ctype(), H);

        let mut c = a_crossing(X);
        
        c.resolve(Bit1);
        assert!(c.is_resolved());
        assert_eq!(c.ctype(), V);

        let mut c = a_crossing(Xm);
        
        c.resolve(Bit0);
        assert!(c.is_resolved());
        assert_eq!(c.ctype(), V);

        let mut c = a_crossing(Xm);
        
        c.resolve(Bit1);
        assert!(c.is_resolved());
        assert_eq!(c.ctype(), H);
    }

    #[test]
    fn crossing_mirror() {
        let c = a_crossing(X).mirror();
        assert_eq!(c.ctype(), Xm);

        let c = a_crossing(Xm).mirror();
        assert_eq!(c.ctype(), X);

        let c = a_crossing(H).mirror();
        assert_eq!(c.ctype(), H);

        let c = a_crossing(V).mirror();
        assert_eq!(c.ctype(), V);
    }

    #[test]
    fn crossing_pass() {
        let c = a_crossing(X);
        assert_eq!(c.pass(0), 2);
        assert_eq!(c.pass(1), 3);
        assert_eq!(c.pass(2), 0);
        assert_eq!(c.pass(3), 1);

        let c = a_crossing(Xm);
        assert_eq!(c.pass(0), 2);
        assert_eq!(c.pass(1), 3);
        assert_eq!(c.pass(2), 0);
        assert_eq!(c.pass(3), 1);

        let c = a_crossing(V);
        assert_eq!(c.pass(0), 3);
        assert_eq!(c.pass(1), 2);
        assert_eq!(c.pass(2), 1);
        assert_eq!(c.pass(3), 0);

        let c = a_crossing(H);
        assert_eq!(c.pass(0), 1);
        assert_eq!(c.pass(1), 0);
        assert_eq!(c.pass(2), 3);
        assert_eq!(c.pass(3), 2);
    }
}