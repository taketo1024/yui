use std::fmt::Display;
use derive_more::Display;
use crate::LinkComp;

use super::{Edge, Resolution};

use Resolution::{Res0, Res1};
use CrossingType::{Xp, Xn, V, H};

#[derive(Debug, Clone, Copy, PartialEq, Display)]
pub enum CrossingType { 
    Xp, Xn, V, H 
}

impl CrossingType { 
    pub fn mirror(self) -> CrossingType {
        match self { 
            Xp => Xn,
            Xn => Xp,
            other => other
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Crossing { 
    ctype: CrossingType,
    edges: [Edge; 4]
}

impl Crossing {
    pub fn new(ctype: CrossingType, edges: [Edge; 4]) -> Self { 
        Crossing { ctype, edges }
    }

    pub fn from_pd_code(edges: [Edge; 4]) -> Self { 
        Crossing::new(CrossingType::Xn, edges)
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
        match self.ctype { 
            V | H => true,
            _ => false
        }
    }

    pub fn resolve(&mut self, r: Resolution) {
        match (self.ctype, r) {
            (Xp, Res0) | (Xn, Res1) => self.ctype = V,
            (Xp, Res1) | (Xn, Res0) => self.ctype = H,
            _ => panic!()
        }
    }

    pub fn mirror(&mut self) { 
        self.ctype = self.ctype.mirror();
    }

    pub fn pass(&self, index:usize) -> usize { 
        debug_assert!((0..4).contains(&index));

        match self.ctype {
            Xp | Xn => (index + 2) % 4,
            V => 3 - index,
            H => (5 - index) % 4
        }
    }

    pub fn arcs(&self) -> (LinkComp, LinkComp) {
        let comp = |i: usize, j: usize| {
            LinkComp::new(vec![self.edges[i], self.edges[j]], false)
        };
        match self.ctype { 
            Xp | 
            Xn => (comp(0, 2), comp(1, 3)),
            V  => (comp(0, 3), comp(1, 2)),
            H  => (comp(0, 1), comp(2, 3))
        }
    }

    pub fn res_arcs(&self, r: Resolution) -> (LinkComp, LinkComp) { 
        let mut x = self.clone();
        x.resolve(r);
        x.arcs()        
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
        let c = a_crossing(Xp);
        assert!(!c.is_resolved());

        let c = a_crossing(Xn);
        assert!(!c.is_resolved());

        let c = a_crossing(H);
        assert!(c.is_resolved());

        let c = a_crossing(V);
        assert!(c.is_resolved());
    }

    #[test]
    fn crossing_resolve() {
        let mut c = a_crossing(Xp);
        
        c.resolve(Res0);
        assert!(c.is_resolved());
        assert_eq!(c.ctype(), V);

        let mut c = a_crossing(Xp);
        
        c.resolve(Res1);
        assert!(c.is_resolved());
        assert_eq!(c.ctype(), H);

        let mut c = a_crossing(Xn);
        
        c.resolve(Res0);
        assert!(c.is_resolved());
        assert_eq!(c.ctype(), H);

        let mut c = a_crossing(Xn);
        
        c.resolve(Res1);
        assert!(c.is_resolved());
        assert_eq!(c.ctype(), V);
    }

    #[test]
    fn crossing_mirror() {
        let mut c = a_crossing(Xp);
        c.mirror();
        assert_eq!(c.ctype(), Xn);

        let mut c = a_crossing(Xn);
        c.mirror();
        assert_eq!(c.ctype(), Xp);

        let mut c = a_crossing(H);
        c.mirror();
        assert_eq!(c.ctype(), H);

        let mut c = a_crossing(V);
        c.mirror();
        assert_eq!(c.ctype(), V);
    }

    #[test]
    fn crossing_pass() {
        let c = a_crossing(Xp);
        assert_eq!(c.pass(0), 2);
        assert_eq!(c.pass(1), 3);
        assert_eq!(c.pass(2), 0);
        assert_eq!(c.pass(3), 1);

        let c = a_crossing(Xn);
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