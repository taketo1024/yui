use std::{fmt::Display, cmp::min};

use itertools::Itertools;

use crate::{Edge, Link};

#[derive(Debug, Clone, Eq)]
pub struct LinkComp { 
    edges:Vec<Edge>,
    closed:bool
}

impl LinkComp { 
    pub fn new(edges: Vec<Edge>, closed: bool) -> Self { 
        assert!(!edges.is_empty());
        Self { edges, closed }
    }

    pub fn arc(edges: Vec<Edge>) -> Self { 
        Self::new(edges, false)
    }

    pub fn circ(edges: Vec<Edge>) -> Self { 
        Self::new(edges, true)
    }

    pub fn contains(&self, e: Edge) -> bool { 
        self.edges.contains(&e)
    }

    pub fn edges(&self) -> &Vec<Edge> { 
        &self.edges
    }

    pub fn min_edge(&self) -> Edge { 
        self.edges.iter().min().unwrap().clone()
    }

    pub fn ends(&self) -> Option<(Edge, Edge)> { 
        if self.is_arc() { 
            let &e0 = self.edges.first().unwrap();
            let &e1 = self.edges.last().unwrap();
            Some((e0, e1))
        } else { 
            None
        }
    }

    pub fn len(&self) -> usize { 
        self.edges.len()
    }
    
    pub fn is_arc(&self) -> bool { 
        !self.closed
    }

    pub fn is_circle(&self) -> bool { 
        self.closed
    }

    pub fn reduce(&mut self) { 
        if self.is_arc() && self.len() > 2 { 
            let e0 = self.edges.remove(0);
            let e1 = self.edges.pop().unwrap();
            let min = self.edges().iter().filter(|&e| e < min(&e0, &e1)).min();

            self.edges = if let Some(&e2) = min {
                vec![e0, e2, e1]
            } else { 
                vec![e0, e1]
            };
        } else if self.is_circle() && self.len() > 1 { 
            let e0 = self.edges.iter().min().unwrap().clone();
            self.edges = vec![e0];
        }
    }

    pub fn is_connectable(&self, other: &Self) -> bool { 
        let Some((e0, e1)) =  self.ends() else { return false };
        let Some((f0, f1)) = other.ends() else { return false };
        
        e0 == f0 || e0 == f1 || e1 == f0 || e1 == f1
    }

    pub fn connect(&mut self, other: Self) { 
        assert!(self.is_connectable(&other), "{self} and {other} is not connectable.");

        let (e0, e1) =  self.ends().unwrap();
        let (f0, f1) = other.ends().unwrap();

        let LinkComp {mut edges, ..} = other;

        if e1 == f0 {        // [.., e1) + [f0, ..)
            edges.remove(0);
            self.edges.append(&mut edges);
        } else if e1 == f1 { // [.., e1) + [.., f1)
            edges.pop();
            edges.reverse();
            self.edges.append(&mut edges);
        } else if e0 == f0 { // [e0, ..) + [f0, ..)
            edges.remove(0);
            edges.reverse();
            edges.append(&mut self.edges);
            self.edges = edges;
        } else if e0 == f1 { // [e0, ..) + [.., f1)
            edges.pop();
            edges.append(&mut self.edges);
            self.edges = edges;
        } else {
            panic!() // unreachable
        }

        let (e0, e1) = self.ends().unwrap();

        if e0 == e1 { 
            self.edges.pop();
            self.closed = true;
        }
    }

    pub fn is_adj(&self, other: &LinkComp, link: &Link) -> bool { 
        // 1) find crossings `x` that touche `self`. 
        // 2) check if `x` also touches `other`.
        
        for x in link.data().iter() { 
            if !x.edges().iter().any(|e| 
                self.edges.contains(e)
            ) { 
                continue
            }

            let Some(e) = x.edges().iter().find(|e| 
                !self.edges.contains(e)
            ) else { 
                continue
            };

            if other.edges.contains(e) { 
                return true
            }
        }

        false
    }
}

impl Display for LinkComp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let c = self.edges.iter().map(|e| e.to_string()).join("-");
        if self.is_circle() { 
            write!(f, "[-{c}-]")
        } else {
            write!(f, "[{c}]")
        }
    }
}

impl PartialEq for LinkComp {
    fn eq(&self, other: &Self) -> bool {
        if self.closed != other.closed { 
            return false;
        }
        if self.edges.len() != other.edges.len() { 
            return false;
        }
        if self.edges == other.edges { 
            return true
        } 
        
        if self.closed { 
            let n = self.edges.len();
            let Some(p) = other.edges.iter().position(|e| e == &self.edges[0]) else { 
                return false
            };
            (0 .. n).all(|i| &self.edges[i] == &other.edges[(p + i) % n]) || 
            (0 .. n).all(|i| &self.edges[i] == &other.edges[(p + n - i) % n])
        } else { 
            self.edges.iter().zip(other.edges.iter().rev()).all(|(e, f)| e == f)
        }
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn reduce() { 
        let mut c = LinkComp::new(vec![0], false);
        c.reduce();
        assert_eq!(c, LinkComp::new(vec![0], false));
        
        let mut c = LinkComp::new(vec![0,1,2,3], false);
        c.reduce();
        assert_eq!(c, LinkComp::new(vec![0, 3], false));
        
        let mut c = LinkComp::new(vec![0], true);
        c.reduce();
        assert_eq!(c, LinkComp::new(vec![0], true));

        let mut c = LinkComp::new(vec![0,1,2,3], true);
        c.reduce();
        assert_eq!(c, LinkComp::new(vec![0], true));
    }

    #[test]
    fn is_connectable() { 
        let c = LinkComp::new(vec![1,2,3,4], false);

        assert_eq!(c.is_connectable(&LinkComp::new(vec![4,5], false)), true);
        assert_eq!(c.is_connectable(&LinkComp::new(vec![5,4], false)), true);
        assert_eq!(c.is_connectable(&LinkComp::new(vec![1,5], false)), true);
        assert_eq!(c.is_connectable(&LinkComp::new(vec![5,1], false)), true);
        assert_eq!(c.is_connectable(&LinkComp::new(vec![5,6], false)), false);
        assert_eq!(c.is_connectable(&LinkComp::new(vec![2,3], false)), false);
        assert_eq!(c.is_connectable(&LinkComp::new(vec![0], true)), false);
    }

    #[test]
    fn connect() { 
        let mut c = LinkComp::new(vec![1,2,3,4], false);
        c.connect(LinkComp::new(vec![4,5], false));
        assert_eq!(c, LinkComp::new(vec![1,2,3,4,5], false));

        let mut c = LinkComp::new(vec![1,2,3,4], false);
        c.connect(LinkComp::new(vec![5,4], false));
        assert_eq!(c, LinkComp::new(vec![1,2,3,4,5], false));

        let mut c = LinkComp::new(vec![1,2,3,4], false);
        c.connect(LinkComp::new(vec![6,1], false));
        assert_eq!(c, LinkComp::new(vec![6,1,2,3,4], false));

        let mut c = LinkComp::new(vec![1,2,3,4], false);
        c.connect(LinkComp::new(vec![1,6], false));
        assert_eq!(c, LinkComp::new(vec![6,1,2,3,4], false));

        let mut c = LinkComp::new(vec![1,2,3,4], false);
        c.connect(LinkComp::new(vec![1], false));
        assert_eq!(c, LinkComp::new(vec![1,2,3,4], false));

        let mut c = LinkComp::new(vec![1,2,3,4], false);
        c.connect(LinkComp::new(vec![4], false));
        assert_eq!(c, LinkComp::new(vec![1,2,3,4], false));
    }

    #[test]
    fn eq_arc() { 
        let c = LinkComp::arc(vec![1,2,3]);

        assert_eq!(c, LinkComp::arc(vec![1,2,3]));
        assert_eq!(c, LinkComp::arc(vec![3,2,1]));

        assert!(c != LinkComp::arc(vec![1,2]));
        assert!(c != LinkComp::arc(vec![1,2,3,4]));
        assert!(c != LinkComp::circ(vec![1,2,3]));
    }

    #[test]
    fn eq_circ() { 
        let c = LinkComp::circ(vec![1,2,3,4]);

        assert_eq!(c, LinkComp::circ(vec![1,2,3,4]));
        assert_eq!(c, LinkComp::circ(vec![2,3,4,1]));
        assert_eq!(c, LinkComp::circ(vec![3,4,1,2]));
        assert_eq!(c, LinkComp::circ(vec![2,3,4,1]));
        assert_eq!(c, LinkComp::circ(vec![4,3,2,1]));
        assert_eq!(c, LinkComp::circ(vec![3,2,1,4]));
        
        assert!(c != LinkComp::circ(vec![1,2,3]));
        assert!(c != LinkComp::circ(vec![1,2,3,4,5]));
        assert!(c != LinkComp::circ(vec![1,2,4,3]));
    }
}