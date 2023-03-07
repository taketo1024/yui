use std::fmt::Display;

use itertools::Itertools;

use crate::{Edge, Link};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Component { 
    edges:Vec<Edge>,
    closed:bool
}

impl Component { 
    pub fn new(edges: Vec<Edge>, closed: bool) -> Self { 
        assert!(!closed || edges.len() > 0);
        Self { edges, closed }
    }

    pub fn empty() -> Self { 
        Self::new(vec![], false)
    }

    pub fn edges(&self) -> &Vec<Edge> { 
        &self.edges
    }

    pub fn ends(&self) -> Option<(Edge, Edge)> { 
        if self.is_empty() || self.is_closed() { 
            None
        } else { 
            let &e0 = self.edges.first().unwrap();
            let &e1 = self.edges.last().unwrap();
            Some((e0, e1))
        }
    }

    pub fn len(&self) -> usize { 
        self.edges.len()
    }
    
    pub fn is_empty(&self) -> bool { 
        self.edges.is_empty()
    }

    pub fn is_closed(&self) -> bool { 
        self.closed
    }

    pub fn append(&mut self, e: Edge) { 
        assert_eq!(self.is_closed(), false);

        if self.edges.len() > 0 && self.edges[0] == e { 
            self.closed = true
        } else { 
            self.edges.push(e)
        }
    }

    pub fn reduce(&mut self) { 
        if self.is_empty() { 
            return
        }

        if self.is_closed() && self.len() > 1 { 
            let e0 = self.edges.remove(0);
            self.edges = vec![e0];
        } else if !self.is_closed() && self.len() > 2 { 
            let e0 = self.edges.remove(0);
            let e1 = self.edges.pop().unwrap();
            self.edges = vec![e0, e1];
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

        let Component {mut edges, ..} = other;

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

    pub fn is_adj(&self, other: &Component, link: &Link) -> bool { 
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

impl Display for Component {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[")?;
        write!(f, "{}", self.edges.iter().map(|e| e.to_string()).join("-"))?;
        if self.closed { 
            write!(f, "-")?;
        }
        write!(f, "]")
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn append() { 
        let mut c = Component::empty();
        c.append(0);
        c.append(1);
        c.append(2);

        assert_eq!(c.edges(), &vec![0,1,2]);
        assert_eq!(c.is_closed(), false);

        c.append(0);

        assert_eq!(c.edges(), &vec![0,1,2]);
        assert_eq!(c.is_closed(), true);
    }

    #[test]
    fn reduce() { 
        let mut c = Component::new(vec![], false);
        c.reduce();
        assert_eq!(c, Component::new(vec![], false));
        
        let mut c = Component::new(vec![0], false);
        c.reduce();
        assert_eq!(c, Component::new(vec![0], false));
        
        let mut c = Component::new(vec![0,1,2,3], false);
        c.reduce();
        assert_eq!(c, Component::new(vec![0, 3], false));
        
        let mut c = Component::new(vec![0], true);
        c.reduce();
        assert_eq!(c, Component::new(vec![0], true));

        let mut c = Component::new(vec![0,1,2,3], true);
        c.reduce();
        assert_eq!(c, Component::new(vec![0], true));
    }

    #[test]
    fn is_connectable() { 
        let c = Component::new(vec![1,2,3,4], false);

        assert_eq!(c.is_connectable(&Component::new(vec![4,5], false)), true);
        assert_eq!(c.is_connectable(&Component::new(vec![5,4], false)), true);
        assert_eq!(c.is_connectable(&Component::new(vec![1,5], false)), true);
        assert_eq!(c.is_connectable(&Component::new(vec![5,1], false)), true);
        assert_eq!(c.is_connectable(&Component::new(vec![5,6], false)), false);
        assert_eq!(c.is_connectable(&Component::new(vec![2,3], false)), false);
        assert_eq!(c.is_connectable(&Component::new(vec![0], true)), false);
        assert_eq!(c.is_connectable(&Component::new(vec![], false)), false);
    }

    #[test]
    fn connect() { 
        let mut c = Component::new(vec![1,2,3,4], false);
        c.connect(Component::new(vec![4,5], false));
        assert_eq!(c, Component::new(vec![1,2,3,4,5], false));

        let mut c = Component::new(vec![1,2,3,4], false);
        c.connect(Component::new(vec![5,4], false));
        assert_eq!(c, Component::new(vec![1,2,3,4,5], false));

        let mut c = Component::new(vec![1,2,3,4], false);
        c.connect(Component::new(vec![6,1], false));
        assert_eq!(c, Component::new(vec![6,1,2,3,4], false));

        let mut c = Component::new(vec![1,2,3,4], false);
        c.connect(Component::new(vec![1,6], false));
        assert_eq!(c, Component::new(vec![6,1,2,3,4], false));

        let mut c = Component::new(vec![1,2,3,4], false);
        c.connect(Component::new(vec![1], false));
        assert_eq!(c, Component::new(vec![1,2,3,4], false));

        let mut c = Component::new(vec![1,2,3,4], false);
        c.connect(Component::new(vec![4], false));
        assert_eq!(c, Component::new(vec![1,2,3,4], false));
    }
}