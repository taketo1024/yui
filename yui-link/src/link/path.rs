use std::{fmt::Display, cmp::min};

use itertools::Itertools;

use crate::{Edge, Link};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Path { 
    edges:Vec<Edge>,
    closed:bool
}

impl Path { 
    pub fn new<I>(edges: I, closed: bool) -> Self
    where I: IntoIterator<Item = Edge> { 
        let edges = edges.into_iter().collect_vec();
        assert!(!edges.is_empty());
        Self { edges, closed }
    }

    pub fn arc<I>(edges: I) -> Self
    where I: IntoIterator<Item = Edge> { 
        Self::new(edges, false)
    }

    pub fn circ<I>(edges: I) -> Self
    where I: IntoIterator<Item = Edge> { 
        Self::new(edges, true)
    }

    pub fn contains(&self, e: Edge) -> bool { 
        self.edges.contains(&e)
    }

    pub fn len(&self) -> usize { 
        self.edges.len()
    }
    
    pub fn edges(&self) -> &Vec<Edge> { 
        &self.edges
    }

    pub fn min_edge(&self) -> Edge { 
        *self.edges.iter().min().unwrap()
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
            let e0 = *self.edges.iter().min().unwrap();
            self.edges = vec![e0];
        }
    }

    pub fn is_connectable(&self, other: &Self) -> bool { 
        let Some((e0, e1)) =  self.ends() else { return false };
        let Some((f0, f1)) = other.ends() else { return false };
        
        e0 == f0 || e0 == f1 || e1 == f0 || e1 == f1
    }

    pub fn is_connectable_bothends(&self, other: &Self) -> bool { 
        let Some((e0, e1)) =  self.ends() else { return false };
        let Some((f0, f1)) = other.ends() else { return false };

        (e0, e1) == (f0, f1) || (e0, e1) == (f1, f0)
    }

    pub fn connect(&mut self, other: Self) { 
        assert!(self.is_connectable(&other), "{self} and {other} is not connectable.");

        let (e0, e1) =  self.ends().unwrap();
        let (f0, f1) = other.ends().unwrap();

        let Path {mut edges, ..} = other.clone();

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

    pub fn is_adj(&self, other: &Path, link: &Link) -> bool { 
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

    fn edge_sum(&self) -> usize { 
        self.edges.iter().sum()
    }

    pub fn unori_eq(&self, other: &Self) -> bool {
        // early failure
        if self.closed != other.closed ||
           self.edges.len() != other.edges.len() || 
           self.edge_sum() != other.edge_sum()
        { 
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
            (0 .. n).all(|i| self.edges[i] == other.edges[(p + i) % n]) || 
            (0 .. n).all(|i| self.edges[i] == other.edges[(p + n - i) % n])
        } else { 
            Iterator::zip(
                self.edges.iter(), 
                other.edges.iter().rev()
            ).all(|(e, f)| 
                e == f
            )
        }
    }
}

impl Display for Path {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let c = self.edges.iter().map(|e| e.to_string()).join("-");
        if self.is_circle() { 
            write!(f, "⚪︎({c})")
        } else {
            write!(f, "[{c}]")
        }
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn reduce() { 
        let mut c = Path::new(vec![0], false);
        c.reduce();
        assert_eq!(c, Path::new(vec![0], false));
        
        let mut c = Path::new(vec![0,1,2,3], false);
        c.reduce();
        assert_eq!(c, Path::new(vec![0, 3], false));
        
        let mut c = Path::new(vec![0], true);
        c.reduce();
        assert_eq!(c, Path::new(vec![0], true));

        let mut c = Path::new(vec![0,1,2,3], true);
        c.reduce();
        assert_eq!(c, Path::new(vec![0], true));
    }

    #[test]
    fn is_connectable() { 
        let c = Path::new(vec![1,2,3,4], false);

        assert!( c.is_connectable(&Path::new(vec![4,5], false)));
        assert!( c.is_connectable(&Path::new(vec![5,4], false)));
        assert!( c.is_connectable(&Path::new(vec![1,5], false)));
        assert!( c.is_connectable(&Path::new(vec![5,1], false)));
        assert!(!c.is_connectable(&Path::new(vec![5,6], false)));
        assert!(!c.is_connectable(&Path::new(vec![2,3], false)));
        assert!(!c.is_connectable(&Path::new(vec![0],   true)));
    }

    #[test]
    fn is_connectable_bothends() { 
        let c = Path::new(vec![1,2,3,4], false);

        assert!(!c.is_connectable_bothends(&Path::new(vec![4,5], false)));
        assert!(!c.is_connectable_bothends(&Path::new(vec![5,4], false)));
        assert!(!c.is_connectable_bothends(&Path::new(vec![1,5], false)));
        assert!(!c.is_connectable_bothends(&Path::new(vec![5,1], false)));
        assert!( c.is_connectable_bothends(&Path::new(vec![1,4], false)));
        assert!( c.is_connectable_bothends(&Path::new(vec![4,1], false)));
        assert!(!c.is_connectable_bothends(&Path::new(vec![0],   true)));
    }

    #[test]
    fn connect() { 
        let mut c = Path::new(vec![1,2,3,4], false);
        c.connect(Path::new(vec![4,5], false));
        assert_eq!(c, Path::new(vec![1,2,3,4,5], false));

        let mut c = Path::new(vec![1,2,3,4], false);
        c.connect(Path::new(vec![5,4], false));
        assert_eq!(c, Path::new(vec![1,2,3,4,5], false));

        let mut c = Path::new(vec![1,2,3,4], false);
        c.connect(Path::new(vec![6,1], false));
        assert_eq!(c, Path::new(vec![6,1,2,3,4], false));

        let mut c = Path::new(vec![1,2,3,4], false);
        c.connect(Path::new(vec![1,6], false));
        assert_eq!(c, Path::new(vec![6,1,2,3,4], false));

        let mut c = Path::new(vec![1,2,3,4], false);
        c.connect(Path::new(vec![1], false));
        assert_eq!(c, Path::new(vec![1,2,3,4], false));

        let mut c = Path::new(vec![1,2,3,4], false);
        c.connect(Path::new(vec![4], false));
        assert_eq!(c, Path::new(vec![1,2,3,4], false));
    }

    #[test]
    fn unori_eq_arc() { 
        let c = Path::arc(vec![1,2,3]);

        assert!(c.unori_eq(&Path::arc(vec![1,2,3])));
        assert!(c.unori_eq(&Path::arc(vec![3,2,1])));

        assert!(!c.unori_eq(&Path::arc(vec![1,2])));
        assert!(!c.unori_eq(&Path::arc(vec![1,2,3,4])));
        assert!(!c.unori_eq(&Path::circ(vec![1,2,3])));
    }

    #[test]
    fn unori_eq_circ() { 
        let c = Path::circ(vec![1,2,3,4]);

        assert!(c.unori_eq(&Path::circ(vec![1,2,3,4])));
        assert!(c.unori_eq(&Path::circ(vec![2,3,4,1])));
        assert!(c.unori_eq(&Path::circ(vec![3,4,1,2])));
        assert!(c.unori_eq(&Path::circ(vec![2,3,4,1])));
        assert!(c.unori_eq(&Path::circ(vec![4,3,2,1])));
        assert!(c.unori_eq(&Path::circ(vec![3,2,1,4])));
        
        assert!(!c.unori_eq(&Path::circ(vec![1,2,3])));
        assert!(!c.unori_eq(&Path::circ(vec![1,2,3,4,5])));
        assert!(!c.unori_eq(&Path::circ(vec![1,2,4,3])));
    }
}