use std::collections::HashSet;
use std::fmt::Display;
use itertools::Itertools;
use yui_link::{Edge, Crossing, Resolution, LinkComp};

#[derive(Clone, Copy, Eq, PartialOrd, Ord, Debug)]
pub enum TngComp { 
    Arc(Edge, Edge, Edge), // the middle edge keeps the min edge-id.
    Circ(Edge)
}

impl TngComp { 
    pub fn arc(e0: Edge, e1: Edge) -> Self { 
        assert!(e0 != e1);
        if e0 < e1 { 
            Self::Arc(e0, e0, e1)
        } else { 
            Self::Arc(e1, e1, e0)
        }
    }

    pub fn circ(e: Edge) -> Self { 
        Self::Circ(e)
    }

    pub fn is_arc(&self) -> bool { 
        matches!(self, TngComp::Arc(..))
    }

    pub fn is_circle(&self) -> bool { 
        matches!(self, TngComp::Circ(..))
    }

    pub fn endpts(&self) -> Option<(Edge, Edge)> { 
        match self { 
            &TngComp::Arc(e0, _, e2) => Some((e0, e2)),
            _ => None
        }
    }

    pub fn is_connectable(&self, other: &Self) -> bool { 
        use TngComp::Arc;
        let Arc(l0, _, l2) = self  else { return false };
        let Arc(r0, _, r2) = other else { return false };
        l0 == r0 || l0 == r2 || l2 == r0 || l2 == r2
    }

    pub fn connect(&self, other: &Self) -> Self { 
        use TngComp::Arc;
        let Arc(l0, l1, l2) = self  else { panic!() };
        let Arc(r0, r1, r2) = other else { panic!() };

        let e1 = std::cmp::min(l1, r1);
        let (e0, e2) = 
        if l2 == r0 {
            (l0, r2)
        } else if l2 == r2 {
            (l0, r0)
        } else if l0 == r0 {
            (l2, r2)
        } else if l0 == r2 { 
            (l2, r0)
        } else {
            panic!()
        };

        if e0 == e2 { 
            TngComp::circ(*e1)
        } else if e0 < e2 { 
            TngComp::Arc(*e0, *e1, *e2)
        } else { 
            TngComp::Arc(*e2, *e1, *e0)
        }
    }
}

impl PartialEq for TngComp {
    fn eq(&self, other: &Self) -> bool {
        use TngComp::{Arc, Circ};
        match (self, other) {
            (Arc(l0, _, l2), Arc(r0, _, r2)) => l0 == r0 && l2 == r2,
            (Circ(l0), Circ(r0)) => l0 == r0,
            _ => false,
        }
    }
}

impl Display for TngComp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use TngComp::{Arc, Circ};
        match self { 
            Arc(e0, _, e2) => write!(f, "-{e0}-{e2}-"),
            Circ(e)        => write!(f, "○{}", yui_utils::subscript(*e as isize))
        }
    }
}

impl From<LinkComp> for TngComp {
    fn from(c: LinkComp) -> Self {
        let e = c.min_edge();
        if c.is_circle() { 
            Self::Circ(e)
        } else { 
            let (l, r) = c.ends().unwrap();
            Self::Arc(l, e, r)
        }
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Tng {
    comps: Vec<TngComp> // arc or circle
}

impl Tng { 
    pub fn new(comps: Vec<TngComp>) -> Self { 
        Self { comps }
    }

    pub fn res(x: &Crossing, r: Resolution) -> Self { 
        let (r0, r1) = x.res_arcs(r);
        let (c0, c1) = (TngComp::from(r0), TngComp::from(r1));
        Self::new(vec![c0, c1])
    }

    pub fn empty() -> Self { 
        Self::new(vec![])
    }

    pub fn is_empty(&self) -> bool { 
        self.comps.is_empty()
    }

    pub fn is_closed(&self) -> bool { 
        self.comps.iter().all(|a| a.is_circle())
    }

    pub fn ncomps(&self) -> usize { 
        self.comps.len()
    }

    pub fn comps(&self) -> &Vec<TngComp> { 
        &self.comps
    }

    pub fn comp(&self, i: usize) -> &TngComp { 
        &self.comps[i]
    }

    pub fn sub<I>(&self, iter: I) -> Self 
    where I: IntoIterator<Item = usize> { 
        Self::new(
            iter.into_iter().map(|i| self.comp(i).clone()).collect()
        )
    }

    pub fn endpts(&self) -> HashSet<Edge> { 
        self.comps.iter().flat_map(|c| 
            c.endpts().map(|(e0, e1)| vec![e0, e1]).unwrap_or(vec![])
        ).collect()
    }

    pub fn contains(&self, c: &TngComp) -> bool { 
        self.comps.contains(c)
    }

    pub fn index_of(&self, c: &TngComp) -> Option<usize> {
        self.comps.iter().position(|c1| c1 == c)
    }

    pub fn remove_at(&mut self, i: usize) -> TngComp {
        self.comps.remove(i)
    }

    pub fn connect(&mut self, other: Self) {
        for c in other.comps.into_iter() {
            if c.is_circle() { 
                self.comps.push(c);
            } else { 
                self.append_arc(c);
            }
        }
    }

    pub fn append_arc(&mut self, arc: TngComp) { 
        assert!(arc.is_arc());

        let n = self.ncomps();

        // If one end of `arc` is connectable:
        if let Some(i) = self.find_connectable(&arc, n) { 
            let mut c = self.comps[i].connect(&arc);

            // If the other end is also connectable to a different component:
            if let Some(j) = self.find_connectable(&c, i) { 
                let arc_j = self.comps.remove(j);
                c = c.connect(&arc_j);
            }

            self.comps[i] = c;
        } else { 
            self.comps.push(arc);
        }
    }

    fn find_connectable(&self, arc: &TngComp, j: usize) -> Option<usize> {
        let n = self.comps.len();

        (0..n).find(|&i| 
            i != j && self.comps[i].is_connectable(arc)
        )
    }

    pub fn find_loop(&self) -> Option<usize> {
        self.comps.iter().enumerate().find(|(_, c)| 
            c.is_circle()
        ).map(|(i, _)|
            i
        )
    }

    pub fn euler_num(&self) -> isize { 
        // NOTE: χ(arc) = 1, χ(circle) = 0. 
        self.comps.iter().filter(|c| 
            c.is_arc()
        ).count() as isize
    }
}

impl Display for Tng {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_empty() { 
            write!(f, "∅")
        } else { 
            let arcs = self.comps.iter().map(|a| 
                a.to_string()
            ).join(", ");
            write!(f, "{{{}}}", arcs)
        }
    }
}

impl From<TngComp> for Tng {
    fn from(c: TngComp) -> Self {
        Self::new(vec![c])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_connectable() { 
        let c0 = TngComp::arc(0, 1);
        let c1 = TngComp::arc(1, 2);
        let c2 = TngComp::arc(2, 3);
        let e  = TngComp::circ(4);

        assert!( c0.is_connectable(&c1));
        assert!(!c0.is_connectable(&c2));
        assert!(!c0.is_connectable(&e));

        assert!( c1.is_connectable(&c0));
        assert!( c1.is_connectable(&c2));
        assert!(!c1.is_connectable(&e));
    }

    #[test]
    fn connect_comp() { 
        let c0 = TngComp::arc(0, 1);
        let c1 = TngComp::arc(1, 2);
        let c2 = TngComp::arc(0, 2);

        assert_eq!(c0.connect(&c1), TngComp::arc(0, 2));
        assert_eq!(c1.connect(&c0), TngComp::arc(0, 2));
        assert_eq!(c0.connect(&c2), TngComp::arc(1, 2));
        assert_eq!(c2.connect(&c0), TngComp::arc(1, 2));
        assert_eq!(c0.connect(&c1).connect(&c2), TngComp::circ(0));
    }

    #[test]
    fn append_arc() { 
        let mut t = Tng::empty();
        assert_eq!(t.ncomps(), 0);

        t.append_arc(TngComp::arc(0,1)); // [0-1]
        assert_eq!(t.ncomps(), 1);

        t.append_arc(TngComp::arc(2,3)); // [0-1] [2-3]
        assert_eq!(t.ncomps(), 2);

        t.append_arc(TngComp::arc(1,2)); // [0-1-2-3]
        assert_eq!(t.ncomps(), 1);

        t.append_arc(TngComp::arc(0,3)); // [-0-1-2-3-]
        assert_eq!(t.ncomps(), 1);
    }

    #[test]
    fn connect() { 
        let mut t0 = Tng::new(vec![
            TngComp::arc(0,1),
            TngComp::arc(2,3),
            TngComp::circ(10),
        ]);

        let t1 = Tng::new(vec![
            TngComp::arc(1,2),
            TngComp::arc(3,4),
            TngComp::circ(11),
        ]);

        t0.connect(t1);

        assert_eq!(t0, Tng::new(vec![
            TngComp::arc(0,4), // [0,1,2,3,4] -> [0,4]
            TngComp::circ(10),
            TngComp::circ(11),
        ]));
    }

    #[test]
    fn deloop() { 
        let mut t = Tng::empty();
        assert_eq!(t.ncomps(), 0);
        assert_eq!(t.find_loop(), None);

        t.append_arc(TngComp::arc(0,1));
        assert_eq!(t.ncomps(), 1);
        assert_eq!(t.find_loop(), None);

        t.append_arc(TngComp::arc(2,3));
        assert_eq!(t.ncomps(), 2);
        assert_eq!(t.find_loop(), None);

        t.append_arc(TngComp::arc(2,3));
        assert_eq!(t.ncomps(), 2);
        assert_eq!(t.find_loop(), Some(1));

        t.remove_at(1);

        assert_eq!(t.ncomps(), 1);
        assert_eq!(t.find_loop(), None);
    }
}