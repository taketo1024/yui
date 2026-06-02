use std::collections::HashSet;
use std::fmt::Display;
use std::hash::Hash;
use delegate::delegate;
use itertools::Itertools;
use smallvec::SmallVec;
use yui_link::{Edge, Node, Path};

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct TngComp(Path);

impl From<Path> for TngComp {
    /// Wraps a [`Path`] and normalizes it, so that two `TngComp`s compare
    /// equal iff their backing paths are unoriented-equivalent.
    fn from(path: Path) -> Self {
        let mut c = Self(path);
        c.normalize();
        c
    }
}

impl TngComp {
    pub fn arc<I>(edges: I) -> Self
    where I: IntoIterator<Item = Edge> {
        Self::from(Path::arc(edges))
    }

    pub fn circ<I>(edges: I) -> Self
    where I: IntoIterator<Item = Edge> {
        Self::from(Path::circ(edges))
    }

    delegate! {
        to self.0 {
            pub fn len(&self) -> usize;
            pub fn is_arc(&self) -> bool;
            pub fn is_circle(&self) -> bool;
            pub fn end_pts(&self) -> Option<(Edge, Edge)>;
            pub fn contains(&self, e: Edge) -> bool;
            pub fn min_edge(&self) -> Edge;
        }
    }

    /// Rewrite the backing [`Path`] in canonical form so two `TngComp`s
    /// compare equal iff their oriented paths agree as *unoriented* sequences.
    /// - `Arc(es)`:  reverse if `es[0] > es[last]`.
    /// - `Circ(es)`: rotate so `es[0] == min(es)`, then reflect the suffix so
    ///               `es[1] ≤ es[last]`.
    fn normalize(&mut self) {
        match &mut self.0 {
            Path::Arc(es) => {
                if es.len() >= 2 && *es.last().unwrap() < es[0] {
                    es.reverse();
                }
            }
            Path::Circ(es) => {
                let n = es.len();
                if n > 1 {
                    let min_idx = (0..n).min_by_key(|&i| es[i]).unwrap();
                    es.rotate_left(min_idx);
                    if n > 2 && es[n - 1] < es[1] {
                        es[1..].reverse();
                    }
                }
            }
        }
    }

    pub fn is_connectable(&self, other: &Self) -> bool {
        let Some((e0, e1)) = self.end_pts() else { return false };
        let Some((f0, f1)) = other.end_pts() else { return false };
        e0 == f0 || e0 == f1 || e1 == f0 || e1 == f1
    }

    /// Endpoint-set match between two arcs. Both `TngComp`s are normalized,
    /// so endpoint pairs are sorted — direct equality suffices.
    pub fn is_connectable_bothends(&self, other: &Self) -> bool {
        let Some(s) = self.end_pts() else { return false };
        let Some(o) = other.end_pts() else { return false };
        s == o
    }

    pub fn connect(&mut self, other: Self) {
        assert!(self.is_connectable(&other), "{self} and {other} are not connectable.");

        let placeholder = Path::Arc(SmallVec::new());
        let this = std::mem::replace(&mut self.0, placeholder);
        let (mut left, right) = (this.into_seq(), other.0.into_seq());

        let (e0, e1) = (left[0], *left.last().unwrap());
        let (f0, f1) = (right[0], *right.last().unwrap());

        let mut combined = if e1 == f0 {
            // self ++ other[1..]
            left.extend_from_slice(&right[1..]);
            left
        } else if e1 == f1 {
            // self ++ reverse(other[..-1])
            let mut r = right;
            r.pop();
            r.reverse();
            left.extend(r);
            left
        } else if e0 == f0 {
            // reverse(other[1..]) ++ self
            let mut r = right;
            r.remove(0);
            r.reverse();
            r.append(&mut left);
            r
        } else {
            // e0 == f1: other[..-1] ++ self
            let mut r = right;
            r.pop();
            r.append(&mut left);
            r
        };

        // If the new endpoints coincide, the result closes into a circle.
        let closes = combined.len() > 1 && combined[0] == *combined.last().unwrap();
        self.0 = if closes {
            combined.pop();
            Path::Circ(combined)
        } else {
            Path::Arc(combined)
        };
        self.normalize();
    }

    pub fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge {
        let mapped = self.0.edges().iter().map(|e| f(*e));
        let path = if self.0.is_circle() {
            Path::circ(mapped)
        } else {
            Path::arc(mapped)
        };
        Self::from(path)
    }
}

impl Display for TngComp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}


#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Tng {
    comps: Vec<TngComp> // arc or circle
}

impl Tng { 
    pub fn new<I>(comps: I) -> Self 
    where I: IntoIterator<Item = TngComp> { 
        let comps = comps.into_iter().sorted().collect_vec();
        Self { comps }
    }

    pub fn from_resolved(x: &Node) -> Self { 
        assert!(x.is_resolved());

        let (r0, r1) = x.arcs();
        let (mut c0, c1) = (
            TngComp::from(r0), 
            TngComp::from(r1)
        );

        if c0.is_connectable(&c1) { 
            c0.connect(c1);
            Self::from(c0)
        } else { 
            Self::new([c0, c1])
        }
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

    pub fn contains_circle(&self) -> bool { 
        self.comps.iter().any(|a| a.is_circle())
    }

    pub fn n_comps(&self) -> usize { 
        self.comps.len()
    }

    pub fn comps(&self) -> impl Iterator<Item = &TngComp> { 
        self.comps.iter()
    }

    pub fn comp(&self, i: usize) -> &TngComp { 
        &self.comps[i]
    }

    pub fn end_pts(&self) -> HashSet<Edge> {
        self.end_pts_iter().collect()
    }

    /// Alloc-free iteration over boundary 0-cells. Walks each arc and yields
    /// its two endpoints; circles contribute nothing.
    pub fn end_pts_iter(&self) -> impl Iterator<Item = Edge> + '_ {
        self.comps.iter()
            .filter_map(|c| c.end_pts())
            .flat_map(|(e0, e1)| [e0, e1])
    }

    /// True iff `e` is an endpoint of some arc in `self`.
    pub fn contains_endpoint(&self, e: Edge) -> bool {
        self.comps.iter().any(|c|
            matches!(c.end_pts(), Some((a, b)) if a == e || b == e)
        )
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
        self.normalize();
    }

    pub fn append_arc(&mut self, arc: TngComp) { 
        assert!(arc.is_arc());

        // If one end of `arc` is connectable:
        if let Some(i) = self.find_comp(|c| c.is_connectable(&arc)) { 
            self.comps[i].connect(arc);

            // If the other end is also connectable to a different component:
            let ci = &self.comps[i];
            if let Some(j) = self.find_comp(|c| c != ci && c.is_connectable(ci)) { 
                let cj = self.comps.remove(j);
                self.comps[i].connect(cj);
            }
        } else { 
            self.comps.push(arc);
        }

        self.normalize();
    }

    pub fn find_comp<F>(&self, pred: F) -> Option<usize>
    where F: Fn(&TngComp) -> bool {
        self.comps.iter().enumerate().find(|(_, c)| 
            pred(c)
        ).map(|(i, _)| i)
    }

    pub fn euler_num(&self) -> isize { 
        // NOTE: χ(arc) = 1, χ(circle) = 0. 
        self.comps.iter().filter(|c| 
            c.is_arc()
        ).count() as isize
    }

    fn normalize(&mut self) { 
        self.comps.sort()
    }

}

impl Display for Tng {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_empty() { 
            write!(f, "∅")
        } else if self.comps.len() == 1 { 
            write!(f, "{}", self.comps[0])
        } else { 
            write!(f, "{{{}}}", self.comps.iter().join(", "))
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
    fn tng_comp_eq() { 
        assert_eq!(TngComp::arc([0, 1, 2]), TngComp::arc([0, 1, 2]));
        assert_eq!(TngComp::arc([0, 1, 2]), TngComp::arc([2, 1, 0]));
        assert_ne!(TngComp::arc([0, 1, 2]), TngComp::arc([0, 2]));
        assert_ne!(TngComp::arc([0, 1, 2]), TngComp::circ([0, 1, 2]));
        assert_eq!(TngComp::circ([0, 1, 2]), TngComp::circ([0, 1, 2]));
        assert_eq!(TngComp::circ([0, 1, 2]), TngComp::circ([1, 2, 0]));
        assert_eq!(TngComp::circ([0, 1, 2]), TngComp::circ([2, 1, 0]));
        assert_ne!(TngComp::circ([0, 1, 2]), TngComp::circ([1, 2, 3]));
        assert_ne!(TngComp::circ([0, 1, 2]), TngComp::circ([0, 1, 2, 3]));
    }

    #[test]
    fn is_connectable() { 
        let c0 = TngComp::arc([0, 1]);
        let c1 = TngComp::arc([1, 2]);
        let c2 = TngComp::arc([2, 3]);
        let e  = TngComp::circ([4]);

        assert!( c0.is_connectable(&c1));
        assert!(!c0.is_connectable(&c2));
        assert!(!c0.is_connectable(&e));

        assert!( c1.is_connectable(&c0));
        assert!( c1.is_connectable(&c2));
        assert!(!c1.is_connectable(&e));
    }

    #[test]
    fn connect_comp() { 
        let mut c0 = TngComp::arc([0, 1]);
        let c1 = TngComp::arc([1, 2]);
        let c2 = TngComp::arc([0, 2]);

        c0.connect(c1);
        assert_eq!(c0, TngComp::arc([0, 1, 2]));

        c0.connect(c2);
        assert_eq!(c0, TngComp::circ([0, 1, 2]));
    }

    #[test]
    fn append_arc() { 
        let mut t = Tng::empty();
        assert_eq!(t.n_comps(), 0);

        t.append_arc(TngComp::arc([0, 1])); // [0-1]
        assert_eq!(t.n_comps(), 1);
        assert!(t.comp(0).is_arc());

        t.append_arc(TngComp::arc([2, 3])); // [0-1] [2-3]
        assert_eq!(t.n_comps(), 2);
        assert!(t.comp(0).is_arc());
        assert!(t.comp(1).is_arc());

        t.append_arc(TngComp::arc([1, 2])); // [0-1-2-3]
        assert_eq!(t.n_comps(), 1);
        assert!(t.comp(0).is_arc());

        t.append_arc(TngComp::arc([0, 3])); // [-0-1-2-3-]
        assert_eq!(t.n_comps(), 1);
        assert!(t.comp(0).is_circle());
    }

    #[test]
    fn connect() { 
        let mut t0 = Tng::new(vec![
            TngComp::arc([0, 1]),
            TngComp::arc([2, 3]),
            TngComp::circ([10]),
        ]);

        let t1 = Tng::new(vec![
            TngComp::arc([1, 2]),
            TngComp::arc([3, 4]),
            TngComp::circ([11]),
        ]);

        t0.connect(t1);

        assert_eq!(t0, Tng::new(vec![
            TngComp::arc([0, 1, 2, 3, 4]),
            TngComp::circ([10]),
            TngComp::circ([11]),
        ]));
    }

    #[test]
    fn tng_eq() { 
        let t0 = Tng::new(vec![
            TngComp::arc([0, 1]),
            TngComp::arc([2, 3]),
        ]);

        let t1 = Tng::new(vec![
            TngComp::arc([2, 3]),
            TngComp::arc([0, 1]),
        ]);

        assert_eq!(t0, t1);
    }

    #[test]
    fn find_loop() { 
        let mut t = Tng::empty();
        assert_eq!(t.n_comps(), 0);
        assert_eq!(t.find_comp(|c| c.is_circle()), None);

        t.append_arc(TngComp::arc([0, 1]));
        assert_eq!(t.n_comps(), 1);
        assert_eq!(t.find_comp(|c| c.is_circle()), None);

        t.append_arc(TngComp::arc([2, 3]));
        assert_eq!(t.n_comps(), 2);
        assert_eq!(t.find_comp(|c| c.is_circle()), None);

        t.append_arc(TngComp::arc([2, 3]));
        assert_eq!(t.n_comps(), 2);
        assert_eq!(t.find_comp(|c| c.is_circle()), Some(1));

        t.remove_at(1);

        assert_eq!(t.n_comps(), 1);
        assert_eq!(t.find_comp(|c| c.is_circle()), None);
    }

}