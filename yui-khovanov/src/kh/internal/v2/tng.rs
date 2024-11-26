use std::collections::HashSet;
use std::fmt::Display;
use std::hash::Hash;
use delegate::delegate;
use itertools::Itertools;
use yui_link::{Edge, Crossing, Path};

#[derive(Debug, Clone, Eq)]
pub struct TngComp(Path);

impl From<Path> for TngComp {
    fn from(path: Path) -> Self {
        Self(path)
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
            #[call(ends)]
            pub fn endpts(&self) -> Option<(Edge, Edge)>;
            pub fn contains(&self, e: Edge) -> bool;
            pub fn min_edge(&self) -> Edge;
        }
    }

    pub fn path(&self) -> &Path { 
        &self.0
    }

    pub fn is_connectable(&self, other: &Self) -> bool { 
        self.0.is_connectable(&other.0)
    }

    pub fn connect(&mut self, other: Self) { 
        self.0.connect(other.0)
    }

    pub fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge { 
        let path = Path::new(
            self.0.edges().iter().map(|e| f(*e)),
            self.0.is_circle()
        );
        Self::from(path)
    }
}

impl Display for TngComp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl PartialEq for TngComp {
    fn eq(&self, other: &Self) -> bool {
        self.0.unori_eq(&other.0)
    }
}

impl PartialOrd for TngComp {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for TngComp {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.is_circle().cmp(&other.is_circle())
        .then_with(|| self.min_edge().cmp(&other.min_edge()))
    }
}

impl Hash for TngComp {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.is_circle().hash(state);
        self.min_edge().hash(state);
        self.len().hash(state);
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

    pub fn from_resolved(x: &Crossing) -> Self { 
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

    pub fn ncomps(&self) -> usize { 
        self.comps.len()
    }

    pub fn comps(&self) -> impl Iterator<Item = &TngComp> { 
        self.comps.iter()
    }

    pub fn comp(&self, i: usize) -> &TngComp { 
        &self.comps[i]
    }

    pub fn endpts(&self) -> HashSet<Edge> { 
        self.comps.iter().flat_map(|c| 
            c.endpts().map(|(e0, e1)| vec![e0, e1]).unwrap_or_default()
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

    pub fn connected(&self, other: &Self) -> Self { 
        let mut res = self.clone();
        res.connect(other.clone());
        res
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

    pub fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge { 
        Self::new(
            self.comps.iter().map(|c| c.convert_edges(&f))
        )
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
        assert_eq!(t.ncomps(), 0);

        t.append_arc(TngComp::arc([0, 1])); // [0-1]
        assert_eq!(t.ncomps(), 1);
        assert!(t.comp(0).is_arc());

        t.append_arc(TngComp::arc([2, 3])); // [0-1] [2-3]
        assert_eq!(t.ncomps(), 2);
        assert!(t.comp(0).is_arc());
        assert!(t.comp(1).is_arc());

        t.append_arc(TngComp::arc([1, 2])); // [0-1-2-3]
        assert_eq!(t.ncomps(), 1);
        assert!(t.comp(0).is_arc());

        t.append_arc(TngComp::arc([0, 3])); // [-0-1-2-3-]
        assert_eq!(t.ncomps(), 1);
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
        assert_eq!(t.ncomps(), 0);
        assert_eq!(t.find_comp(|c| c.is_circle()), None);

        t.append_arc(TngComp::arc([0, 1]));
        assert_eq!(t.ncomps(), 1);
        assert_eq!(t.find_comp(|c| c.is_circle()), None);

        t.append_arc(TngComp::arc([2, 3]));
        assert_eq!(t.ncomps(), 2);
        assert_eq!(t.find_comp(|c| c.is_circle()), None);

        t.append_arc(TngComp::arc([2, 3]));
        assert_eq!(t.ncomps(), 2);
        assert_eq!(t.find_comp(|c| c.is_circle()), Some(1));

        t.remove_at(1);

        assert_eq!(t.ncomps(), 1);
        assert_eq!(t.find_comp(|c| c.is_circle()), None);
    }

    #[test]
    fn convert_edges() { 
        let t = Tng::new(vec![
            TngComp::arc([0, 1]),
            TngComp::arc([2, 3]),
            TngComp::circ([10]),
        ]);
        let t = t.convert_edges(|e| 100 + e);
        assert_eq!(t, Tng::new(vec![
            TngComp::arc([100, 101]),
            TngComp::arc([102, 103]),
            TngComp::circ([110]),
        ]));
    }
}