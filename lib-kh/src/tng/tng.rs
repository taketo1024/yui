//! [`Tng`]: a planar tangle as a set of [`TngComp`]s — arcs and circles carrying
//! the edges they pass through.

use std::fmt::Display;
use itertools::Itertools;
use yui_core::bitmap::BitMap;
use yui_core::ext::CloneAnd;
use yui_link::{Edge, Node, Path};
use crate::util::CachedHash;

// Edge presence as a packed bitmap, `MAX_EDGE` the largest label it can hold.
// With the conventional 1-indexed consecutive numbering, u128 covers ≤ 63
// crossings (a 64-crossing link has edge 128, which does not fit).
cfg_if::cfg_if! {
    if #[cfg(feature = "big-link")] {
        use yui_core::u256::U256;
        type EdgeSet = BitMap<Edge, U256>;
        pub(crate) const MAX_EDGE: Edge = 255;
    } else {
        type EdgeSet = BitMap<Edge, u128>;
        pub(crate) const MAX_EDGE: Edge = 127;
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum TngCompKind {
    /// Arc with endpoint edges `e0 ≤ e1` (canonical).
    Arc { e0: Edge, e1: Edge },
    Circ,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct TngComp {
    kind: TngCompKind,
    edges: EdgeSet,
    marked: bool,
}

impl From<Path> for TngComp {
    fn from(path: Path) -> Self {
        Self::from_path(path, false)
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

    pub fn from_path(path: Path, marked: bool) -> Self {
        debug_assert!(path.edges().iter().all(|&e| e <= MAX_EDGE), "edge label exceeds the EdgeSet capacity ({MAX_EDGE})");
        let edges: EdgeSet = path.edges().iter().copied().collect();
        let kind = match &path {
            Path::Arc(es) => {
                let first = es[0];
                let last = *es.last().unwrap();
                TngCompKind::Arc { e0: first.min(last), e1: first.max(last) }
            }
            Path::Circ(_) => TngCompKind::Circ,
        };
        Self { kind, edges, marked }
    }

    pub fn is_marked(&self) -> bool {
        self.marked
    }

    pub fn is_arc(&self) -> bool {
        matches!(self.kind, TngCompKind::Arc { .. })
    }

    pub fn is_circle(&self) -> bool {
        matches!(self.kind, TngCompKind::Circ)
    }

    pub fn end_pts(&self) -> Option<(Edge, Edge)> {
        match self.kind {
            TngCompKind::Arc { e0, e1 } => Some((e0, e1)),
            TngCompKind::Circ => None,
        }
    }

    pub fn len(&self) -> usize {
        self.edges.len()
    }

    pub fn contains(&self, e: Edge) -> bool {
        self.edges.contains(e)
    }

    pub fn min_edge(&self) -> Edge {
        self.edges.iter().next().expect("empty TngComp has no min edge")
    }

    /// Iterate the edges this comp contains, in ascending order.
    pub fn edges(&self) -> impl ExactSizeIterator<Item = Edge> {
        self.edges.iter()
    }

    pub fn is_connectable(&self, other: &Self) -> bool {
        let Some((a, b)) = self.end_pts() else { return false };
        let Some((c, d)) = other.end_pts() else { return false };
        a == c || a == d || b == c || b == d
    }

    /// True iff both arcs have the same unordered endpoint pair (would close
    /// into a circle on `connect`).
    pub fn is_connectable_bothends(&self, other: &Self) -> bool {
        let Some(s) = self.end_pts() else { return false };
        let Some(o) = other.end_pts() else { return false };
        s == o
    }

    pub fn connect(&self, other: &Self) -> Self {
        debug_assert!(self.is_connectable(other), "{self} and {other} are not connectable.");

        let (TngCompKind::Arc { e0: a, e1: b }, TngCompKind::Arc { e0: c, e1: d }) =
            (self.kind, other.kind)
            else { panic!("connect requires two arcs") };

        let meets = (a == c) as u8 + (a == d) as u8 + (b == c) as u8 + (b == d) as u8;
        let new_kind = match meets {
            // Both endpoints match → closes into a circle.
            2 => TngCompKind::Circ,
            // One endpoint matches → arc with the two unmatched endpoints.
            1 => {
                let (free_self, free_other) =
                    if a == c { (b, d) }
                    else if a == d { (b, c) }
                    else if b == c { (a, d) }
                    else { /* b == d */ (a, c) };
                TngCompKind::Arc {
                    e0: free_self.min(free_other),
                    e1: free_self.max(free_other),
                }
            }
            _ => unreachable!("is_connectable guarantees meets ∈ {{1, 2}}"),
        };

        Self {
            kind: new_kind,
            edges: self.edges | other.edges,
            marked: self.marked || other.marked,
        }
    }

    pub fn connect_mut(&mut self, other: &Self) {
        *self = self.connect(other)
    }

    pub(crate) fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge {
        let edges: EdgeSet = self.edges().map(&f).collect();
        let kind = match self.kind {
            TngCompKind::Arc { e0, e1 } => {
                let (a, b) = (f(e0), f(e1));
                TngCompKind::Arc { e0: a.min(b), e1: a.max(b) }
            }
            TngCompKind::Circ => TngCompKind::Circ,
        };
        // `marked` preserved: InvLink's base_pt is on-axis (`f(b) == b`).
        Self { kind, edges, marked: self.marked }
    }
}

impl Display for TngComp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.marked {
            write!(f, "*")?;
        }
        let body = self.edges.iter().map(|e| e.to_string()).join("-");
        match self.kind {
            TngCompKind::Arc { .. } => write!(f, "[{body}]"),
            TngCompKind::Circ       => write!(f, "⚪︎({body})"),
        }
    }
}


// `CachedHash` caches the structural hash — `Tng` is hashed constantly via interning and `Arc<Tng>`
// in `CobComp`. Mutate only via `inner_mut`. `comps` stays sorted (`normalize`) so the derived `Ord`
// (required since `CobComp` derives `Ord` over `Arc<Tng>`) is canonical.
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug, Default)]
pub struct Tng {
    comps: CachedHash<Vec<TngComp>>, // arc or circle
}

impl Tng {
    pub fn new<I>(comps: I) -> Self
    where I: IntoIterator<Item = TngComp> {
        let comps = comps.into_iter().sorted().collect_vec();
        Self { comps: CachedHash::new(comps) }
    }

    pub fn from_resolved(x: &Node, base_pt: Option<Edge>) -> Self {
        assert!(x.is_resolved());

        let mk = |r: Path| -> TngComp {
            let marked = base_pt.map(|b| r.contains(b)).unwrap_or(false);
            TngComp::from_path(r, marked)
        };
        let (r0, r1) = x.arcs();
        let (c0, c1) = (mk(r0), mk(r1));

        if c0.is_connectable(&c1) {
            Self::from(c0.connect(&c1))
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

    pub fn end_pts(&self) -> impl Iterator<Item = Edge> + '_ {
        self.comps.iter()
            .filter_map(|c| c.end_pts())
            .flat_map(|(e0, e1)| [e0, e1])
    }

    pub fn contains(&self, c: &TngComp) -> bool {
        self.comps.binary_search(c).is_ok() // comps kept sorted (see normalize)
    }

    pub fn index_of(&self, c: &TngComp) -> Option<usize> {
        self.comps.binary_search(c).ok() // comps kept sorted (see normalize)
    }

    pub fn remove(&mut self, c: &TngComp) -> TngComp {
        let i = self.index_of(c).expect("component not found");
        self.comps.inner_mut().remove(i)
    }

    pub fn connect(&self, other: &Self) -> Self {
        self.clone_and(|t| t.connect_mut(other))
    }

    pub fn connect_mut(&mut self, other: &Self) {
        for c in other.comps.iter() {
            if c.is_circle() {
                self.comps.inner_mut().push(c.clone());
            } else {
                self.append_arc(c.clone());
            }
        }
        self.normalize();
    }

    pub fn append_arc(&mut self, arc: TngComp) {
        debug_assert!(arc.is_arc());

        // If one end of `arc` is connectable:
        if let Some(i) = self.find_comp(|c| c.is_connectable(&arc)) {
            self.comps.inner_mut()[i].connect_mut(&arc);

            // If the other end is also connectable to a different component:
            let ci = self.comps[i].clone();
            if let Some(j) = self.find_comp(|c| *c != ci && c.is_connectable(&ci)) {
                let cj = self.comps.inner_mut().remove(j);
                self.comps.inner_mut()[i].connect_mut(&cj);
            }
        } else {
            self.comps.inner_mut().push(arc);
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
        self.comps.inner_mut().sort()
    }

    pub(crate) fn convert_edges<F>(&self, f: F) -> Self
    where F: Fn(Edge) -> Edge {
        Self::new(self.comps.iter().map(|c| c.convert_edges(&f)))
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
    use std::hash::{Hash, Hasher};
    use rustc_hash::FxHasher;
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
        let c0 = TngComp::arc([0, 1]);
        let c1 = TngComp::arc([1, 2]);
        let c2 = TngComp::arc([0, 2]);

        let c0 = c0.connect(&c1);
        assert_eq!(c0, TngComp::arc([0, 1, 2]));

        let c0 = c0.connect(&c2);
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

        t0.connect_mut(&t1);

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
    fn tng_hash_cache_invalidates() {
        fn hash_of(t: &Tng) -> u64 {
            let mut h = FxHasher::default();
            t.hash(&mut h);
            h.finish()
        }

        // remove: force the cache, mutate, then the hash must match a freshly-built Tng.
        let mut t = Tng::new([TngComp::arc([0, 1]), TngComp::arc([2, 3])]);
        let _ = hash_of(&t); // populate cache
        t.remove(&TngComp::arc([0, 1]));
        let fresh = Tng::new([TngComp::arc([2, 3])]);
        assert_eq!(hash_of(&t), hash_of(&fresh), "stale hash after remove");
        assert_eq!(t, fresh);

        // connect_mut: ends in normalize → invalidate.
        let mut a = Tng::new([TngComp::arc([0, 1])]);
        let _ = hash_of(&a);
        a.connect_mut(&Tng::new([TngComp::arc([2, 3])]));
        let fresh2 = Tng::new([TngComp::arc([0, 1]), TngComp::arc([2, 3])]);
        assert_eq!(hash_of(&a), hash_of(&fresh2), "stale hash after connect_mut");
        assert_eq!(a, fresh2);
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

        let c = t.comp(1).clone();
        t.remove(&c);

        assert_eq!(t.n_comps(), 1);
        assert_eq!(t.find_comp(|c| c.is_circle()), None);
    }

}