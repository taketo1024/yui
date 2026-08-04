use core::panic;
use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use itertools::Itertools;

use super::{Node, NodeOri, Path};

#[cfg(not(feature = "big-link"))]
pub type Edge = u8;
#[cfg(not(feature = "big-link"))]
pub type StateRepr = u64;

// `big-link`: diagrams with up to 128 crossings (256 edge labels).
#[cfg(feature = "big-link")]
pub type Edge = u16;
#[cfg(feature = "big-link")]
pub type StateRepr = u128;

pub type State = yui_core::bitseq::BitSeq<StateRepr>;

#[derive(Debug, Clone)]
pub struct Link {
    nodes: Vec<Node>,
    loops: Vec<Edge>,
    base_pt: Option<Edge>,
}

impl Link {
    /// The maximum number of crossings a `Link` can carry, bounded by the `State` width
    /// (64 by default, 128 under the `big-link` feature).
    pub const MAX_CROSSING: usize = State::MAX_LEN;

    pub fn new(
        nodes: impl IntoIterator<Item = Node>,
        loops: impl IntoIterator<Item = Edge>,
    ) -> Self {
        let nodes = nodes.into_iter().collect_vec();
        let loops = loops.into_iter().collect_vec();

        assert!(
            nodes.len() <= Self::MAX_CROSSING,
            "too many crossings: {} > MAX_CROSSING = {} (enable the `big-link` feature for up to 128)",
            nodes.len(), Self::MAX_CROSSING
        );

        let edge_counts = nodes.iter().flat_map(|x| x.edges()).cloned().counts();
        assert!(
            edge_counts.values().all(|&c| c == 2),
            "Invalid data: each edge in the diagram must appear exactly twice."
        );

        let node_edges: HashSet<Edge> = edge_counts.into_keys().collect();
        let mut loop_set: HashSet<Edge> = HashSet::new();
        for &e in &loops {
            assert!(!node_edges.contains(&e), "loop edge {e} is already used in a node");
            assert!(loop_set.insert(e), "duplicate loop edge: {e}");
        }

        // Default base_pt to the minimal edge (if any).
        let base_pt = node_edges.iter().chain(loops.iter()).copied().min();

        Self { nodes, loops, base_pt }
    }

    pub fn from_nodes(nodes: impl IntoIterator<Item = Node>) -> Self {
        Self::new(nodes, [])
    }

    // Re-derive each crossing's orientation by traversing components. `is_incoming(i, j)` tells
    // whether port j of node i is known to receive an incoming strand (PD codes: j == 0). The first
    // claimed port met by a tentative traversal fixes the component's direction; a component claiming
    // no port is undetermined (cf. `unlink2`) and the whole link is left unoriented. A fixed direction
    // contradicting `is_incoming` (an odd PD code) panics. Returns whether the link is now oriented.
    pub(crate) fn reorient<F>(&mut self, is_incoming: F) -> bool
    where F: Fn(usize, usize) -> bool {
        use crate::NodeOri::None;

        let mut incoming: Vec<Vec<usize>> = vec![vec![]; self.n_nodes()];
        let mut remain: HashSet<Edge> = self.nodes.iter().flat_map(|x| x.edges().iter().copied()).collect();
        let mut undetermined = false;

        while !remain.is_empty() {
            // start at a claimed port of an untraversed component, so the direction is correct
            // from the outset. Components claiming no port are undetermined (cf. `unlink2`).
            let Some(start) = self.find_port(|i, j|
                remain.contains(&self.node(i).edge(j)) && is_incoming(i, j)
            ) else {
                undetermined = true;
                break;
            };

            self.traverse_from(start, |i, j| {
                remain.remove(&self.node(i).edge(j));
                assert!(
                    is_incoming(i, j) || !is_incoming(i, (j + 2) % 4),
                    "inconsistent orientation: the strand through node {i} exits at port {}, which is claimed incoming", (j + 2) % 4
                );
                incoming[i].push(j);
            });
        }

        // a crossing's two incoming ports fix the orientation (see `NodeOri::from_in_ports`);
        // if any node is incoherent, or some component is undetermined, the whole link is unoriented.
        let oris = incoming.iter().map(|ports| match ports[..] {
            [p, q] => NodeOri::from_in_ports(p, q),
            _ => None,
        }).collect_vec();
        let coherent = !undetermined && !oris.contains(&None);
        
        self.nodes.iter_mut().zip(oris).for_each(|(n, o)| 
            n.ori = if coherent { o } else { None }
        );

        coherent
    }

    // note: builder links may have an edge with both ends at port 2, so no port is excluded here;
    // PD-specific constraints (never enter at 2) belong in the caller's predicate.
    fn find_port(&self, f: impl Fn(usize, usize) -> bool) -> Option<(usize, usize)> {
        let n = self.n_nodes();
        (0..n).flat_map(|i|
            (0..4).map(move |j| (i, j))
        ).find(|&(i, j)|
            f(i, j)
        )
    }

    pub fn with_base_pt(mut self, e: Edge) -> Self {
        let exists = self.nodes.iter().any(|x| x.edges().contains(&e))
            || self.loops.contains(&e);
        assert!(exists, "base_pt {e} is not an edge of this link");
        self.base_pt = Some(e);
        self
    }

    pub fn base_pt(&self) -> Option<Edge> {
        self.base_pt
    }

    pub fn empty() -> Link {
        Self::new([], [])
    }

    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty() && self.loops.is_empty()
    }

    pub fn unknot() -> Link {
        Self::unlink(1)
    }

    pub fn unlink(n: usize) -> Link {
        Self::new([], (1..=n).map(|e| e as Edge))
    }

    pub fn is_knot(&self) -> bool { 
        self.n_comps() == 1
    }

    pub fn is_oriented(&self) -> bool { 
        self.nodes().all(|n| n.is_oriented())
    }

    pub fn writhe(&self) -> i32 { 
        let (p, n) = self.n_signed_crossings();
        (p as i32) - (n as i32)
    }

    pub fn n_nodes(&self) -> usize { 
        self.nodes.len()
    }

    pub fn nodes(&self) -> impl Iterator<Item = &Node> { 
        self.nodes.iter()
    }

    pub fn node(&self, i: usize) -> &Node { 
        &self.nodes[i]
    }

    pub(crate) fn node_mut(&mut self, i: usize) -> &mut Node { 
        &mut self.nodes[i]
    }

    pub fn crossings(&self) -> impl Iterator<Item = &Node> { 
        self.nodes.iter().filter(|x| x.is_crossing())
    }

    pub fn n_crossings(&self) -> usize { 
        self.nodes.iter()
            .filter(|x| x.is_crossing())
            .count()
    }

    pub fn n_signed_crossings(&self) -> (usize, usize) {
        let mut pos = 0;
        let mut neg = 0;
        for n in self.nodes.iter() { 
            if n.is_pos() { pos += 1 } 
            else if n.is_neg() { neg += 1}
        }
        (pos, neg)
    }

    pub fn loops(&self) -> &[Edge] {
        &self.loops
    }

    pub fn n_loops(&self) -> usize {
        self.loops.len()
    }

    pub fn n_edges(&self) -> usize {
        self.nodes.len() * 2 + self.loops.len()
    }

    pub fn edges(&self) -> Vec<Edge> {
        let mut edges: Vec<Edge> = self.nodes.iter()
            .flat_map(|x| x.edges().iter().copied())
            .chain(self.loops.iter().copied())
            .collect();
        edges.sort();
        edges.dedup();
        edges
    }

    pub fn n_comps(&self) -> usize {
        let mut count = 0;
        self.traverse_comps(|c, _, _|
            if count <= c { count = c + 1 }
        );
        count + self.loops.len()
    }

    pub fn comps(&self) -> Vec<Path> {
        let mut comps = vec![];

        self.traverse_comps(|c, i, j| {
            if c == comps.len() {
                comps.push(vec![]);
            }

            let e = self.node(i).edge(j);
            comps[c].push(e);
        });

        let mut result: Vec<Path> = comps.into_iter().map(Path::circ).collect();
        for &e in &self.loops {
            result.push(Path::circ(vec![e]));
        }
        result
    }

    pub fn traverse_comps<F>(&self, mut f: F) where 
    F: FnMut(usize, usize, usize) { 
        let mut c = 0; // component counter
        let mut remain: HashSet<Edge> = self.nodes.iter().flat_map(|x| x.edges().iter().copied()).collect();

        while !remain.is_empty() {
            // Take minimal edge-id. 
            let e0 = remain.iter().min().cloned().unwrap();

            // Find node & point having edge e0. 
            let (i0, j0) = self.find_port(|i, j| 
                self.node(i).edge(j) == e0
            ).unwrap();

            self.traverse_from((i0, j0), |i, j| { 
                let e = self.node(i).edge(j);
                remain.remove(&e);
                f(c, i, j);
            });

            // Onto next component.
            c += 1;
        }
    }

    pub fn traverse_from<F>(&self, start: (usize, usize), mut f:F) where
        F: FnMut(usize, usize)
    {
        let (mut i, mut j) = start;

        f(i, j); // call starting point

        loop {
            let c = self.node(i);
            let k = c.counter_pos(j);
            let next = self.traverse_outer(i, k);

            if next == start {
                break
            }

            (i, j) = next;

            f(i, j)
        }
    }

    fn traverse_outer(&self, n_index: usize, e_index: usize) -> (usize, usize) {
        let e = self.nodes[n_index].edge(e_index);

        for (i, c) in self.nodes.iter().enumerate() {
            for (j, &f) in c.edges().iter().enumerate() {
                if e == f && (n_index != i || (n_index == i && e_index != j)) {
                    return (i, j)
                }
            }
        }

        panic!("Broken data")
    }

    // Renumber the edges base..base+n in the order they are met traversing from `start_edge`
    // (a knot's one traversal covers every edge), keeping the diagram. base_pt becomes `base`
    // (base = 1 gives the usual 1-based numbering of knot theory).
    pub fn reindexed(&self, start_edge: Edge, base: Edge) -> Link {
        assert!(self.is_knot() && self.loops.is_empty(), "reindexed expects a knot");
        let start = self.find_port(|i, j|
            self.node(i).edge(j) == start_edge
        ).expect("start_edge must be an edge of the link");

        let mut map: HashMap<Edge, Edge> = HashMap::new();
        let mut next = base;
        self.traverse_from(start, |i, j| {
            map.entry(self.node(i).edge(j)).or_insert_with(|| {
                let id = next;
                next += 1;
                id
            });
        });

        let nodes = self.nodes.iter().map(|x|
            x.convert_edges(|e| map[&e])
        );
        Link::new(nodes, []).with_base_pt(base)
    }

    // The two (node, slot) ends of edge `e`. When `directed`, they are ordered as (tail, head)
    // along the orientation — the strand exits at the tail and enters at the head (cf.
    // `NodeOri::in_ports`); otherwise the order carries no meaning.
    pub(crate) fn edge_ends(&self, e: Edge, directed: bool) -> ((usize, usize), (usize, usize)) {
        assert!(!directed || self.is_oriented(), "directed edge_ends requires an oriented link");

        let (x, y) = self.nodes().enumerate().flat_map(|(i, n)|
            (0..4).filter(move |&s| n.edge(s) == e).map(move |s| (i, s))
        ).collect_tuple().unwrap_or_else(||
            panic!("edge {e} must appear exactly twice")
        );
        if !directed {
            return (x, y);
        }

        let is_in = |(i, s): (usize, usize)| {
            let ports = self.node(i).ori().in_ports().expect("directed edge_ends requires an oriented link");
            ports.contains(&s)
        };
        debug_assert!(is_in(x) != is_in(y), "edge {e} must have one head and one tail");
        if is_in(x) { (y, x) } else { (x, y) }
    }
}

impl Display for Link {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "L[{}]", self.nodes.iter().map(|x| x.to_string()).join(", "))
    }
}

#[cfg(test)]
mod tests {
    use yui_core::bitseq::Bit;

    use super::*;

    #[test]
    fn link_init() {
        let l = Link::from_nodes(vec![]);
        assert_eq!(l.nodes.len(), 0);
    }

    #[test]
    fn link_is_empty() {
        let l = Link::empty();
        assert!(l.is_empty());

        let l = Link::test_data("unknot_l_twist");
        assert!(!l.is_empty());
    }

    #[test]
    fn link_crossing_num() {
        let l = Link::empty();
        assert_eq!(l.n_crossings(), 0);

        let l = Link::test_data("unknot_l_twist");
        assert_eq!(l.n_crossings(), 1);
        
        let l = Link::test_data("3_1");
        assert_eq!(l.n_crossings(), 3);
    }

    #[test]
    fn link_next() {
        let l = Link::test_data("unknot_l_twist");

        assert_eq!(l.traverse_outer(0, 0), (0, 1));
        assert_eq!(l.traverse_outer(0, 1), (0, 0));
        assert_eq!(l.traverse_outer(0, 2), (0, 3));
        assert_eq!(l.traverse_outer(0, 3), (0, 2));
    }

    #[test]
    fn link_traverse() {
        let traverse = |l: &Link, (i0, j0)| { 
            let mut queue = vec![];
            l.traverse_from((i0, j0), |i, j| queue.push((i, j)));
            queue
        };

        let l = Link::test_data("unknot_l_twist");
        let path = traverse(&l, (0, 0));
        
        assert_eq!(path, [(0, 0), (0, 3)]); // loop
    }

    #[test]
    fn link_crossing_signs() {
        let l = Link::test_data("unknot_l_twist");
        assert_eq!(l.n_signed_crossings(), (1, 0));

        let l = Link::test_data("unknot_r_twist");
        assert_eq!(l.n_signed_crossings(), (0, 1));

        let l = Link::test_data("unknot_l_twist").resolve_at(0, Bit::Bit0);
        assert_eq!(l.n_signed_crossings(), (0, 0));
    }

    #[test]
    fn link_writhe() {
        let l = Link::test_data("unknot_l_twist");
        assert_eq!(l.writhe(), 1);

        let l = Link::test_data("unknot_r_twist");
        assert_eq!(l.writhe(), -1);

        let l = Link::test_data("unknot_l_twist").resolve_at(0, Bit::Bit0);
        assert_eq!(l.writhe(), 0);
    }

    #[test]
    fn link_components() {
        let l = Link::test_data("unknot_l_twist");
        let comps = l.comps();
        assert_eq!(comps, vec![ Path::circ(vec![1, 2])]);
    }



    #[test]
    fn empty_link() {
        let l = Link::empty();
        assert_eq!(l.n_crossings(), 0);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.n_comps(), 0);
    }

    #[test]
    fn unknot() {
        let l = Link::unknot();

        assert!(!l.is_empty());
        assert!(l.is_oriented());
        assert!(l.is_knot());

        assert_eq!(l.n_crossings(), 0);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.n_edges(), 1);
        assert_eq!(l.n_comps(), 1);
        assert_eq!(l.n_loops(), 1);

        assert_eq!(l.loops(), &[1]);
        assert_eq!(l.comps(), vec![Path::circ(vec![1])]);
    }

    #[test]
    fn unlink_zero() {
        let l = Link::unlink(0);

        assert!(l.is_empty());
        assert!(l.is_oriented());
        assert!(!l.is_knot());

        assert_eq!(l.n_crossings(), 0);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.n_edges(), 0);
        assert_eq!(l.n_comps(), 0);
        assert_eq!(l.n_loops(), 0);

        assert_eq!(l.loops(), &[] as &[Edge]);
        assert_eq!(l.comps(), vec![] as Vec<Path>);
    }

    #[test]
    fn unlink_n() {
        let l = Link::unlink(3);

        assert!(!l.is_empty());
        assert!(l.is_oriented());
        assert!(!l.is_knot());

        assert_eq!(l.n_crossings(), 0);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.n_edges(), 3);
        assert_eq!(l.n_comps(), 3);
        assert_eq!(l.n_loops(), 3);

        assert_eq!(l.loops(), &[1, 2, 3]);
        assert_eq!(
            l.comps(),
            vec![Path::circ(vec![1]), Path::circ(vec![2]), Path::circ(vec![3])]
        );
    }


    #[test]
    fn trefoil() {
        let l = Link::test_data("3_1");
        assert_eq!(l.n_crossings(), 3);
        assert_eq!(l.writhe(), -3);
        assert_eq!(l.n_comps(), 1);
    }

    #[test]
    fn figure8() {
        let l = Link::test_data("4_1");
        assert_eq!(l.n_crossings(), 4);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.n_comps(), 1);
    }

    #[test]
    fn hopf_link() {
        let l = Link::test_data("L2a1");
        assert_eq!(l.n_crossings(), 2);
        assert_eq!(l.writhe(), -2);
        assert_eq!(l.n_comps(), 2);
    }

    #[test]
    fn unlink_2() {
        // the over-component has no under-anchor, so the PD code leaves the link unoriented.
        let l = Link::test_data("unlink2");
        assert_eq!(l.n_crossings(), 2);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.n_comps(), 2);
        assert!(!l.is_oriented());
    }

    #[test]
    fn unlink_2_r2() {
        // the R2 pair alternates over/under, so both components are under-anchored.
        let l = Link::test_data("unlink2_r2");
        assert_eq!(l.n_crossings(), 2);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.n_comps(), 2);
        assert!(l.is_oriented());
    }

    #[test]
    fn l2x4() {
        let l = Link::test_data("L4a1");
        assert_eq!(l.n_crossings(), 4);
        assert_eq!(l.writhe(), -4);
        assert_eq!(l.n_comps(), 2);
    }


    #[test]
    fn base_pt_default_min_edge() {
        // Defaults to the minimal edge of the link.
        let l = Link::test_data("3_1");
        assert_eq!(l.base_pt(), Some(1));

        // Empty link has no edge, so base_pt is None.
        assert_eq!(Link::empty().base_pt(), None);
    }

    #[test]
    fn with_base_pt_sets_base_pt() {
        let l = Link::test_data("3_1").with_base_pt(1);
        assert_eq!(l.base_pt(), Some(1));
    }

    #[test]
    fn with_base_pt_on_loop() {
        let l = Link::unlink(3).with_base_pt(2);
        assert_eq!(l.base_pt(), Some(2));
    }


    #[test]
    #[should_panic]
    fn with_base_pt_invalid_panics() {
        // Edge 99 is not in the trefoil's edge set.
        let _ = Link::test_data("3_1").with_base_pt(99);
    }


}