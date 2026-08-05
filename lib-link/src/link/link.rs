use core::panic;
use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use itertools::Itertools;

use super::{Node, Path, Slot};

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

#[derive(Debug, Clone, PartialEq, Eq)]
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
        let bad = edge_counts.iter()
            .filter(|&(_, &c)| c != 2)
            .map(|(&e, &c)| (e, c))
            .sorted()
            .collect_vec();
        assert!(bad.is_empty(), "each edge must appear exactly twice; (edge, count) = {bad:?}");

        let node_edges: HashSet<Edge> = edge_counts.into_keys().collect();
        let mut loop_set: HashSet<Edge> = HashSet::new();
        for &e in &loops {
            assert!(!node_edges.contains(&e), "loop edge {e} is already used in a node");
            assert!(loop_set.insert(e), "duplicate loop edge: {e}");
        }

        // Default base_pt to the minimal edge (if any).
        let base_pt = node_edges.iter().chain(loops.iter()).copied().min();

        let l = Self { nodes, loops, base_pt };
        l.verify_ori();
        l
    }

    pub fn from_nodes(nodes: impl IntoIterator<Item = Node>) -> Self {
        Self::new(nodes, [])
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

    // An orientation belongs to the whole diagram: once any node has lost it, drop it everywhere.
    pub(crate) fn normalize_ori(&mut self) {
        if !self.is_oriented() {
            self.nodes.iter_mut().for_each(|n|
                n.set_incoming(None)
            );
        }
    }

    // Oriented throughout or not at all, every edge from an outgoing slot to an incoming one.
    pub fn verify_ori(&self) {
        let n_ori = self.nodes.iter().filter(|x| x.is_oriented()).count();
        if n_ori == 0 {
            return;
        }
        assert_eq!(n_ori, self.n_nodes(), "some nodes are oriented and some are not");

        self.nodes.iter().flat_map(|x| {
            let (p, q) = x.incoming().unwrap();
            Slot::ALL.map(move |s| (x.edge(s), s == p || s == q))
        }).into_group_map().into_iter().for_each(|(e, ins)|
            assert!(
                matches!(ins[..], [a, b] if a != b),
                "edge {e} does not run from an outgoing slot to an incoming one"
            )
        );
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

        self.traverse_comps(|c, i, s| {
            if c == comps.len() {
                comps.push(vec![]);
            }
            comps[c].push(self.node(i).edge(s));
        });

        let mut result: Vec<Path> = comps.into_iter().map(Path::circ).collect();
        for &e in &self.loops {
            result.push(Path::circ(vec![e]));
        }
        result
    }

    pub fn traverse_comps<F>(&self, mut f: F) where 
    F: FnMut(usize, usize, Slot) { 
        let mut c = 0; // component counter
        let mut remain: HashSet<Edge> = self.nodes.iter().flat_map(|x| x.edges().iter().copied()).collect();

        while !remain.is_empty() {
            // Take minimal edge-id. 
            let e0 = remain.iter().min().cloned().unwrap();

            // Find node & point having edge e0, entering at its head so the walk runs forward.
            let (i0, j0) = if self.is_oriented() {
                self.edge_ends(e0, true).1
            } else {
                self.find_port(|i, s|
                    self.node(i).edge(s) == e0
                ).unwrap()
            };

            self.traverse_from((i0, j0), |i, s| { 
                remain.remove(&self.node(i).edge(s));
                f(c, i, s);
            });

            // Onto next component.
            c += 1;
        }
    }

    pub fn traverse_from<F>(&self, start: (usize, Slot), mut f: F) where
        F: FnMut(usize, Slot)
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

    fn traverse_outer(&self, n_index: usize, slot: Slot) -> (usize, Slot) {
        let e = self.nodes[n_index].edge(slot);
        self.nodes.iter().enumerate().flat_map(|(i, _)|
            Slot::ALL.map(move |s| (i, s))
        ).find(|&(i, s)|
            self.nodes[i].edge(s) == e && (i, s) != (n_index, slot)
        ).expect("Broken data")
    }

    // Re-derive each crossing's orientation by traversing components. `is_incoming(i, j)` tells
    // whether port j of node i is known to receive an incoming strand (PD codes: j == 0). The first
    // claimed port met by a tentative traversal fixes the component's direction; a component claiming
    // no port is undetermined (cf. `unlink2`) and the whole link is left unoriented. A fixed direction
    // contradicting `is_incoming` (an odd PD code) panics. Returns whether the link is now oriented.
    pub(crate) fn reorient<F>(&mut self, is_incoming: F) -> bool
    where F: Fn(usize, Slot) -> bool {
        let mut incoming: Vec<Vec<Slot>> = vec![vec![]; self.n_nodes()];
        let mut remain: HashSet<Edge> = self.nodes.iter().flat_map(|x| x.edges().iter().copied()).collect();
        let mut undetermined = false;

        while !remain.is_empty() {
            // start at a claimed port of an untraversed component, so the direction is correct
            // from the outset. Components claiming no port are undetermined (cf. `unlink2`).
            let Some(start) = self.find_port(|i, s|
                remain.contains(&self.node(i).edge(s)) && is_incoming(i, s)
            ) else {
                undetermined = true;
                break;
            };

            self.traverse_from(start, |i, s| {
                remain.remove(&self.node(i).edge(s));
                let out = s.shift(2);
                assert!(
                    is_incoming(i, s) || !is_incoming(i, out),
                    "inconsistent orientation: the strand through node {i} exits at slot {out}, which is claimed incoming"
                );
                incoming[i].push(s);
            });
        }

        // a node is oriented by its two incoming slots, provided they lie on different strands;
        // if any node fails that, or some component is undetermined, the whole link is unoriented.
        let oris = Iterator::zip(self.nodes.iter(), incoming.iter()).map(|(n, slots)|
            match slots[..] {
                [p, q] => Node::orientable(n.node_type(), p, q).then_some((p, q)),
                _ => None,
            }
        ).collect_vec();
        let coherent = !undetermined && oris.iter().all(Option::is_some);

        self.nodes.iter_mut().zip(oris).for_each(|(n, o)| 
            n.set_incoming(if coherent { o } else { None })
        );

        coherent
    }

    pub fn unoriented(&self) -> Self {
        if !self.is_oriented() {
            return self.clone();
        }
        let mut l = self.clone();
        l.nodes.iter_mut().for_each(|n|
            n.set_incoming(None)
        );
        l
    }

    // Renumber the edges base..base+n in the order they are met traversing from `start_edge`
    // (a knot's one traversal covers every edge), keeping the diagram. base_pt becomes `base`
    // (base = 1 gives the usual 1-based numbering of knot theory).
    pub fn reindexed(&self, start_edge: Edge, base: Edge) -> Link {
        assert!(self.is_knot() && self.loops.is_empty(), "reindexed expects a knot");
        assert!(self.is_oriented(), "reindexed needs an orientation to traverse in");

        // note: `traverse_from` runs forward from an in-port, backward from an out-port.
        let (_, start) = self.edge_ends(start_edge, true);

        let mut map: HashMap<Edge, Edge> = HashMap::new();
        let mut next = base;
        self.traverse_from(start, |i, s| {
            map.entry(self.node(i).edge(s)).or_insert_with(|| {
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

    // The canonical relabelling: the least `reindexed(e, 1)` over all start edges. Rebuilt from that
    // PD code so the node order is canonical too, making `==` equality up to relabelling.
    pub fn reindexed_canon(&self) -> Link {
        let pd = self.edges().into_iter()
            .map(|e| self.reindexed(e, 1).pd_code())
            .min()
            .expect("a knot has at least one edge");
        Self::from_pd_code(pd)
    }

    // The two (node, slot) ends of edge `e`. When `directed`, they are ordered as (tail, head)
    // along the orientation — the strand exits at the tail and enters at the head (cf.
    // `Node::incoming`); otherwise the order carries no meaning.
    pub(crate) fn edge_ends(&self, e: Edge, directed: bool) -> ((usize, Slot), (usize, Slot)) {
        assert!(!directed || self.is_oriented(), "directed edge_ends requires an oriented link");

        let (x, y) = self.nodes().enumerate().flat_map(|(i, n)|
            Slot::ALL.into_iter().filter(move |&s| n.edge(s) == e).map(move |s| (i, s))
        ).collect_tuple().unwrap_or_else(||
            panic!("edge {e} must appear exactly twice")
        );
        if !directed {
            return (x, y);
        }

        let is_in = |(i, s): (usize, Slot)| {
            let (p, q) = self.node(i).incoming().expect("directed edge_ends requires an oriented link");
            s == p || s == q
        };
        debug_assert!(is_in(x) != is_in(y), "edge {e} must have one head and one tail");
        if is_in(x) { (y, x) } else { (x, y) }
    }

    fn find_port(&self, f: impl Fn(usize, Slot) -> bool) -> Option<(usize, Slot)> {
        (0..self.n_nodes()).flat_map(|i|
            Slot::ALL.map(move |s| (i, s))
        ).find(|&(i, s)| f(i, s))
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
    #[should_panic(expected = "(edge, count) = [(4, 1), (5, 3)]")]
    fn new_names_the_miscounted_edges() {
        // the trefoil's symmetric PD with edge 4 mistyped as 5
        let _ = Link::from_pd_code([[1, 5, 2, 4], [3, 1, 5, 6], [5, 3, 6, 2]]);
    }

    #[test]
    #[should_panic(expected = "does not run from an outgoing slot")]
    fn new_rejects_disagreeing_orientation() {
        use crate::NodeType::XL;
        // both nodes take edge 1 as incoming, so it would enter at both of its ends.
        let a = Node::new(XL, Some((Slot::SW, Slot::SE)), [1, 2, 3, 4]);
        let b = Node::new(XL, Some((Slot::NE, Slot::NW)), [3, 4, 1, 2]);
        let _ = Link::from_nodes([a, b]);
    }

    #[test]
    #[should_panic(expected = "some nodes are oriented and some are not")]
    fn new_rejects_partial_orientation() {
        use crate::NodeType::XL;
        let a = Node::new(XL, Some((Slot::SW, Slot::SE)), [1, 2, 3, 4]);
        let b = Node::new(XL, None, [3, 4, 1, 2]);
        let _ = Link::from_nodes([a, b]);
    }

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

        let s = |i: usize| Slot::from(i);
        assert_eq!(l.traverse_outer(0, s(0)), (0, s(1)));
        assert_eq!(l.traverse_outer(0, s(1)), (0, s(0)));
        assert_eq!(l.traverse_outer(0, s(2)), (0, s(3)));
        assert_eq!(l.traverse_outer(0, s(3)), (0, s(2)));
    }

    #[test]
    fn link_traverse() {
        let traverse = |l: &Link, start: (usize, Slot)| { 
            let mut queue = vec![];
            l.traverse_from(start, |i, s| queue.push((i, s.index())));
            queue
        };

        let l = Link::test_data("unknot_l_twist");
        let path = traverse(&l, (0, Slot::SW));
        
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
        assert_eq!(l.writhe(), 3);
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

    #[test]
    fn reindexed_numbers_along_the_orientation() {
        // Renumbering depends only on (diagram, start edge), so the least code is a canonical form.
        let canon = |k: &Link| k.reindexed_canon();

        for l in [Link::test_data("3_1"), Link::test_data("6_1"), Link::pretzel(1, 3, 5)] {
            let c = canon(&l);
            for e in l.edges() {
                assert_eq!(canon(&l.reindexed(e, 1)), c, "relabelling from edge {e} changed the canonical form");
            }
        }
    }
}