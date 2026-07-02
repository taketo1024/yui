use core::panic;
use std::collections::HashSet;
use std::fmt::Display;
use itertools::Itertools;
use yui_core::{CloneAnd, Sign};
use yui_core::bitseq::Bit;

use petgraph::Graph;

use super::{Node, Path};

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
pub type PDCodeX = [Edge; 4];

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

    // Planer Diagram code, represented by sequence of crossings of the form:
    //
    //     d   c
    //      \ /
    //       \     = [a, b, c, d]
    //      / \
    //     a   b
    //
    // The lower edge is always oriented a -> c.
    // see: http://katlas.math.toronto.edu/wiki/Planar_Diagrams

    pub fn from_pd_code<I>(pd_code: I) -> Self
    where I: IntoIterator<Item = PDCodeX> {
        let nodes = pd_code.into_iter().map(Node::from_pd_code).collect_vec();
        let mut l = Self::from_nodes(nodes); // unoriented
        l.reorient(|_, j| j == 0); // PD convention: the under-strand enters at pos 0.
        l
    }

    // Re-derive each crossing's orientation by traversing components. `is_incoming(i, j)` tells
    // whether port j of node i is known to receive an incoming strand (PD codes: j == 0). The first
    // claimed port met by a tentative traversal fixes the component's direction; a component claiming
    // no port is undetermined (cf. `unlink2`) and the whole link is left unoriented. A fixed direction
    // contradicting `is_incoming` (an odd PD code) panics. Returns whether the link is now oriented.
    pub(crate) fn reorient<F>(&mut self, is_incoming: F) -> bool
    where F: Fn(usize, usize) -> bool {
        use crate::NodeOri::{Up, Down, Left, Right, None};

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

        // a crossing's two incoming ports are an adjacent pair, which fixes the orientation;
        // if any node is incoherent, or some component is undetermined, the whole link is unoriented.
        let oris = incoming.iter().map(|ports| match ports[..] {
            [0, 1] | [1, 0] => Up,
            [1, 2] | [2, 1] => Left,
            [2, 3] | [3, 2] => Down,
            [3, 0] | [0, 3] => Right,
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

    pub fn load(name: &str) -> Result<Link, Box<dyn std::error::Error>> {
        let json = yui_core::util::data_dir::load_json("links", name)?;
        let data: Vec<PDCodeX> = serde_json::from_str(&json)?;
        Ok(Link::from_pd_code(data))
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

    pub fn mirror(&self) -> Self {
        let mut l = Self::new(
            self.nodes().map(|x| x.mirror()),
            self.loops.iter().copied(),
        );
        l.base_pt = self.base_pt;
        l
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

    pub fn node_mut(&mut self, i: usize) -> &mut Node { 
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

    pub fn cc_at(&self, i: usize) -> Self {
        assert!(self.node(i).is_crossing());
        self.clone_and(|l| 
            *l.node_mut(i) = l.node(i).mirror()
        )
    }

    pub fn resolve_at(&self, i: usize, r: Bit) -> Self {
        assert!(self.node(i).is_crossing());
        self.clone_and(|l| 
            *l.node_mut(i) = l.node(i).resolve(r)
        )
    }

    pub fn resolve_by(&self, s: &State) -> Self {
        assert!(s.len() == self.n_crossings());

        let n = self.nodes.len();
        let itr = (0..n).filter(|&i| self.node(i).is_crossing());

        self.clone_and(|l| {
            for (i, r) in Iterator::zip(itr, s.iter()) {
                *l.node_mut(i) = self.node(i).resolve(r); 
            }
        })
    }

    pub fn seifert_state(&self) -> State {
        assert!(self.is_oriented());

        let seq = self.crossings().map(|x|
            match x.sign() {
                Some(Sign::Pos) => 0,
                Some(Sign::Neg) => 1,
                None => panic!("Impossible.")
            }
        ); 
        State::from_iter(seq)
    }

    pub fn seifert_circles(&self) -> Vec<Path> {
        self.resolve_by(&self.seifert_state()).comps()
    }

    pub fn seifert_graph(&self) -> Graph<Path, usize> {
        assert!(self.is_oriented());

        use crate::NodeType;
        type G = Graph<Path, usize>;

        let s0 = self.seifert_state();
        let l0 = self.resolve_by(&s0);
        let mut graph = Graph::new();

        // Vertices = Seifert circles (and free loops contribute their own circles).
        for c in l0.comps() {
            graph.add_node(c);
        }

        let find_node = |graph: &G, e| {
            graph.node_indices().find(|&i|
                graph[i].contains(e)
            )
        };

        // Edges = one per original crossing (now resolved into a V/H smoothing).
        // Free loops have no nodes, so they remain isolated vertices.
        for (i, x) in l0.nodes().enumerate() {
            let (e1, e2) = if x.node_type() == NodeType::V {
                (x.edge(0), x.edge(1))
            } else {
                (x.edge(0), x.edge(2))
            };
            let n1 = find_node(&graph, e1).unwrap();
            let n2 = find_node(&graph, e2).unwrap();
            graph.add_edge(n1, n2, i);
        }

        graph
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
}

impl Display for Link {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "L[{}]", self.nodes.iter().map(|x| x.to_string()).join(", "))
    }
}

#[cfg(test)]
mod tests { 
    use crate::NodeType::{XL, XR};

    use super::*;

    #[test]
    fn link_init() { 
        let l = Link::from_nodes(vec![]);
        assert_eq!(l.nodes.len(), 0);
    }

    #[test]
    fn link_from_pd_code() { 
        let l = Link::test_data("unknot_l_twist");
        assert_eq!(l.nodes.len(), 1);
        assert_eq!(l.node(0).node_type(), XL);
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
        assert_eq!(comps, vec![ Path::new(vec![1, 2], true)]);
    }

    #[test]
    fn link_mirror() { 
        let l = Link::test_data("unknot_l_twist");
        assert_eq!(l.node(0).node_type(), XL);

        let l = l.mirror();
        assert_eq!(l.node(0).node_type(), XR);
    }

    #[test]
    fn link_resolve() {
        let s = State::from([0, 0, 0]);
        let l = Link::test_data("3_1").resolve_by(&s);

        let comps = l.comps();
        assert_eq!(comps.len(), 3);
        assert!(comps.iter().all(|c| c.is_circle()));

        let s = State::from([1, 1, 1]);
        let l = Link::test_data("3_1").resolve_by(&s);

        let comps = l.comps();
        assert_eq!(comps.len(), 2);
        assert!(comps.iter().all(|c| c.is_circle()));
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
    fn mirror_preserves_loops() {
        let l = Link::unlink(1).mirror();
        assert_eq!(l.n_loops(), 1);
        assert_eq!(l.loops(), &[1]);
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
    fn crossing_change() {
        use crate::link::node::NodeOri;

        let l = Link::test_data("3_1");
        let l2 = l.cc_at(1);

        assert_eq!(l.node(1),  &Node::new(XL, NodeOri::Up, [3,6,4,1]));
        assert_eq!(l2.node(1), &Node::new(XR, NodeOri::Up, [3,6,4,1]));
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
    fn mirror_preserves_base_pt() {
        let l = Link::test_data("3_1").with_base_pt(3).mirror();
        assert_eq!(l.base_pt(), Some(3));
    }

    #[test]
    #[should_panic]
    fn with_base_pt_invalid_panics() {
        // Edge 99 is not in the trefoil's edge set.
        let _ = Link::test_data("3_1").with_base_pt(99);
    }

    #[test]
    fn seifert_graph_trefoil() {
        let l = Link::test_data("3_1");
        let g = l.seifert_graph();
        assert_eq!(g.node_count(), 2);
        assert_eq!(g.edge_count(), 3);
    }

    #[test]
    fn seifert_graph_unlink() {
        // Free loops contribute isolated vertices and no edges.
        let l = Link::unlink(3);
        let g = l.seifert_graph();
        assert_eq!(g.node_count(), 3);
        assert_eq!(g.edge_count(), 0);
    }
}