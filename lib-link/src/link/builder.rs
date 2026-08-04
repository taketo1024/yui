use std::collections::{BTreeMap, HashMap};
use std::fmt::Display;

use petgraph::stable_graph::{StableDiGraph, NodeIndex, EdgeIndex};
use yui_core::algo::UnionFind;

use crate::{Link, Node, NodeType, NodeOri, Edge};

// A port is slot `s` (0..4, CCW) of vertex `v`. For a crossing the slots are
//     3   2
//      \ /          under-strand 0 → 2, over-strand 1 → 3  (PD convention)
//      / \
//     0   1
pub type Port = (NodeIndex, usize);

// Two-phase builder for `Link`: add nodes (vertices weighted by `NodeType`), `connect` their ports
// pairwise (graph edges, weighted by the slot at each end), then `build()` validates and orients.
#[derive(Debug, Default)]
pub struct LinkBuilder {
    graph: StableDiGraph<NodeType, (usize, usize)>,
    loops: usize,
}

impl LinkBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn n_nodes(&self) -> usize {
        self.graph.node_count()
    }

    // Add a node of the given type; returns its vertex. Its ports are (vertex, 0..4) in CCW order.
    pub fn add_node(&mut self, node_type: NodeType) -> NodeIndex {
        self.graph.add_node(node_type)
    }

    pub fn add_crossing(&mut self, node_type: NodeType) -> NodeIndex {
        assert!(matches!(node_type, NodeType::XL | NodeType::XR), "add_crossing expects XL or XR, got {node_type}");
        self.add_node(node_type)
    }

    // Join two ports into one edge, recording the slot at each end.
    pub fn connect(&mut self, a: Port, b: Port) {
        assert!(a.1 < 4 && b.1 < 4, "slots must be 0..4, got {} and {}", a.1, b.1);
        self.graph.add_edge(a.0, b.0, (a.1, b.1));
    }

    // Add a free loop (a closed component with no crossings).
    pub fn add_loop(&mut self) {
        self.loops += 1;
    }

    // A column of `k ≥ 1` crossings of type `tt`, chained bottom→top; returns the column's four
    // corner ports in CCW slot order (SW, SE, NE, NW) — for k = 1 the crossing's own ports.
    pub fn add_v_twist(&mut self, tt: NodeType, k: usize) -> (Port, Port, Port, Port) {
        assert!(k >= 1, "a twist needs at least one crossing");
        let v: Vec<_> = (0..k).map(|_| self.add_crossing(tt)).collect();
        for w in v.windows(2) {
            self.connect((w[0], 3), (w[1], 0));
            self.connect((w[0], 2), (w[1], 1));
        }
        ((v[0], 0), (v[0], 1), (v[k - 1], 2), (v[k - 1], 3))
    }

    // A row of `k ≥ 1` crossings of type `tt`, chained left→right; returns the row's four
    // corner ports in CCW slot order (SW, SE, NE, NW) — for k = 1 the crossing's own ports.
    pub fn add_h_twist(&mut self, tt: NodeType, k: usize) -> (Port, Port, Port, Port) {
        assert!(k >= 1, "a twist needs at least one crossing");
        let h: Vec<_> = (0..k).map(|_| self.add_crossing(tt)).collect();
        for w in h.windows(2) {
            self.connect((w[0], 1), (w[1], 0));
            self.connect((w[0], 2), (w[1], 3));
        }
        ((h[0], 0), (h[k - 1], 1), (h[k - 1], 2), (h[0], 3))
    }

    // Absorb all nodes, edges and free loops of `l`; returns its node-index → builder-vertex map (so
    // port `(verts[i], s)` is node i's slot s), for later rewiring.
    pub fn add_link(&mut self, l: &Link) -> Vec<NodeIndex> {
        let verts: Vec<NodeIndex> = l.nodes().map(|x| self.add_node(x.node_type())).collect();

        // group the two occurrences of each edge, sorted by edge id for deterministic numbering.
        let occ = l.nodes().enumerate().flat_map(|(i, x)| {
            let v = verts[i];
            (0..4).map(move |s| (x.edge(s), (v, s)))
        }).fold(BTreeMap::<Edge, Vec<Port>>::new(), |mut occ, (e, p)| {
            occ.entry(e).or_default().push(p);
            occ
        });

        for ports in occ.values() {
            self.connect(ports[0], ports[1]);
        }
        for _ in l.loops() {
            self.add_loop();
        }
        verts
    }

    // Remove the edge incident to port `p` (a connected port has exactly one).
    pub fn disconnect(&mut self, p: Port) {
        let e = self.graph.edge_indices().find(|&e| {
            let (a, b) = self.ends(e);
            a == p || b == p
        }).unwrap_or_else(||
            panic!("port (node {}, slot {}) has no edge to disconnect", p.0.index(), p.1)
        );
        self.graph.remove_edge(e);
    }

    // The `Edge` id `build()` will give the edge incident to port `p` (build() numbers edges 1.. in
    // `edge_indices()` order). Call after the final wiring so the numbering matches.
    pub fn edge_at(&self, p: Port) -> Option<Edge> {
        self.graph.edge_indices().enumerate().find_map(|(i, e)| {
            let (a, b) = self.ends(e);
            (a == p || b == p).then_some(i as Edge + 1)
        })
    }

    pub fn build(self) -> Result<Link, LinkError> {
        self.build_with(|_, _| true) // any port may be incoming: each direction is a free choice
    }

    // Build with orientation knowledge: `is_incoming(i, s)` tells whether slot `s` of node `i`
    // (in insertion order) receives an incoming strand — see `Link::reorient`.
    pub fn build_with<F>(self, is_incoming: F) -> Result<Link, LinkError>
    where F: Fn(usize, usize) -> bool {
        self.validate()?;

        // number the edges 1.., recording each id at its two ports.
        let edge_at: HashMap<Port, Edge> = self.graph.edge_indices().enumerate().flat_map(|(i, e)| {
            let (a, b) = self.ends(e);
            let id = i as Edge + 1;
            [(a, id), (b, id)]
        }).collect();

        let nodes = self.graph.node_indices().map(|v|
            Node::new(self.graph[v], NodeOri::None, [0, 1, 2, 3].map(|s| edge_at[&(v, s)]))
        );

        let e0 = self.graph.edge_count() as Edge + 1;
        let loops = e0 .. e0 + self.loops as Edge;

        let mut l = Link::new(nodes.collect::<Vec<_>>(), loops);
        l.reorient(is_incoming);
        Ok(l)
    }

    // All build-blocking defects, most specific first: doubly-used / open ports, then genus.
    fn validate(&self) -> Result<(), LinkError> {
        let alpha = self.edge_pairing()?;
        if self.is_planar_with(&alpha) {
            Ok(())
        } else {
            Err(LinkError::NonPlanar)
        }
    }

    // False if any port is open or doubly connected.
    pub fn is_planar(&self) -> bool {
        self.edge_pairing().map(|alpha|
            self.is_planar_with(&alpha)
        ).unwrap_or(false)
    }

    // A planar (genus-0) wiring satisfies V − E + F = 2·#components (Euler, componentwise). The
    // faces are the orbits of φ = rotate ∘ α on the 4V ports, where α pairs the two ends of each
    // edge and rotate steps to the next CCW slot.
    fn is_planar_with(&self, alpha: &HashMap<Port, Port>) -> bool {
        let (v, e) = (self.graph.node_count(), self.graph.edge_count());

        let dart = |(n, s): &Port| 4 * n.index() + s;
        let phi = |p: &Port| {
            let (n, s) = alpha[p];
            (n, (s + 1) % 4)
        };

        // count orbits by union-find: faces are the orbits of φ, components those of α (per node).
        let mut faces = UnionFind::new(4 * v);
        let mut comps = UnionFind::new(v);
        for p in alpha.keys() {
            faces.union(dart(p), dart(&phi(p)));
            comps.union(p.0.index(), alpha[p].0.index());
        }
        let (f, c) = (faces.into_disjoint().len(), comps.into_disjoint().len());

        (v + f) as isize - e as isize == 2 * c as isize
    }

    // α: pairs each port with the port at the other end of its edge; errs unless every port is
    // connected exactly once.
    fn edge_pairing(&self) -> Result<HashMap<Port, Port>, LinkError> {
        let mut alpha = HashMap::new();
        for e in self.graph.edge_indices() {
            let (a, b) = self.ends(e);
            for (p, q) in [(a, b), (b, a)] {
                if alpha.insert(p, q).is_some() {
                    return Err(LinkError::DuplicatePort(p));
                }
            }
        }

        let open = self.graph.node_indices().flat_map(|v|
            (0..4).map(move |s| (v, s))
        ).find(|p|
            !alpha.contains_key(p)
        );
        match open {
            Some(p) => Err(LinkError::UnconnectedPort(p)),
            None => Ok(alpha),
        }
    }

    // The two ports of a graph edge.
    fn ends(&self, e: EdgeIndex) -> (Port, Port) {
        let (s, t) = self.graph.edge_endpoints(e).unwrap();
        let &(fs, ts) = self.graph.edge_weight(e).unwrap();
        ((s, fs), (t, ts))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LinkError {
    // A port was never connected.
    UnconnectedPort(Port),
    // A port was claimed by two edges (e.g. a self-edge with both ends on the same slot).
    DuplicatePort(Port),
    // The wiring has positive genus (fails Euler's formula for its ribbon structure).
    NonPlanar,
}

impl Display for LinkError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LinkError::UnconnectedPort((v, s)) =>
                write!(f, "port (node {}, slot {s}) is not connected", v.index()),
            LinkError::DuplicatePort((v, s)) =>
                write!(f, "port (node {}, slot {s}) is connected more than once", v.index()),
            LinkError::NonPlanar =>
                write!(f, "the wiring is non-planar (genus > 0)"),
        }
    }
}

impl std::error::Error for LinkError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::misc::jones_polynomial;

    fn rebuild(l: &Link) -> Link {
        let mut b = LinkBuilder::new();
        b.add_link(l);
        b.build().unwrap()
    }

    #[test]
    fn round_trip_preserves_knot() {
        for name in ["3_1", "4_1", "5_2", "6_2"] {
            let l = Link::test_data(name);
            let r = rebuild(&l);
            assert_eq!(r.n_comps(), l.n_comps());
            assert!(r.is_oriented());
            assert_eq!(jones_polynomial(&r), jones_polynomial(&l), "round-trip changed {name}");
        }
    }

    #[test]
    fn round_trip_orients_unlink2() {
        // unlink2 loads unoriented (its over-component has no under-anchor in the PD code), but the
        // builder orients freely; it also exercises reorient's coverage of never-under components.
        let l = Link::test_data("unlink2");
        assert!(!l.is_oriented());

        let r = rebuild(&l);
        assert_eq!(r.n_comps(), 2);
        assert!(r.is_oriented());
        assert_eq!(jones_polynomial(&r), jones_polynomial(&Link::unlink(2)));
    }

    #[test]
    fn is_planar_detects_genus() {
        // a kink (adjacent self-connections) is planar; connecting opposite ports is genus 1.
        let mut planar = LinkBuilder::new();
        let x = planar.add_crossing(NodeType::XR);
        planar.connect((x, 0), (x, 1));
        planar.connect((x, 2), (x, 3));
        assert!(planar.is_planar());

        let mut torus = LinkBuilder::new();
        let y = torus.add_crossing(NodeType::XR);
        torus.connect((y, 0), (y, 2));
        torus.connect((y, 1), (y, 3));
        assert!(!torus.is_planar());
        assert!(matches!(torus.build(), Err(LinkError::NonPlanar)));
    }

    #[test]
    fn is_planar_accepts_real_knots() {
        // every from_pd_code diagram is genuinely planar, so is_planar must accept all of them.
        for name in ["3_1", "4_1", "5_1", "5_2", "6_1", "6_2", "6_3", "7_1", "7_2"] {
            let mut b = LinkBuilder::new();
            b.add_link(&Link::test_data(name));
            assert!(b.is_planar(), "{name} is planar but is_planar returned false");
        }
    }

    #[test]
    fn build_with_free_loop() {
        let mut b = LinkBuilder::new();
        b.add_loop();
        let l = b.build().unwrap();
        assert_eq!(l.n_comps(), 1);
        assert_eq!(l.n_crossings(), 0);
    }

    #[test]
    fn unconnected_port_errors() {
        let mut b = LinkBuilder::new();
        b.add_crossing(NodeType::XR);
        assert!(matches!(b.build(), Err(LinkError::UnconnectedPort(_))));
    }
}
