use std::collections::HashMap;

use delegate::delegate;
use itertools::Itertools;
use crate::{Node, Edge, Link, Path, Slot, State, PDCodeX};

// Involutive link
#[derive(Debug, Clone)]
pub struct InvLink {
    inner: Link,
    e_map: HashMap<Edge, Edge>,
    x_map: HashMap<Node, Node>
}

impl InvLink {
    pub fn new<I>(inner: Link, e_map: I) -> InvLink
    where I: IntoIterator<Item = (Edge, Edge)> {
        assert!(inner.loops().is_empty(), "free loops are not supported");
        let e_map: HashMap<Edge, Edge> = e_map.into_iter().collect();

        // `e_map` must be an involution of the whole edge set.
        let link_edges = inner.edges();
        assert_eq!(e_map.len(), link_edges.len(), "e_map must cover every edge exactly once");
        for &e in &link_edges {
            let f = e_map[&e];
            assert_eq!(e_map[&f], e, "e_map is not involutive at edge {e}");
        }

        // ... and must carry each node to a node of the same type.
        let mut x_map = HashMap::new();
        for x in inner.nodes() {
            let edges = x.edges().map(|e| e_map[&e]);
            let y = inner.nodes()
                .find(|y| edges.iter().all(|e| y.edges().contains(e)))
                .unwrap_or_else(|| panic!("e_map sends node {x} to no node of the link"));
            assert_eq!(x.node_type(), y.node_type(), "e_map changes the type of node {x}");

            x_map.insert(x.clone(), y.clone());
            if x != y {
                x_map.insert(y.clone(), x.clone());
            }
        }
        assert_eq!(x_map.len(), inner.n_nodes(), "e_map is not a bijection on nodes");

        Self { inner, e_map, x_map }
    }

    pub fn from_symmetric_pd_code<I1>(pd_code: I1) -> Self
    where I1: IntoIterator<Item = PDCodeX> { 
        // the base point defaults to the least edge, which the symmetric convention puts on the axis.
        Self::si_knot_from(Link::from_pd_code(pd_code))
    }

    // A strongly invertible knot, given a diagram based on its axis. τ reverses the traversal, so
    // walking both ways from the base point pairs each edge with its image.
    pub(super) fn si_knot_from(inner: Link) -> InvLink {
        assert!(inner.is_knot(), "a strongly invertible knot must be a knot");
        let base = inner.base_pt().expect("the diagram needs a base point on the axis");
        let comps = inner.comps();
        let seq = comps[0].edges();
        let n = seq.len();
        let k = seq.iter().position(|&e| e == base).expect("the base point is not on a strand");

        let e_map = (0..n).map(|i| (seq[(k + i) % n], seq[(k + n - i) % n])).collect_vec();
        let l = Self::new(inner.clone(), e_map);
        assert!(l.is_strongly_invertible(), "the based diagram is not a strongly invertible knot");
        l
    }

    pub fn inner(&self) -> &Link {
        &self.inner
    }

    /// The symmetric PD code — the inner diagram's PD, with edges numbered so the strong inversion is
    /// the standard `e ↦ (n+1-e)%n+1`. Feeding this to `from_symmetric_pd_code` (or the `ykh` CLI)
    /// reconstructs the same `InvLink`. The builders (`from_symmetric_pd_code`, `conn_sum`,
    /// `whitehead_double`) already reindex to this involution, so an arbitrary `e_map` is rejected.
    pub fn pd_code(&self) -> Vec<PDCodeX> {
        let n = self.inner.n_edges() as Edge;
        assert!(
            self.inner.edges().into_iter().all(|e| self.inv_edge(e) == (n + 1 - e) % n + 1),
            "pd_code requires the standard involution `e ↦ (n+1-e)%n+1`; reindex the InvLink first"
        );
        self.inner.pd_code()
    }

    // delegate methods from Link

    delegate! {
        to self.inner {
            pub fn is_empty(&self) -> bool;
            pub fn is_knot(&self) -> bool;
            pub fn is_oriented(&self) -> bool;
            pub fn writhe(&self) -> i32;
            pub fn n_nodes(&self) -> usize;
            pub fn nodes(&self) -> impl Iterator<Item = &Node>;
            pub fn node(&self, i: usize) -> &Node;
            pub fn crossings(&self) -> impl Iterator<Item = &Node>;
            pub fn n_crossings(&self) -> usize;
            pub fn n_signed_crossings(&self) -> (usize, usize);
            pub fn n_edges(&self) -> usize;
            pub fn edges(&self) -> Vec<Edge>;
            pub fn n_comps(&self) -> usize;
            pub fn comps(&self) -> Vec<Path>;
            pub fn seifert_state(&self) -> State;
            pub fn seifert_circles(&self) -> Vec<Path>;
            pub fn base_pt(&self) -> Option<Edge>;
        }
    }

    pub fn with_base_pt(mut self, e: Edge) -> Self {
        assert_eq!(self.inv_edge(e), e, "base_pt {e} must be on-axis (fixed by involution)");
        self.inner = self.inner.with_base_pt(e);
        self
    }

    pub fn inv_edge(&self, e: Edge) -> Edge { 
        self.e_map.get(&e).cloned().unwrap()
    }

    pub fn inv_node(&self, x: &Node) -> &Node {
        self.x_map.get(x).unwrap()
    }

    pub fn on_axis_edges(&self) -> Vec<Edge> {
        self.edges().into_iter().filter(|&e| self.inv_edge(e) == e).collect()
    }

    // A strong inversion reverses the orientation.
    pub fn is_strongly_invertible(&self) -> bool {
        self.is_oriented()
            && self.inner.nodes().all(|x| self.preserves_dir_at(x) == Some(false))
    }

    // A 2-periodic link carries an involution preserving the orientation.
    pub fn is_2periodic(&self) -> bool {
        self.is_oriented()
            && self.inner.nodes().all(|x| self.preserves_dir_at(x) == Some(true))
    }

    // Whether τ keeps the strands running the same way at `x`. The rotation reverses the cyclic
    // order of a node's slots, sending both incoming slots to incoming ones, or both to outgoing.
    fn preserves_dir_at(&self, x: &Node) -> Option<bool> {
        let y = self.inv_node(x);
        let (p, q) = x.incoming()?;
        let k = Slot::ALL.into_iter().find(|k|
            Slot::ALL.into_iter().all(|s| y.edge(k.shift(4 - s.index())) == self.inv_edge(x.edge(s)))
        )?;

        let is_in = |s: Slot| {
            let img = k.shift(4 - s.index());
            y.incoming().is_some_and(|(a, b)| img == a || img == b)
        };
        match (is_in(p), is_in(q)) {
            (true, true) => Some(true),
            (false, false) => Some(false),
            _ => None,
        }
    }

    pub fn mirror(&self) -> Self {
        Self {
            inner: self.inner.mirror(),
            e_map: self.e_map.clone(),
            x_map: self.x_map.iter().map(|(x, y)|
                (x.mirror(), y.mirror())
            ).collect(),
        }
    }

    // Equivariant connected sum: splice self's other on-axis edge to other's base point, so self's
    // base point survives as the sum's.
    pub fn conn_sum(&self, other: &InvLink) -> InvLink {
        let base = self.base_pt().expect("self needs a base point");
        let self_e = self.on_axis_edges().into_iter()
            .find(|&e| e != base)
            .expect("self needs a second on-axis edge");
        let other_e = other.base_pt().expect("other needs a base point");
        self.conn_sum_at(other, self_e, other_e)
    }

    // Equivariant connected sum: splice along on-axis edges (`inv_edge(e) == e`) of each summand,
    // then recover the combined strong inversion by reindexing to the standard involution.
    pub fn conn_sum_at(&self, other: &InvLink, self_e: Edge, other_e: Edge) -> InvLink {
        assert!(self.is_knot() && other.is_knot(), "connected sum requires knots");
        assert!(self.is_strongly_invertible() && other.is_strongly_invertible(),
            "connected sum requires strongly invertible knots");
        assert_ne!(Some(self_e), self.base_pt(), "the splice must not consume self's base point");
        assert_eq!(self.inv_edge(self_e), self_e, "self_e {self_e} must be on-axis");
        assert_eq!(other.inv_edge(other_e), other_e, "other_e {other_e} must be on-axis");

        // `Link::conn_sum_at` carries self's base point through the splice, and it is on-axis.
        let inner = self.inner.conn_sum_at(other.inner(), self_e, other_e);
        Self::si_knot_from(inner)
    }
}

impl InvLink {
    pub fn load(name: &str) -> Result<InvLink, Box<dyn std::error::Error>> {
        let json = yui_core::util::data_dir::load_json("inv_link", name)?;
        let data: Vec<PDCodeX> = serde_json::from_str(&json)?;
        Ok(InvLink::from_symmetric_pd_code(data))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Slot;
    use crate::misc::det;

    #[test]
    fn reindexed_keeps_strong_inversion() {
        // reindex the symmetric trefoil from edge 1; the standard e↦(n+1-e)%n+1 must stay a valid τ.
        let il = InvLink::test_data("3_1");
        let r = il.inner().reindexed(1, 1);
        assert_eq!(r.edges(), (1..=6).collect::<Vec<Edge>>());

        let n = r.n_edges() as Edge;
        let e_map: Vec<_> = r.edges().into_iter().map(|e| (e, (n + 1 - e) % n + 1)).collect();
        InvLink::new(r, e_map);  // panics if the involution is invalid
    }

    #[test]
    fn pd_code_roundtrip() {
        // 5_1 as a symmetric PD (standard involution by construction); emit + reparse must recover
        // the same diagram and the same strong inversion.
        let l = InvLink::from_symmetric_pd_code([[1,7,2,6],[3,9,4,8],[5,1,6,10],[7,3,8,2],[9,5,10,4]]);
        let l2 = InvLink::from_symmetric_pd_code(l.pd_code());
        assert_eq!(l.pd_code(), l2.pd_code());
        for e in l.inner().edges() {
            assert_eq!(l.inv_edge(e), l2.inv_edge(e), "involution differs at edge {e}");
        }
    }

    #[test]
    fn load_3_1() {
        let l = InvLink::test_data("3_1");
        assert_eq!(l.n_crossings(), 3);
    }

    #[test]
    fn load_4_1() { 
        let l = InvLink::test_data("4_1");
        assert_eq!(l.n_crossings(), 4);
    }
    
    #[test]
    fn inv_edge() { 
        let l = InvLink::test_data("3_1");

        assert_eq!(l.inv_edge(1), 1);
        assert_eq!(l.inv_edge(2), 6);
        assert_eq!(l.inv_edge(3), 5);
        assert_eq!(l.inv_edge(4), 4);
        assert_eq!(l.inv_edge(5), 3);
        assert_eq!(l.inv_edge(6), 2);
    }
    
    #[test]
    fn inv_node() {
        let l = InvLink::test_data("3_1");
        let nodes = l.inner.nodes().collect_vec();

        assert_eq!(l.inv_node(&nodes[0]), nodes[1]);
        assert_eq!(l.inv_node(&nodes[1]), nodes[0]);
        assert_eq!(l.inv_node(&nodes[2]), nodes[2]);
    }

    #[test]
    fn mirror_keeps_the_involution() {
        let k = InvLink::test_data("3_1");
        let m = k.mirror();
        assert_eq!(m.n_crossings(), k.n_crossings());
        for e in k.edges() {
            assert_eq!(m.inv_edge(e), k.inv_edge(e), "involution differs at edge {e}");
        }
        for (x, y) in k.inner().nodes().zip(m.inner().nodes()) {
            assert_eq!(m.inv_node(y), &k.inv_node(x).mirror(), "τ differs at node {x:?}");
        }
    }

    #[test]
    fn twist_unknot() {
        // the 1-crossing unknot: the identity is the π-rotation through the crossing, a strong
        // inversion; swapping the two edges is a rotation missing the knot, so 2-periodic.
        let l = Link::from_pd_code([[0, 1, 1, 0]]);

        let si = InvLink::new(l.clone(), [(0, 0), (1, 1)]);
        assert_eq!(si.on_axis_edges(), vec![0, 1]);
        assert!(si.is_strongly_invertible());
        assert!(!si.is_2periodic());

        let per = InvLink::new(l, [(0, 1), (1, 0)]);
        assert!(per.on_axis_edges().is_empty());
        assert!(!per.is_strongly_invertible());
        assert!(per.is_2periodic());
    }

    #[test]
    fn involution_is_a_strong_inversion() {
        // `new` only checks that the edge map is an involution; these are the conditions making it
        // a π-rotation about an axis in the plane.
        fn check(name: &str, l: &InvLink) {
            let inner = l.inner();
            assert_eq!(l.on_axis_edges().len(), 2, "{name}: the axis must meet the knot twice");
            assert!(l.is_strongly_invertible(), "{name}: τ must reverse the orientation");

            for x in inner.nodes() {
                let y = l.inv_node(x);
                assert_eq!(y.node_type(), x.node_type(), "{name}: τ changed a crossing type");
                assert_eq!(l.inv_node(y), x, "{name}: τ is not an involution on nodes");

                // the rotation reverses the cyclic order of a crossing's slots: s ↦ (k - s) mod 4, k odd.
                let k = (0..4).find(|&k|
                    Slot::ALL.iter().all(|&s|
                        y.edge(Slot::from((k + 4 - s.index()) % 4)) == l.inv_edge(x.edge(s))
                    )
                );
                assert!(matches!(k, Some(1) | Some(3)), "{name}: τ does not reverse the slot order at {x}");
            }
        }

        for name in ["3_1", "4_1", "6_3"] {
            check(name, &InvLink::test_data(name));
        }
        check("sym_pretzel(-3,3,-3)", &InvLink::sym_pretzel(-3, 3, -3));
        check("sym_wh+(3_1)", &InvLink::whitehead_double(&InvLink::test_data("3_1"), true, 0));
    }

    #[test]
    fn conn_sum_is_equivariant() {
        // construction succeeding ⟺ from_standard_reindex found a valid τ on the sum.
        let k1 = InvLink::test_data("3_1");
        let k2 = InvLink::test_data("4_1");
        let cs = k1.conn_sum(&k2);
        assert!(cs.is_knot());
        assert!(cs.is_strongly_invertible());
        let base = cs.base_pt().expect("the sum keeps a base point");
        assert_eq!(cs.inv_edge(base), base, "the sum's base point is on-axis");
        assert_eq!(det(cs.inner()), 3 * 5, "det is multiplicative under conn sum");
    }
}