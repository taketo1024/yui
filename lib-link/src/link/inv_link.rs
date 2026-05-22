use std::collections::HashMap;

use delegate::delegate;
use itertools::Itertools;
use num_integer::Integer;
use crate::{Node, Edge, Link, Path, State, XCode};

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
        let e_map: HashMap<Edge, Edge> = e_map.into_iter().collect();

        // Validate: domain = link edges, image ⊆ link edges, e_map is involutive.
        let link_edges = inner.edges();
        assert_eq!(
            e_map.len(), link_edges.len(),
            "e_map must have one entry per link edge (got {}, expected {})",
            e_map.len(), link_edges.len()
        );
        for &e in &link_edges {
            let f = *e_map.get(&e).unwrap_or_else(|| panic!("e_map missing edge {e}"));
            let g = *e_map.get(&f).unwrap_or_else(|| panic!("e_map({e}) = {f} is not a link edge"));
            assert_eq!(g, e, "e_map is not an involution: {e} ↦ {f} ↦ {g}");
        }

        let mut x_map = HashMap::new();

        for x in inner.nodes() {
            let edges = x.edges().map(|e| e_map.get(&e).unwrap());
            let find = inner.nodes().find_position(|y|
                edges.iter().all(|e| y.edges().contains(e))
            );

            assert!(find.is_some(), "no match for x: {x} -> {edges:?}");

            let j = find.unwrap().0;
            let y = inner.node(j);

            assert_eq!(x.ntype(), y.ntype()); 

            x_map.insert(x.clone(), y.clone());
            if x != y {
                x_map.insert(y.clone(), x.clone());
            }
        }

        assert_eq!(x_map.len(), inner.n_nodes());

        Self { inner, e_map, x_map }
    }

    pub fn from_symmetric_pd_code<I1>(pd_code: I1) -> Self
    where I1: IntoIterator<Item = XCode> { 
        let code = pd_code.into_iter().collect_vec();
        let l = Link::from_pd_code(code);

        assert!(l.is_knot(), "currently only supports strongly invertible knots");

        let n = l.n_edges();

        assert!(n.is_even(), "number of edges must be even.");
        assert_eq!(l.edges().first(), Some(&1), "edge must start from index 1.");
        assert_eq!(l.edges().last(), Some(&n), "edges must have sequential indexing.");

        let e_map = l.edges().into_iter()
            .map(|e| (e, (n + 1 - e) % n + 1));
        
        Self::new(l, e_map).with_base_pt(1)
    }

    pub fn inner(&self) -> &Link { 
        &self.inner
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
            pub fn min_edge(&self) -> Option<Edge>;
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

    pub fn mirror(&self) -> Self {
        Self {
            inner: self.inner.mirror(),
            e_map: self.e_map.clone(),
            x_map: self.x_map.iter().map(|(x, y)|
                (x.mirror(), y.mirror())
            ).collect(),
        }
    }
}

impl InvLink {
    pub fn load(name: &str) -> Result<InvLink, Box<dyn std::error::Error>> {
        let json = yui_core::util::data_dir::load_json("inv_link", name)?;
        let data: Vec<XCode> = serde_json::from_str(&json)?;
        Ok(InvLink::from_symmetric_pd_code(data))
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

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
    fn from_sinv() { 
        let l = InvLink::test_data("3_1");

        assert_eq!(l.inv_edge(1), 1);
        assert_eq!(l.inv_edge(2), 6);
        assert_eq!(l.inv_edge(3), 5);
        assert_eq!(l.inv_edge(4), 4);
        assert_eq!(l.inv_edge(5), 3);
        assert_eq!(l.inv_edge(6), 2);
    }
}