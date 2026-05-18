use core::panic;
use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use itertools::Itertools;
use yui_core::{hashmap, CloneAnd, Sign};
use yui_core::bitseq::Bit;

use crate::NodeType;
use crate::link::node::NodeOri;

use super::{Node, Path};

pub type Edge = usize;
pub type State = yui_core::bitseq::BitSeq;
pub type XCode = [Edge; 4];

#[derive(Debug, Clone)]
pub struct Link { 
    nodes: Vec<Node>,
    edges: HashSet<Edge>
}

impl Link {
    pub fn from_nodes(nodes: impl IntoIterator<Item = Node>) -> Self { 
        let nodes = nodes.into_iter().collect_vec();
        let edges = nodes.iter().flat_map(|x| x.edges()).cloned().collect();
        let l = Self { nodes, edges };
        l.validate();
        l
    }

    fn validate(&self) { 
        assert_eq!(self.edges.len(), self.nodes.len() * 2, "Invalid data.");
        self.traverse(|_, _, _| ());
    }

    // Planer Diagram code, represented by crossings:
    //
    //     3   2
    //      \ /
    //       \      = (0, 1, 2, 3)
    //      / \
    //     0   1
    //
    // The lower edge has direction 0 -> 2.
    // The crossing is +1 if the upper goes 3 -> 1.
    // see: http://katlas.math.toronto.edu/wiki/Planar_Diagrams

    pub fn from_pd_code<I>(pd_code: I) -> Self
    where I: IntoIterator<Item = XCode> { 
        let nodes = pd_code.into_iter().map(Node::from_pd_code).collect_vec();
        Self::from_nodes(nodes)
    }

    pub fn load(name: &str) -> Result<Link, Box<dyn std::error::Error>> {
        let json = yui_core::util::data_dir::load_json("links", name)?;
        let data: Vec<XCode> = serde_json::from_str(&json)?;
        Ok(Link::from_pd_code(data))
    }

    pub fn empty() -> Link {
        Link { nodes: vec![], edges: HashSet::new() }
    }

    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    pub fn unknot() -> Link {
        let n = Node::new(NodeType::H, NodeOri::None, [1, 2, 2, 1]);
        Link::from_nodes([n])
    }

    pub fn is_knot(&self) -> bool { 
        self.n_comps() == 1
    }

    pub fn writhe(&self) -> i32 { 
        let (p, n) = self.count_signed_crossings();
        (p as i32) - (n as i32)
    }

    pub fn mirror(&self) -> Self {
        Self::from_nodes(self.nodes().map(|x| x.mirror()))
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

    pub fn count_signed_crossings(&self) -> (usize, usize) {
        let signs = self.iter_signed_crossings();
        let pos = signs.iter().filter(|(_, s)| s.is_positive()).count();
        let neg = signs.len() - pos;
        (pos, neg)
    }

    pub fn iter_signed_crossings(&self) -> HashMap<usize, Sign> {
        // TODO replace this. 
        use super::node::NodeType::*;

        let mut result = hashmap!{};

        self.traverse(|_, i, j| {
            match (self.node(i).ntype(), j) { 
                (XR, 1) | (XL, 3) => result.insert(i, Sign::Pos),
                (XR, 3) | (XL, 1) => result.insert(i, Sign::Neg),
                _ => None
            };
        });

        result
    }

    pub fn n_edges(&self) -> usize { 
        self.edges.len()
    }
    
    pub fn edges(&self) -> impl Iterator<Item = &Edge> {
        self.edges.iter()
    }

    pub fn min_edge(&self) -> Option<Edge> { 
        self.nodes.first().map(|x| x.min_edge())
    }

    pub fn n_comps(&self) -> usize { 
        let mut count = 0;
        self.traverse(|c, _, _| 
            if count <= c { count = c + 1 } 
        );
        count
    }

    pub fn comps(&self) -> Vec<Path> {
        let mut comps = vec![];

        self.traverse(|c, i, j| { 
            if c == comps.len() { 
                comps.push(vec![]);
            }

            let e = self.node(i).edge(j);
            comps[c].push(e);
        });

        comps.into_iter().map(|edges| 
            Path::circ(edges)
        ).collect()
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
        let signs = self.iter_signed_crossings();
        let seq = signs.iter().sorted_by_key(|(&i, _)| i).map(|(_, s)| 
            match s { 
                Sign::Pos => 0, 
                Sign::Neg => 1
            }
        ); 
        State::from_iter(seq)
    }

    pub fn seifert_circles(&self) -> Vec<Path> { 
        self.resolve_by(&self.seifert_state()).comps()
    }

    pub fn traverse<F>(&self, mut f: F) where 
    F: FnMut(usize, usize, usize) { 
        let n = self.n_nodes();

        let mut c = 0; // component counter
        let mut remain = self.edges.clone();

        // MEMO: For a link obtained from a PD-code, 
        // the following loop should break after the first iteration: j0 = 0.

        for j0 in [0, 1, 2] { 
            for i0 in 0..n {
                let e0 = self.node(i0).edge(j0);
                if !remain.remove(&e0) { 
                    continue 
                }

                self.traverse_from((i0, j0), |i, j| { 
                    let e = self.node(i).edge(j);
                    remain.remove(&e);

                    f(c, i, j);
                });

                c += 1;
            }

            if remain.is_empty() { 
                break
            }
        }

        assert!(remain.is_empty())
    }

    pub fn traverse_from<F>(&self, start: (usize, usize), mut f:F) where
        F: FnMut(usize, usize)
    {
        let (mut i, mut j) = start;

        f(i, j); // call starting point

        loop {
            let c = self.node(i);
            let k = c.traverse_inner(j);
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
    use yui_core::hashmap;
    use crate::NodeType::{XL, XR};

    use super::*;

    #[test]
    fn link_init() { 
        let l = Link::from_nodes(vec![]);
        assert_eq!(l.nodes.len(), 0);
    }

    #[test]
    fn link_from_pd_code() { 
        let l = Link::unknot_l_twist();
        assert_eq!(l.nodes.len(), 1);
        assert_eq!(l.node(0).ntype(), XL);
    }

    #[test]
    fn link_is_empty() {
        let l = Link::empty();
        assert!(l.is_empty());

        let l = Link::unknot_l_twist();
        assert!(!l.is_empty());
    }

    #[test]
    fn link_crossing_num() {
        let l = Link::empty();
        assert_eq!(l.n_crossings(), 0);

        let l = Link::unknot_l_twist();
        assert_eq!(l.n_crossings(), 1);
        
        let l = Link::test_data("3_1").unwrap();
        assert_eq!(l.n_crossings(), 3);
    }

    #[test]
    fn link_next() {
        let l = Link::unknot_l_twist();

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

        let l = Link::unknot_l_twist();
        let path = traverse(&l, (0, 0));
        
        assert_eq!(path, [(0, 0), (0, 3)]); // loop
    }

    #[test]
    fn link_crossing_signs() {
        let l = Link::unknot_l_twist();
        assert_eq!(l.iter_signed_crossings(), hashmap!{ 0 => Sign::Pos});

        let l = Link::unknot_r_twist();
        assert_eq!(l.iter_signed_crossings(), hashmap!{ 0 => Sign::Neg} );

        let l = Link::unknot_l_twist().resolve_at(0, Bit::Bit0);
        assert_eq!(l.iter_signed_crossings(), hashmap!{});
    }

    #[test]
    fn link_writhe() {
        let l = Link::unknot_l_twist();
        assert_eq!(l.writhe(), 1);

        let l = Link::unknot_r_twist();
        assert_eq!(l.writhe(), -1);

        let l = Link::unknot_l_twist().resolve_at(0, Bit::Bit0);
        assert_eq!(l.writhe(), 0);
    }

    #[test]
    fn link_components() {
        let l = Link::unknot_l_twist();
        let comps = l.comps();
        assert_eq!(comps, vec![ Path::new(vec![1, 2], true)]);
    }

    #[test]
    fn link_mirror() { 
        let l = Link::unknot_l_twist();
        assert_eq!(l.node(0).ntype(), XL);

        let l = l.mirror();
        assert_eq!(l.node(0).ntype(), XR);
    }

    #[test]
    fn link_resolve() {
        let s = State::from([0, 0, 0]);
        let l = Link::test_data("3_1").unwrap().resolve_by(&s);

        let comps = l.comps();
        assert_eq!(comps.len(), 3);
        assert!(comps.iter().all(|c| c.is_circle()));

        let s = State::from([1, 1, 1]);
        let l = Link::test_data("3_1").unwrap().resolve_by(&s);

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
        assert_eq!(l.n_crossings(), 0);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.n_comps(), 1);
    }

    #[test]
    fn trefoil() {
        let l = Link::test_data("3_1").unwrap();
        assert_eq!(l.n_crossings(), 3);
        assert_eq!(l.writhe(), -3);
        assert_eq!(l.n_comps(), 1);
    }

    #[test]
    fn figure8() {
        let l = Link::test_data("4_1").unwrap();
        assert_eq!(l.n_crossings(), 4);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.n_comps(), 1);
    }

    #[test]
    fn hopf_link() {
        let l = Link::test_data("L2a1").unwrap();
        assert_eq!(l.n_crossings(), 2);
        assert_eq!(l.writhe(), -2);
        assert_eq!(l.n_comps(), 2);
    }

    #[test]
    fn unlink_2() {
        let pd_code = [[1,2,3,4], [3,2,1,4]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.n_crossings(), 2);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.n_comps(), 2);
    }


    #[test]
    fn l2x4() {
        let pd_code = [[1,5,2,8],[5,3,6,2],[3,7,4,6],[7,1,8,4]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.n_crossings(), 4);
        assert_eq!(l.writhe(), 4);
        assert_eq!(l.n_comps(), 2);
    }

    #[test]
    fn crossing_change() {
        use crate::link::node::NodeOri;

        let l = Link::test_data("3_1").unwrap();
        let l2 = l.cc_at(1);

        // TODO must change. 
        assert_eq!(l.node(1),  &Node::new(XL, NodeOri::None, [3,6,4,1]));
        assert_eq!(l2.node(1), &Node::new(XR, NodeOri::None, [3,6,4,1]));
    }
}