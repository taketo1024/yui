use core::panic;
use std::collections::HashSet;
use std::fmt::Display;
use itertools::Itertools;
use yui_core::{CloneAnd, Sign};
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
        let edge_counts = nodes.iter().flat_map(|x| x.edges()).cloned().counts();

        assert!(
            edge_counts.values().all(|&c| c == 2),
            "Invalid data: each edge must appear exactly twice."
        );

        let edges: HashSet<_> = edge_counts.into_keys().collect();
        Self { nodes, edges }
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
    where I: IntoIterator<Item = XCode> { 
        use crate::NodeOri::{Up, Right, None};
        
        let nodes = pd_code.into_iter().map(Node::from_pd_code).collect_vec();
        let mut l = Self::from_nodes(nodes); // unoriented
        
        let n = l.n_crossings();
        let mut ori = vec![None; l.n_nodes()];
        let mut remain = l.edges.clone();

        while !remain.is_empty() {
            // Take minimal edge-id. 
            let e0 = remain.iter().min().cloned().unwrap();

            // Find node & point where edge-id increases. 
            let (i0, j0) = (0..n).flat_map(|i| 
                [0usize, 1, 3].map(move |j| (i, j)) // no incoming from index 2
            ).find(|&(i, j)| 
                l.node(i).edge(j) == e0 && (l.node(i).counter_edge(j) == e0 + 1)
            ).unwrap();

            l.traverse_from((i0, j0), |i, j| { 
                let e = l.node(i).edge(j);
                if !remain.remove(&e) { 
                    panic!("Invalid data");
                }

                if j == 1 { 
                    ori[i] = Up 
                } else if j == 3 { 
                    ori[i] = Right
                }
            });
        }

        for (i, o) in ori.into_iter().enumerate() { 
            l.node_mut(i).ori = o;
        }

        assert!(l.is_oriented());

        l
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
        let (p, n) = self.n_signed_crossings();
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

    pub fn n_signed_crossings(&self) -> (usize, usize) {
        let mut pos = 0;
        let mut neg = 0;
        for n in self.nodes.iter() { 
            if n.is_pos() { pos += 1 } 
            else if n.is_neg() { neg += 1}
        }
        (pos, neg)
    }

    pub fn is_oriented(&self) -> bool { 
        self.nodes().all(|n| n.is_oriented())
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
        self.traverse_comps(|c, _, _| 
            if count <= c { count = c + 1 } 
        );
        count
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

        comps.into_iter().map(|edges| 
            Path::circ(edges)
        ).collect()
    }

    fn traverse_comps<F>(&self, mut f: F) where 
    F: FnMut(usize, usize, usize) { 
        let n = self.n_nodes();

        let mut c = 0; // component counter
        let mut remain = self.edges.clone();

        while !remain.is_empty() {
            // Take minimal edge-id. 
            let e0 = remain.iter().min().cloned().unwrap();

            // Find node & point having edge e0. 
            let (i0, j0) = (0..n).flat_map(|i| 
                (0usize..4).map(move |j| (i, j))
            ).find(|&(i, j)| 
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
        // MEMO: no assertion here since `unknot` is not oriented. 
        // assert!(self.is_oriented());

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
        assert_eq!(l.node(0).ntype(), XL);
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
        assert_eq!(l.node(0).ntype(), XL);

        let l = l.mirror();
        assert_eq!(l.node(0).ntype(), XR);
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
        assert_eq!(l.n_crossings(), 0);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.n_comps(), 1);
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
        let l = Link::test_data("unlink2");
        assert_eq!(l.n_crossings(), 2);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.n_comps(), 2);
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
}