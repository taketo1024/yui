use core::panic;
use std::collections::HashSet;
use std::fmt::Display;
use itertools::Itertools;
use yui::{CloneAnd, Sign};
use yui::bitseq::Bit;

use super::{Node, NodeType, Path};

pub type Edge = usize;
pub type State = yui::bitseq::BitSeq;
pub type XCode = [Edge; 4];

#[derive(Debug, Clone)]
pub struct Link { 
    nodes: Vec<Node>
}

impl Link {
    pub fn new(nodes: Vec<Node>) -> Self { 
        let l = Self { nodes };
        l.validate();
        l
    }

    fn validate(&self) { 
        assert_eq!(self.edges().len(), self.nodes.len() * 2, "Invalid data.");
        assert!(self.components().iter().all(|c| c.is_circle()), "Non-closed components.")
    }

    pub fn empty() -> Link {
        Link { nodes: vec![] }
    }

    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    pub fn is_knot(&self) -> bool { 
        self.components().len() == 1
    }

    pub fn writhe(&self) -> i32 { 
        let (p, n) = self.count_signed_crossings();
        (p as i32) - (n as i32)
    }

    pub fn mirror(&self) -> Self {
        self.clone_and(|l|
            l.nodes.iter_mut().for_each(|x| x.cc())
        )
    }

    pub fn nodes(&self) -> &Vec<Node> { 
        &self.nodes
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

    pub fn count_crossings(&self) -> usize { 
        self.nodes.iter()
            .filter(|x| x.is_crossing())
            .count()
    }

    pub fn count_signed_crossings(&self) -> (usize, usize) {
        let signs = self.collect_crossing_signs().into_iter().counts();
        let pos = signs.get(&Sign::Pos).cloned().unwrap_or(0);
        let neg = signs.get(&Sign::Neg).cloned().unwrap_or(0);
        (pos, neg)
    }

    pub fn collect_crossing_signs(&self) -> Vec<Sign> {
        use NodeType::{X, Xm};

        let n = self.nodes.len();
        let mut signs = vec![None; n];
        let mut remain = self.edges();

        for j0 in [0, 1, 2] { 
            for i0 in 0..n {
                let start_edge = self.nodes[i0].edge(j0);
                if !remain.remove(&start_edge) { 
                    continue 
                }

                self.traverse_from((i0, j0), |i, j| { 
                    let c = &self.nodes[i];
                    let e = c.edge(j);
                    remain.remove(&e);

                    let sign = match (c.ntype(), j) { 
                        (Xm, 1) | (X, 3) => Some(Sign::Pos),
                        (Xm, 3) | (X, 1) => Some(Sign::Neg),
                        _                => None
                    };
                    if sign.is_some() { 
                        signs[i] = sign;
                    }
                });
            }

            if remain.is_empty() { 
                break
            }
        }

        assert!(remain.is_empty());

        signs.into_iter().flatten().collect_vec()
    }
    
    pub fn edges(&self) -> HashSet<Edge> {
        self.nodes.iter().flat_map(|x| x.edges()).cloned().collect()
    }

    pub fn first_edge(&self) -> Option<Edge> { 
        let x = self.nodes.first()?;
        x.edges().iter().min().cloned()
    }

    pub fn components(&self) -> Vec<Path> {
        let n = self.nodes.len();

        let mut comps = vec![];
        let mut remain = self.edges();

        for j0 in [0, 1, 2] {
            for i0 in 0..n {
                let start_edge = self.node(i0).edge(j0);
                if !remain.remove(&start_edge) { 
                    continue 
                }

                let mut edges = vec![];

                self.traverse_from((i0, j0), |i, j| { 
                    let e = self.node(i).edge(j);
                    edges.push(e);
                    remain.remove(&e);
                });

                let c = Path::circ(edges);

                comps.push(c);
            }

            if remain.is_empty() { 
                break
            }
        }

        assert!(remain.is_empty());

        comps
    }

    pub fn crossing_change(&self, i: usize) -> Self { 
        assert!(self.node(i).is_crossing());
        self.clone_and(|l| l.node_mut(i).cc())
    }

    pub fn resolve_crossing(&self, i: usize, r: Bit) -> Self {
        assert!(self.node(i).is_crossing());
        self.clone_and(|l| l.node_mut(i).resolve(r))
    }

    pub fn resolve_all_crossings(&self, s: &State) -> Self {
        assert!(s.len() == self.count_crossings());

        let n = self.nodes.len();
        let itr = (0..n).filter(|&i| self.node(i).is_crossing());

        self.clone_and(|l| {
            for (i, r) in Iterator::zip(itr, s.iter()) {
                l.node_mut(i).resolve(r); 
            }
        })
    }

    pub fn ori_pres_state(&self) -> State { 
        let signs = self.collect_crossing_signs(); 
        State::from_iter(signs.into_iter().map( |s|
            if s.is_positive() { 0 } else { 1 }
        ))
    }

    pub fn seifert_circles(&self) -> Vec<Path> { 
        self.resolve_all_crossings(&self.ori_pres_state()).components()
    }

    fn traverse_from<F>(&self, start: (usize, usize), mut f:F) where
        F: FnMut(usize, usize)
    {
        let n = self.nodes.len();
        assert!((0..n).contains(&start.0));
        assert!((0..4).contains(&start.1));

        f(start.0, start.1); // call starting point

        let (mut i, mut j) = start;

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

    fn traverse_outer(&self, c_index:usize, e_index:usize) -> (usize, usize) {
        let e = self.nodes[c_index].edge(e_index);

        for (i, c) in self.nodes.iter().enumerate() { 
            for (j, &f) in c.edges().iter().enumerate() { 
                if e == f && (c_index != i || (c_index == i && e_index != j)) { 
                    return (i, j)
                }
            }
        }

        panic!("Broken data")
    }
}

impl Link { 
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
        let nodes = pd_code.into_iter().map(Node::from_pd_code).collect();
        Self::new(nodes)
    }

    pub fn is_valid_name(str: &str) -> bool { 
        use regex::Regex;
        let r1 = Regex::new(r"^([1-9]|10)_[0-9]+$").unwrap();
        let r2 = Regex::new(r"^(K|L)?[1-9]+(a|n)_?[0-9]+$").unwrap(); // FIXME tmp
        r1.is_match(str) || r2.is_match(str)
    }

    pub fn load(name_or_path: &str) -> Result<Link, Box<dyn std::error::Error>> {
        const RESOURCE_DIR: &str = "resources/links/";
        
        if Self::is_valid_name(name_or_path) { 
            let dir = std::env!("CARGO_MANIFEST_DIR");
            let path = format!("{dir}/{RESOURCE_DIR}{name_or_path}.json");
            Self::_load(&path)
        } else { 
            Self::_load(name_or_path)
        }
    }

    fn _load(path: &str) -> Result<Link, Box<dyn std::error::Error>> {
        let json = std::fs::read_to_string(path)?;
        let data: Vec<XCode> = serde_json::from_str(&json)?;
        let l = Link::from_pd_code(data);
        Ok(l)
    }

    pub fn unknot() -> Link { 
        Link::from_pd_code([[0, 1, 1, 0]]).resolve_crossing(0, Bit::Bit0)
    }

    pub fn trefoil() -> Link { 
        Link::from_pd_code([[1,4,2,5],[3,6,4,1],[5,2,6,3]])
    }

    pub fn figure8() -> Link { 
        Link::from_pd_code([[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]])
    }

    pub fn hopf_link() -> Link { 
        Link::from_pd_code([[4,1,3,2],[2,3,1,4]])
    }
}

impl Display for Link {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "L[{}]", self.nodes.iter().map(|x| x.to_string()).join(", "))
    }
}

#[cfg(test)]
mod tests { 
    use super::*;
    use super::NodeType::{X, Xm};

    #[test]
    fn link_init() { 
        let l = Link { nodes: vec![] };
        assert_eq!(l.nodes.len(), 0);
    }

    #[test]
    fn link_from_pd_code() { 
        let pd_code = [[0,0,1,1]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.nodes.len(), 1);
        assert_eq!(l.node(0).ntype(), X);
    }

    #[test]
    fn link_is_empty() {
        let l = Link::empty();
        assert!(l.is_empty());

        let pd_code = [[0,0,1,1]];
        let l = Link::from_pd_code(pd_code);
        assert!(!l.is_empty());
    }

    #[test]
    fn link_crossing_num() {
        let l = Link::empty();
        assert_eq!(l.count_crossings(), 0);

        let pd_code = [[0,0,1,1]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.count_crossings(), 1);
        
        let pd_code = [[1,4,2,5],[3,6,4,1],[5,2,6,3]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.count_crossings(), 3);
    }

    #[test]
    fn link_next() {
        let pd_code = [[0,0,1,1]];
        let l = Link::from_pd_code(pd_code);

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

        let pd_code = [[0,0,1,1]];
        let l = Link::from_pd_code(pd_code);
        let path = traverse(&l, (0, 0));
        
        assert_eq!(path, [(0, 0), (0, 3)]); // loop
    }

    #[test]
    fn link_crossing_signs() {
        let pd_code = [[0,0,1,1]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.collect_crossing_signs(), vec![Sign::Pos]);

        let pd_code = [[0,1,1,0]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.collect_crossing_signs(), vec![Sign::Neg]);

        let pd_code = [[0,0,1,1]];
        let l = Link::from_pd_code(pd_code).resolve_crossing(0, Bit::Bit0);
        assert_eq!(l.collect_crossing_signs(), vec![]);
    }

    #[test]
    fn link_writhe() {
        let pd_code = [[0,0,1,1]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.writhe(), 1);

        let pd_code = [[0,1,1,0]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.writhe(), -1);

        let pd_code = [[0,0,1,1]];
        let l = Link::from_pd_code(pd_code).resolve_crossing(0, Bit::Bit0);
        assert_eq!(l.writhe(), 0);

    }

    #[test]
    fn link_components() {
        let pd_code = [[0,0,1,1]];
        let l = Link::from_pd_code(pd_code);
        let comps = l.components();
        assert_eq!(comps, vec![ Path::new(vec![0, 1], true)]);
    }

    #[test]
    fn link_mirror() { 
        let pd_code = [[0,0,1,1]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.node(0).ntype(), X);

        let l = l.mirror();
        assert_eq!(l.node(0).ntype(), Xm);
    }

    #[test]
    fn link_resolve() {
        let s = State::from([0, 0, 0]);
        let l = Link::from_pd_code([[1,4,2,5],[3,6,4,1],[5,2,6,3]]) // trefoil
            .resolve_all_crossings(&s);

        let comps = l.components();
        assert_eq!(comps.len(), 3);
        assert!(comps.iter().all(|c| c.is_circle()));

        let s = State::from([1, 1, 1]);
        let l = Link::from_pd_code([[1,4,2,5],[3,6,4,1],[5,2,6,3]]) // trefoil
            .resolve_all_crossings(&s);

        let comps = l.components();
        assert_eq!(comps.len(), 2);
        assert!(comps.iter().all(|c| c.is_circle()));
    }

    #[test]
    fn empty_link() {
        let l = Link::empty();
        assert_eq!(l.count_crossings(), 0);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.components().len(), 0);
    }

    #[test]
    fn unknot() { 
        let l = Link::unknot();
        assert_eq!(l.count_crossings(), 0);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.components().len(), 1);
    }

    #[test]
    fn trefoil() { 
        let l = Link::trefoil();
        assert_eq!(l.count_crossings(), 3);
        assert_eq!(l.writhe(), -3);
        assert_eq!(l.components().len(), 1);
    }

    #[test]
    fn figure8() { 
        let l = Link::figure8();
        assert_eq!(l.count_crossings(), 4);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.components().len(), 1);
    }

    #[test]
    fn hopf_link() { 
        let l = Link::hopf_link();
        assert_eq!(l.count_crossings(), 2);
        assert_eq!(l.writhe(), -2);
        assert_eq!(l.components().len(), 2);
    }

    #[test]
    fn unlink_2() {
        let pd_code = [[1,2,3,4], [3,2,1,4]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.count_crossings(), 2);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.components().len(), 2);
    }


    #[test]
    fn l2x4() {
        let pd_code = [[1,5,2,8],[5,3,6,2],[3,7,4,6],[7,1,8,4]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.count_crossings(), 4);
        assert_eq!(l.writhe(), 4);
        assert_eq!(l.components().len(), 2);
    }

    #[test]
    fn load() { 
        let l = Link::load("3_1");
        assert!(l.is_ok());

        let l = l.unwrap();
        assert_eq!(l.count_crossings(), 3);
    }

    #[test]
    fn crossing_change() { 
        let l = Link::from_pd_code([[1,4,2,5],[3,6,4,1],[5,2,6,3]]);
        let l2 = l.crossing_change(1);

        assert_eq!(l.node(1),  &Node::new(NodeType::X, [3,6,4,1]));
        assert_eq!(l2.node(1), &Node::new(NodeType::Xm, [3,6,4,1]));
    }
}