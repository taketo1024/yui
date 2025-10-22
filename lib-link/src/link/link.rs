use std::collections::HashSet;
use std::fmt::Display;
use itertools::Itertools;
use yui::{CloneAnd, Sign};
use yui::bitseq::Bit;

use super::{Crossing, CrossingType, Path};

pub type Edge = usize;
pub type State = yui::bitseq::BitSeq;
pub type XCode = [Edge; 4];

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

#[derive(Debug, Clone)]
pub struct Link { 
    data: Vec<Crossing>
}

impl Link {
    pub fn new(data: Vec<Crossing>) -> Self { 
        let l = Self { data };
        l.validate();
        l
    }

    fn validate(&self) { 
        assert_eq!(self.edges().len(), self.data.len() * 2, "Invalid data.");
        assert!(self.components().iter().all(|c| c.is_circle()), "Non-closed components.")
    }

    pub fn from_pd_code<I>(pd_code: I) -> Self
    where I: IntoIterator<Item = XCode> { 
        let data = pd_code.into_iter().map(Crossing::from_pd_code).collect();
        Self::new(data)
    }

    pub fn empty() -> Link {
        Link { data: vec![] }
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn is_knot(&self) -> bool { 
        self.components().len() == 1
    }

    pub fn data(&self) -> &Vec<Crossing> { 
        &self.data
    }

    pub fn crossing_num(&self) -> usize { 
        self.data.iter()
            .filter(|x| !x.is_resolved())
            .count()
    }

    pub fn signed_crossing_nums(&self) -> (usize, usize) {
        let signs = self.crossing_signs().into_iter().counts();
        let pos = signs.get(&Sign::Pos).cloned().unwrap_or(0);
        let neg = signs.get(&Sign::Neg).cloned().unwrap_or(0);
        (pos, neg)
    }

    pub fn crossing_signs(&self) -> Vec<Sign> {
        use CrossingType::{X, Xm};

        let n = self.data.len();

        let mut signs = vec![None; n];
        let mut passed: HashSet<Edge> = HashSet::new();

        let mut traverse = |signs: &mut Vec<Option<Sign>>, j0: usize| {
            for i0 in 0..n {
                let start_edge = self.data[i0].edge(j0);
                if passed.contains(&start_edge) { 
                    continue 
                }

                self.traverse_from((i0, j0), |i, j| { 
                    let c = &self.data[i];
                    let e = c.edge(j);
                    passed.insert(e);

                    let sign = match (c.ctype(), j) { 
                        (Xm, 1) | (X, 3) => Some(Sign::Pos),
                        (Xm, 3) | (X, 1) => Some(Sign::Neg),
                        _                => None
                    };
                    if sign.is_some() { 
                        signs[i] = sign;
                    }
                });
            }
        };

        traverse(&mut signs, 0);

        if (0..n).any(|i| !self.data[i].is_resolved() && signs[i].is_none()) { 
            for j in [1,2] {
                traverse(&mut signs, j);
            }
        };

        let signs = signs.into_iter().flatten().collect_vec();

        assert_eq!(signs.len(), self.crossing_num());

        signs
    }
    
    pub fn writhe(&self) -> i32 { 
        let (p, n) = self.signed_crossing_nums();
        (p as i32) - (n as i32)
    }

    pub fn components(&self) -> Vec<Path> {
        let n = self.data.len();

        let mut comps = vec![];
        let mut passed: HashSet<Edge> = HashSet::new();

        let mut traverse = |j0: usize| {
            for i0 in 0..n {
                let start_edge = self.data[i0].edge(j0);
                if passed.contains(&start_edge) { 
                    continue 
                }

                let mut edges = vec![];

                self.traverse_from((i0, j0), |i, j| { 
                    let e = self.data[i].edge(j);
                    edges.push(e);
                    passed.insert(e);
                });

                let c =Path::circ(edges);

                comps.push(c);
            }
        };

        for i in [0, 1, 2] { 
            traverse(i);
        }

        comps
    }

    // index of data for the i-th `actual' crossing.
    fn crossing_index(&self, i: usize) -> usize {
        assert!(i < self.crossing_num());

        let mut i = i;
        for (j, x) in self.data.iter().enumerate() {
            if !x.is_resolved() {
                if i > 0 { 
                    i -= 1;
                } else { 
                    return j
                }
            }
        }
        
        panic!()
    }

    pub fn crossing_at(&self, i: usize) -> &Crossing {
        let j = self.crossing_index(i);
        &self.data[j]
    }

    pub fn crossing_at_mut(&mut self, i: usize) -> &mut Crossing {
        let j = self.crossing_index(i);
        &mut self.data[j]
    }

    pub fn crossing_changed_at(&self, i: usize) -> Self { 
        let j = self.crossing_index(i);

        let data = self.data.clone_and(|data|
            data[j] = data[j].mirror()
        );
        Link { data }
    }

    pub fn resolved_at(&self, i: usize, r: Bit) -> Self {
        debug_assert!(i < self.crossing_num());

        self.clone_and(|l|
            l.crossing_at_mut(i).resolve(r)
        )
    }

    pub fn resolved_by(&self, s: &State) -> Self {
        debug_assert!(s.len() == self.crossing_num());

        self.clone_and(|l|
            for r in s.iter() {
                l.crossing_at_mut(0).resolve(r)
            }
        )
    }

    pub fn mirror(&self) -> Self {
        let data = self.data.iter().map(|x| x.mirror()).collect();
        Link { data }
    }

    pub fn edges(&self) -> HashSet<Edge> {
        self.data.iter().flat_map(|x| x.edges()).cloned().collect()
    }

    pub fn first_edge(&self) -> Option<Edge> { 
        let x = self.data.first()?;
        x.edges().iter().min().cloned()
    }

    pub fn ori_pres_state(&self) -> State { 
        let signs = self.crossing_signs(); 
        State::from_iter(signs.into_iter().map( |s|
            if s.is_positive() { 0 } else { 1 }
        ))
    }

    pub fn seifert_circles(&self) -> Vec<Path> { 
        self.resolved_by(&self.ori_pres_state()).components()
    }

    fn traverse_from<F>(&self, start: (usize, usize), mut f:F) where
        F: FnMut(usize, usize)
    {
        let n = self.data.len();
        assert!((0..n).contains(&start.0));
        assert!((0..4).contains(&start.1));

        f(start.0, start.1); // call starting point

        let (mut i, mut j) = start;

        loop {
            let c = &self.data[i];
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
        let e = &self.data[c_index].edge(e_index);

        for (i, c) in self.data.iter().enumerate() { 
            for (j, f) in c.edges().iter().enumerate() { 
                if e == f && (c_index != i || (c_index == i && e_index != j)) { 
                    return (i, j)
                }
            }
        }

        panic!("Broken data")
    }
}

impl Link { 
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
        Link::from_pd_code([[0, 1, 1, 0]]).resolved_at(0, Bit::Bit0)
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
        write!(f, "L[{}]", self.data.iter().map(|x| x.to_string()).join(", "))
    }
}

#[cfg(test)]
mod tests { 
    use super::*;
    use super::CrossingType::{X, Xm};

    #[test]
    fn link_init() { 
        let l = Link { data: vec![] };
        assert_eq!(l.data.len(), 0);
    }

    #[test]
    fn link_from_pd_code() { 
        let pd_code = [[0,0,1,1]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.data.len(), 1);
        assert_eq!(l.data[0].ctype(), X);
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
        assert_eq!(l.crossing_num(), 0);

        let pd_code = [[0,0,1,1]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.crossing_num(), 1);
        
        let pd_code = [[1,4,2,5],[3,6,4,1],[5,2,6,3]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.crossing_num(), 3);
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
        assert_eq!(l.crossing_signs(), vec![Sign::Pos]);

        let pd_code = [[0,1,1,0]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.crossing_signs(), vec![Sign::Neg]);

        let pd_code = [[0,0,1,1]];
        let l = Link::from_pd_code(pd_code).resolved_at(0, Bit::Bit0);
        assert_eq!(l.crossing_signs(), vec![]);
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
        let l = Link::from_pd_code(pd_code).resolved_at(0, Bit::Bit0);
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
        assert_eq!(l.data[0].ctype(), X);

        let l = l.mirror();
        assert_eq!(l.data[0].ctype(), Xm);
    }

    #[test]
    fn link_resolve() {
        let s = State::from([0, 0, 0]);
        let l = Link::from_pd_code([[1,4,2,5],[3,6,4,1],[5,2,6,3]]) // trefoil
            .resolved_by(&s);

        let comps = l.components();
        assert_eq!(comps.len(), 3);
        assert!(comps.iter().all(|c| c.is_circle()));

        let s = State::from([1, 1, 1]);
        let l = Link::from_pd_code([[1,4,2,5],[3,6,4,1],[5,2,6,3]]) // trefoil
            .resolved_by(&s);

        let comps = l.components();
        assert_eq!(comps.len(), 2);
        assert!(comps.iter().all(|c| c.is_circle()));
    }

    #[test]
    fn empty_link() {
        let l = Link::empty();
        assert_eq!(l.crossing_num(), 0);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.components().len(), 0);
    }

    #[test]
    fn unknot() { 
        let l = Link::unknot();
        assert_eq!(l.crossing_num(), 0);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.components().len(), 1);
    }

    #[test]
    fn trefoil() { 
        let l = Link::trefoil();
        assert_eq!(l.crossing_num(), 3);
        assert_eq!(l.writhe(), -3);
        assert_eq!(l.components().len(), 1);
    }

    #[test]
    fn figure8() { 
        let l = Link::figure8();
        assert_eq!(l.crossing_num(), 4);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.components().len(), 1);
    }

    #[test]
    fn hopf_link() { 
        let l = Link::hopf_link();
        assert_eq!(l.crossing_num(), 2);
        assert_eq!(l.writhe(), -2);
        assert_eq!(l.components().len(), 2);
    }

    #[test]
    fn unlink_2() {
        let pd_code = [[1,2,3,4], [3,2,1,4]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.crossing_num(), 2);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.components().len(), 2);
    }


    #[test]
    fn l2x4() {
        let pd_code = [[1,5,2,8],[5,3,6,2],[3,7,4,6],[7,1,8,4]];
        let l = Link::from_pd_code(pd_code);
        assert_eq!(l.crossing_num(), 4);
        assert_eq!(l.writhe(), 4);
        assert_eq!(l.components().len(), 2);
    }

    #[test]
    fn load() { 
        let l = Link::load("3_1");
        assert!(l.is_ok());

        let l = l.unwrap();
        assert_eq!(l.crossing_num(), 3);
    }

    #[test]
    fn crossing_change() { 
        let l = Link::from_pd_code([[1,4,2,5],[3,6,4,1],[5,2,6,3]]);
        let l2 = l.crossing_changed_at(1);

        assert_eq!(l.data[1],  Crossing::new(CrossingType::X, [3,6,4,1]));
        assert_eq!(l2.data[1], Crossing::new(CrossingType::Xm, [3,6,4,1]));
    }
}