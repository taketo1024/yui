use std::{collections::HashSet, ops::Index};
use std::fmt;
use itertools::{join, Itertools};
use Resolution::{Res0, Res1};
use CrossingType::{Xp, Xn, V, H};

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
    pub fn empty() -> Link {
        Link { data: vec![] }
    }

    pub fn unknot() -> Link { 
        Link::from(&[[0, 1, 1, 0]]).resolved_at(0, Res0)
    }

    pub fn trefoil() -> Link { 
        Link::from(&[[1,4,2,5],[3,6,4,1],[5,2,6,3]])
    }

    pub fn figure8() -> Link { 
        Link::from(&[[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]])
    }

    pub fn hopf_link() -> Link { 
        Link::from(&[[4,1,3,2],[2,3,1,4]])
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn crossing_num(&self) -> u32 { 
        self.data.iter()
            .filter(|x| !x.is_resolved())
            .count() as u32
    }

    pub fn signed_crossing_nums(&self) -> (u32, u32) {
        let signs = self.crossing_signs();
        let pos = signs.iter().filter(|e| e.is_positive()).count() as u32;
        let neg = signs.iter().filter(|e| e.is_negative()).count() as u32;
        (pos, neg)
    }

    pub fn writhe(&self) -> i32 { 
        self.crossing_signs().iter().sum()
    }

    pub fn components(&self) -> Vec<Component> {
        let n = self.data.len();

        let mut comps = vec![];
        let mut passed: HashSet<Edge> = HashSet::new();

        let mut traverse = |e0: usize| {
            for i0 in 0..n {
                let start_edge = self.data[i0].edges[e0];
                if passed.contains(&start_edge) { 
                    continue 
                }

                let mut edges = vec![];
                let mut closed = false;

                self.traverse_edges(i0, e0, |i, j| { 
                    let edge = self.data[i].edges[j];
                    if edges.first() == Some(&edge) { 
                        closed = true;
                    } else { 
                        edges.push(edge);
                        passed.insert(edge);
                    }
                });

                let comp = Component{edges, closed};
                comps.push(comp);
            }
        };

        for i in [0, 1, 2] { 
            traverse(i);
        }

        comps
    }

    pub fn mirror(&self) -> Self {
        let mut data = self.data.clone();
        data.iter_mut().for_each(|x| x.mirror());
        Link { data }
    }

    pub fn resolved_at(&self, i: usize, r: Resolution) -> Self {
        debug_assert!(i < self.data.len());
        let mut data = self.data.clone();
        data[i].resolve(r);
        Link { data }
    }

    pub fn resolved_by(&self, s: &State) -> Self {
        debug_assert!(s.len() <= self.data.len());
        let mut data = self.data.clone();
        for (i, &r) in s.values.iter().enumerate() {
            data[i].resolve(r);
        }
        Link { data }
    }

    pub fn ori_pres_state(&self) -> State { 
        State::from_iter(self.crossing_signs().into_iter().map( |e|
            if e.is_positive() { 0 } else { 1 }
        ))
    }

    pub fn seifert_circles(&self) -> Vec<Component> { 
        self.resolved_by(&self.ori_pres_state()).components()
    }

    pub fn first_edge(&self) -> Option<&Edge> { 
        let Some(x) = self.data.first() else { return None };
        x.edges.iter().min()
    }

    // -- internal methods -- //
    
    fn next(&self, c_index:usize, e_index:usize) -> Option<(usize, usize)> {
        let n = self.data.len();
        debug_assert!((0..n).contains(&c_index));
        debug_assert!((0..4).contains(&e_index));

        let e = &self.data[c_index].edges[e_index];

        for (i, c) in self.data.iter().enumerate() { 
            for (j, f) in c.edges.iter().enumerate() { 
                if e == f && (c_index != i || (c_index == i && e_index != j)) { 
                    return Some((i, j))
                }
            }
        }
        None
    }

    fn traverse_edges<F>(&self, c_index:usize, e_index:usize, mut f:F) where
        F: FnMut(usize, usize)
    {
        let n = self.data.len();
        debug_assert!((0..n).contains(&c_index));
        debug_assert!((0..4).contains(&e_index));

        let mut i = c_index;
        let mut j = e_index;

        loop {
            f(i, j);

            let c = &self.data[i];
            let k = c.pass(j);
            
            if let Some(next) = self.next(i, k) { 
                if next != (c_index, e_index) {
                    (i, j) = next;
                    continue
                } else {
                    f(c_index, e_index); // final call
                }
            } else { 
                f(i, k); // end edge
            }
            break
        }
    }

    fn crossing_signs(&self) -> Vec<i32> {
        let n = self.data.len();
        let mut signs = vec![0; n];

        let mut traverse = |e0: usize| {
            for i0 in 0..n {
                if signs[i0] != 0 { 
                    continue 
                }
                self.traverse_edges(i0, e0, |i, j| { 
                    let c = &self.data[i];
                    let sign = match (c.ctype, j) { 
                        (Xp, 1) | (Xn, 3) =>  1,
                        (Xp, 3) | (Xn, 1) => -1,
                        _ => 0
                    };
                    if sign != 0 {
                        signs[i] = sign;
                    }
                })
            }
        };

        traverse(0);
        traverse(1); // in case 

        signs
    }
}

impl<const N: usize> From<&[[Edge; 4]; N]> for Link { 
    fn from(pd_code: &[[Edge; 4]; N]) -> Link {
        let data = pd_code.iter().cloned().collect_vec(); 
        Self::from(&data)
    }
}

impl From<&Vec<[Edge; 4]>> for Link { 
    fn from(pd_code: &Vec<[Edge; 4]>) -> Link { 
        let data = pd_code.into_iter().map( |x| 
            Crossing { ctype: Xn, edges: x.clone() }
        ).collect();
        Link{ data }
    }
}

impl Link { 
    pub fn load(name: &str) -> Result<Link, Box<dyn std::error::Error>> {
        let path = format!("{}/links/{}.json", crate::RESOURCE_DIR, name);
        let json = std::fs::read_to_string(path)?;
        let data: Vec<[Edge; 4]> = serde_json::from_str(&json)?;
        let l = Link::from(&data);
        Ok(l)
    }
}

pub type Edge = u8;

#[derive(Debug, Clone, Copy, PartialEq)]
enum CrossingType { 
    Xp, Xn, V, H 
}

impl CrossingType { 
    fn mirror(self) -> CrossingType {
        match self { 
            Xp => Xn,
            Xn => Xp, 
            other => other
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Crossing { 
    ctype: CrossingType,
    edges: [Edge; 4]
}

impl Crossing {
    fn is_resolved(&self) -> bool { 
        match self.ctype { 
            V | H => true,
            _ => false
        }
    }

    fn resolve(&mut self, r: Resolution) {
        match (self.ctype, r) {
            (Xp, Res0) | (Xn, Res1) => self.ctype = V,
            (Xp, Res1) | (Xn, Res0) => self.ctype = H,
            _ => panic!()
        }
    }

    fn mirror(&mut self) { 
        self.ctype = self.ctype.mirror();
    }

    fn pass(&self, index:usize) -> usize { 
        debug_assert!((0..4).contains(&index));

        match self.ctype {
            Xp | Xn => (index + 2) % 4,
            V => 3 - index,
            H => (5 - index) % 4
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Component { 
    edges:Vec<Edge>,
    closed:bool
}

impl Component { 
    pub fn edges(&self) -> &Vec<Edge> { 
        &self.edges
    }

    pub fn is_adj(&self, other: &Component, link: &Link) -> bool { 
        // 1) find crossings `x` that touche `self`. 
        // 2) check if `x` also touches `other`.
        
        for x in &link.data { 
            if !x.edges.iter().any(|e| 
                self.edges.contains(e)
            ) { 
                continue
            }

            let Some(e) = x.edges.iter().find(|e| 
                !self.edges.contains(e)
            ) else { 
                continue
            };

            if other.edges.contains(e) { 
                return true
            }
        }

        false
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Resolution { 
    Res0, Res1
}

impl Resolution { 
    pub fn as_u8(&self) -> u8 {
        match self { 
            Res0 => 0,
            Res1 => 1
        }
    }

    pub fn is_zero(&self) -> bool {
        self == &Res0
    }
}

impl From<u8> for Resolution {
    fn from(a: u8) -> Self {
        match a { 
            0 => Res0,
            1 => Res1,
            _ => panic!()
        }
    }
}

impl fmt::Display for Resolution {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self { 
            Res0 => f.write_str("0"),
            Res1 => f.write_str("1")
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct State { 
    values: Vec<Resolution>
}

impl State { 
    pub fn empty() -> Self { 
        State { values: vec![] }
    }

    pub fn from_iter<I>(iter: I) -> Self
    where I: Iterator<Item = u8> {
        let values = iter.map(|a| Resolution::from(a)).collect();
        State{ values }
    }

    pub fn from_bseq(mut bseq: usize, length: usize) -> Self {
        let seq = (0..length)
            .map(|_| {
                let a = (bseq & 1) as u8;
                bseq >>= 1;
                Resolution::from(a)
            })
            .rev()
            .collect();
        State{ values: seq }
    }

    pub fn values(&self) -> core::slice::Iter<Resolution> {
        self.values.iter()
    }

    pub fn weight(&self) -> usize { 
        self.values.iter().filter(|r| !r.is_zero()).count()
    }
    
    pub fn len(&self) -> usize { 
        self.values.len()
    }

    pub fn targets(&self) -> Vec<State> { 
        let n = self.len();
        (0..n).filter(|&i| self[i] == Res0 ).map(|i| { 
            let mut t = self.clone();
            t.values[i] = Res1;
            t
        }).collect_vec()
    }

    pub fn append(&mut self, mut other: Self) {
        self.values.append(&mut other.values);
    }
}

impl<const N: usize> From<[u8; N]> for State {
    fn from(values: [u8; N]) -> Self {
        Self::from_iter(values.into_iter())
    }
}

impl From<Vec<u8>> for State { 
    fn from(v: Vec<u8>) -> Self {
        Self::from_iter(v.into_iter())
    }
}

impl Index<usize> for State { 
    type Output = Resolution;
    fn index(&self, index: usize) -> &Resolution {
        &self.values[index]
    }
}

impl fmt::Display for State {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{}]", join(self.values.iter(), ", "))
    }
}


#[cfg(test)]
mod tests { 
    use super::*;
    use super::Resolution::{Res0, Res1};
    use super::CrossingType::{Xp, Xn, H, V};

    fn a_crossing(ctype:CrossingType) -> Crossing {
        Crossing{
            ctype, 
            edges: [0,1,2,3]
        }
    }

    #[test]
    fn crossing_is_resolved() {
        let c = a_crossing(Xp);
        assert!(!c.is_resolved());

        let c = a_crossing(Xn);
        assert!(!c.is_resolved());

        let c = a_crossing(H);
        assert!(c.is_resolved());

        let c = a_crossing(V);
        assert!(c.is_resolved());
    }

    #[test]
    fn crossing_resolve() {
        let mut c = a_crossing(Xp);
        
        c.resolve(Res0);
        assert!(c.is_resolved());
        assert_eq!(c.ctype, V);

        let mut c = a_crossing(Xp);
        
        c.resolve(Res1);
        assert!(c.is_resolved());
        assert_eq!(c.ctype, H);

        let mut c = a_crossing(Xn);
        
        c.resolve(Res0);
        assert!(c.is_resolved());
        assert_eq!(c.ctype, H);

        let mut c = a_crossing(Xn);
        
        c.resolve(Res1);
        assert!(c.is_resolved());
        assert_eq!(c.ctype, V);
    }

    #[test]
    fn crossing_mirror() {
        let mut c = a_crossing(Xp);
        c.mirror();
        assert_eq!(c.ctype, Xn);

        let mut c = a_crossing(Xn);
        c.mirror();
        assert_eq!(c.ctype, Xp);

        let mut c = a_crossing(H);
        c.mirror();
        assert_eq!(c.ctype, H);

        let mut c = a_crossing(V);
        c.mirror();
        assert_eq!(c.ctype, V);
    }

    #[test]
    fn crossing_pass() {
        let c = a_crossing(Xp);
        assert_eq!(c.pass(0), 2);
        assert_eq!(c.pass(1), 3);
        assert_eq!(c.pass(2), 0);
        assert_eq!(c.pass(3), 1);

        let c = a_crossing(Xn);
        assert_eq!(c.pass(0), 2);
        assert_eq!(c.pass(1), 3);
        assert_eq!(c.pass(2), 0);
        assert_eq!(c.pass(3), 1);

        let c = a_crossing(V);
        assert_eq!(c.pass(0), 3);
        assert_eq!(c.pass(1), 2);
        assert_eq!(c.pass(2), 1);
        assert_eq!(c.pass(3), 0);

        let c = a_crossing(H);
        assert_eq!(c.pass(0), 1);
        assert_eq!(c.pass(1), 0);
        assert_eq!(c.pass(2), 3);
        assert_eq!(c.pass(3), 2);
    }

    #[test]
    fn link_init() { 
        let l = Link { data: vec![] };
        assert_eq!(l.data.len(), 0);
    }

    #[test]
    fn link_from_pd_code() { 
        let pd_code = [[1,2,3,4]];
        let l = Link::from(&pd_code);
        assert_eq!(l.data.len(), 1);
        assert_eq!(l.data[0].ctype, Xn);
    }

    #[test]
    fn link_is_empty() {
        let l = Link::empty();
        assert!(l.is_empty());

        let pd_code = [[1,2,3,4]];
        let l = Link::from(&pd_code);
        assert!(!l.is_empty());
    }

    #[test]
    fn link_crossing_num() {
        let l = Link::empty();
        assert_eq!(l.crossing_num(), 0);

        let pd_code = [[1,2,3,4]];
        let l = Link::from(&pd_code);
        assert_eq!(l.crossing_num(), 1);
    }

    #[test]
    fn link_next() {
        let pd_code = [[0,0,1,1]];
        let l = Link::from(&pd_code);

        assert_eq!(l.next(0, 0), Some((0, 1)));
        assert_eq!(l.next(0, 1), Some((0, 0)));
        assert_eq!(l.next(0, 2), Some((0, 3)));
        assert_eq!(l.next(0, 3), Some((0, 2)));

        let pd_code = [[0,3,1,4],[3,2,2,1]];
        let l = Link::from(&pd_code);

        assert_eq!(l.next(0, 0), None);
        assert_eq!(l.next(0, 2), Some((1, 3)));
        assert_eq!(l.next(1, 1), Some((1, 2)));
        assert_eq!(l.next(1, 0), Some((0, 1)));
        assert_eq!(l.next(0, 3), None);
    }

    #[test]
    fn link_traverse() {
        let pd_code = [[0,0,1,1]];
        let l = Link::from(&pd_code);
        let mut queue = vec![];
        l.traverse_edges(0, 0, |i, j| queue.push((i, j)));
        
        assert_eq!(queue, [(0, 0), (0, 3), (0, 0)]); // loop

        let pd_code = [[0,3,1,4],[3,2,2,1]];
        let l = Link::from(&pd_code);
        let mut queue = vec![];
        l.traverse_edges(0, 0, |i, j| queue.push((i, j)));

        assert_eq!(queue, [(0, 0), (1, 3), (1, 2), (0, 1), (0, 3)]); // no loop
    }

    #[test]
    fn link_crossing_signs() {
        let pd_code = [[0,0,1,1]];
        let l = Link::from(&pd_code);
        assert_eq!(l.crossing_signs(), vec![1]);

        let pd_code = [[0,1,1,0]];
        let l = Link::from(&pd_code);
        assert_eq!(l.crossing_signs(), vec![-1]);
    }

    #[test]
    fn link_writhe() {
        let pd_code = [[0,0,1,1]];
        let l = Link::from(&pd_code);
        assert_eq!(l.writhe(), 1);

        let pd_code = [[0,1,1,0]];
        let l = Link::from(&pd_code);
        assert_eq!(l.writhe(), -1);
    }

    #[test]
    fn link_components() {
        let pd_code = [[0,0,1,1]];
        let l = Link::from(&pd_code);
        let comps = l.components();
        assert_eq!(comps, vec![Component{ edges: vec![0, 1], closed: true }]);

        let pd_code = [[0,3,1,4],[3,2,2,1]];
        let l = Link::from(&pd_code);
        let comps = l.components();
        assert_eq!(comps, vec![Component{ edges: vec![0,1,2,3,4], closed: false }]);
    }

    #[test]
    fn link_mirror() { 
        let pd_code = [[0,0,1,1]];
        let l = Link::from(&pd_code);
        assert_eq!(l.data[0].ctype, Xn);

        let l = l.mirror();
        assert_eq!(l.data[0].ctype, Xp);
    }

    #[test]
    fn link_resolve() {
        let s = State::from([0, 0, 0]);
        let l = Link::from(&[[1,4,2,5],[3,6,4,1],[5,2,6,3]]) // trefoil
            .resolved_by(&s);

        let comps = l.components();
        assert_eq!(comps.len(), 3);
        assert!(comps.iter().all(|c| c.closed));

        let s = State::from([1, 1, 1]);
        let l = Link::from(&[[1,4,2,5],[3,6,4,1],[5,2,6,3]]) // trefoil
            .resolved_by(&s);

        let comps = l.components();
        assert_eq!(comps.len(), 2);
        assert!(comps.iter().all(|c| c.closed));
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
        let l = Link::from(&pd_code);
        assert_eq!(l.crossing_num(), 2);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.components().len(), 2);
    }

    #[test]
    fn state_targets() {
        let s = State::from(vec![0, 0, 0]);
        assert_eq!(s.targets(), vec![
            State::from(vec![1, 0, 0]),
            State::from(vec![0, 1, 0]),
            State::from(vec![0, 0, 1]),
        ]);

        let s = State::from(vec![0, 1, 0]);
        assert_eq!(s.targets(), vec![
            State::from(vec![1, 1, 0]),
            State::from(vec![0, 1, 1]),
        ]);

        let s = State::from(vec![1, 1, 0]);
        assert_eq!(s.targets(), vec![
            State::from(vec![1, 1, 1]),
        ]);

        let s = State::from(vec![1, 1, 1]);
        assert_eq!(s.targets(), vec![]);
    }
}