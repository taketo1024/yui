use std::collections::HashSet;
use itertools::Itertools;
use super::{Crossing, CrossingType, Resolution, LinkComp, State};

pub type Edge = u8;

use CrossingType::{Xp, Xn};

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
        Link::from(&[[0, 1, 1, 0]]).resolved_at(0, Resolution::Res0)
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

    pub fn components(&self) -> Vec<LinkComp> {
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

                self.traverse_edges((i0, j0), |i, j| { 
                    let e = self.data[i].edge(j);
                    edges.push(e);
                    passed.insert(e);
                });

                let c = if edges.len() > 1 && edges.first() == edges.last() { 
                    edges.pop();
                    LinkComp::circ(edges)
                } else { 
                    LinkComp::arc(edges)
                };

                comps.push(c);
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
        for (i, &r) in s.values().iter().enumerate() {
            data[i].resolve(r);
        }
        Link { data }
    }

    pub fn ori_pres_state(&self) -> State { 
        State::from_iter(self.crossing_signs().into_iter().map( |e|
            if e.is_positive() { 0 } else { 1 }
        ))
    }

    pub fn seifert_circles(&self) -> Vec<LinkComp> { 
        self.resolved_by(&self.ori_pres_state()).components()
    }

    pub fn first_edge(&self) -> Option<&Edge> { 
        let Some(x) = self.data.first() else { return None };
        x.edges().iter().min()
    }

    pub fn data(&self) -> &Vec<Crossing> { 
        &self.data
    }

    // -- internal methods -- //
    
    fn next(&self, c_index:usize, e_index:usize) -> Option<(usize, usize)> {
        let n = self.data.len();
        debug_assert!((0..n).contains(&c_index));
        debug_assert!((0..4).contains(&e_index));

        let e = &self.data[c_index].edge(e_index);

        for (i, c) in self.data.iter().enumerate() { 
            for (j, f) in c.edges().iter().enumerate() { 
                if e == f && (c_index != i || (c_index == i && e_index != j)) { 
                    return Some((i, j))
                }
            }
        }
        None
    }

    fn traverse_edges<F>(&self, start: (usize, usize), mut f:F) where
        F: FnMut(usize, usize)
    {
        let n = self.data.len();
        debug_assert!((0..n).contains(&start.0));
        debug_assert!((0..4).contains(&start.1));

        let (mut i, mut j) = start;

        loop {
            f(i, j);

            let c = &self.data[i];
            let k = c.pass(j);
            
            let Some(next) = self.next(i, k) else {
                // reached end
                f(i, k);
                break
            } ;

            if next == start {
                // returned to start
                f(start.0, start.1);
                break
            }

            (i, j) = next;
        }
    }

    fn crossing_signs(&self) -> Vec<i32> {
        let n = self.data.len();
        let mut signs = vec![0; n];

        let mut traverse = |j0: usize| {
            for i0 in 0..n {
                if signs[i0] != 0 { 
                    continue 
                }

                self.traverse_edges((i0, j0), |i, j| { 
                    let c = &self.data[i];
                    let sign = match (c.ctype(), j) { 
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
            Crossing::new(Xn, x.clone())
        ).collect();
        Link{ data }
    }
}

impl Link { 
    pub fn load(path: &str) -> Result<Link, Box<dyn std::error::Error>> {
        let json = std::fs::read_to_string(path)?;
        let data: Vec<[Edge; 4]> = serde_json::from_str(&json)?;
        let l = Link::from(&data);
        Ok(l)
    }
}

#[cfg(test)]
mod tests { 
    use super::*;
    use super::CrossingType::{Xp, Xn};

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
        assert_eq!(l.data[0].ctype(), Xn);
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
        let traverse = |l: &Link, (i0, j0)| { 
            let mut queue = vec![];
            l.traverse_edges((i0, j0), |i, j| queue.push((i, j)));
            queue
        };

        // single crossing
        let pd_code = [[0,1,2,3]];
        let l = Link::from(&pd_code);
        let path = traverse(&l, (0, 0));
        assert_eq!(path, [(0, 0), (0, 2)]);

        // unknot with one crossing
        let pd_code = [[0,0,1,1]];
        let l = Link::from(&pd_code);
        let path = traverse(&l, (0, 0));
        
        assert_eq!(path, [(0, 0), (0, 3), (0, 0)]); // loop

        // open arc with one crossing
        let pd_code = [[0,3,1,4],[3,2,2,1]];
        let l = Link::from(&pd_code);
        let path = traverse(&l, (0, 0));

        assert_eq!(path, [(0, 0), (1, 3), (1, 2), (0, 1), (0, 3)]); // no loop
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
        assert_eq!(comps, vec![ LinkComp::new(vec![0, 1], true)]);

        let pd_code = [[0,3,1,4],[3,2,2,1]];
        let l = Link::from(&pd_code);
        let comps = l.components();
        assert_eq!(comps, vec![ LinkComp::new(vec![0,1,2,3,4], false) ]);
    }

    #[test]
    fn link_mirror() { 
        let pd_code = [[0,0,1,1]];
        let l = Link::from(&pd_code);
        assert_eq!(l.data[0].ctype(), Xn);

        let l = l.mirror();
        assert_eq!(l.data[0].ctype(), Xp);
    }

    #[test]
    fn link_resolve() {
        let s = State::from([0, 0, 0]);
        let l = Link::from(&[[1,4,2,5],[3,6,4,1],[5,2,6,3]]) // trefoil
            .resolved_by(&s);

        let comps = l.components();
        assert_eq!(comps.len(), 3);
        assert!(comps.iter().all(|c| c.is_circle()));

        let s = State::from([1, 1, 1]);
        let l = Link::from(&[[1,4,2,5],[3,6,4,1],[5,2,6,3]]) // trefoil
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
        let l = Link::from(&pd_code);
        assert_eq!(l.crossing_num(), 2);
        assert_eq!(l.writhe(), 0);
        assert_eq!(l.components().len(), 2);
    }
}