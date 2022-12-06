use std::collections::HashMap;
use itertools::Itertools;
use num_traits::Pow;
use crate::links::links::{Link, State, Component, Resolution};
use crate::math::traits::{Ring, RingOps, PowMod2};
use crate::math::sign::Sign;
use super::algebra::{KhAlgGen, KhAlgStr};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct KhEnhState { 
    state: State,
    label: Vec<KhAlgGen>
}

impl KhEnhState {
    pub fn new(state: State, label: Vec<KhAlgGen>) -> KhEnhState { 
        KhEnhState { state, label }
    }
}

#[derive(Debug)]
struct KhCubeVertex { 
    state: State,
    circles: Vec<Component>,
    generators: Vec<KhEnhState>
}

impl KhCubeVertex { 
    pub fn new(l: &Link, s: &State) -> Self {
        use super::algebra::KhAlgGen::{X, I};
        let circles = l.clone().resolve(&s).components();
        let r = circles.len();
        let generators = (0..2.pow(r)).map(|mut i| { 
            let labels = (0..r).map(|_j| { 
                let x = if i & 1 == 1 { I } else { X };
                i >>= 1;
                x
            }).collect_vec();
            KhEnhState::new( s.clone(), labels )
        }).collect_vec();

        KhCubeVertex { state: s.clone(), circles, generators }
    }

    pub fn generators(&self) -> Vec<&KhEnhState> { 
        self.generators.iter().collect()
    }
}

#[derive(Debug)]
enum KhCubeEdgeTrans { 
    Merge((usize, usize), usize),
    Split(usize, (usize, usize))
}

#[derive(Debug)]
struct KhCubeEdge { 
    from: State,
    to: State,
    trans: KhCubeEdgeTrans,
    sign: Sign
}

impl KhCubeEdge { 
    fn between(from: &KhCubeVertex, to: &KhCubeVertex) -> Self { 
        debug_assert!(from.state.weight() + 1 == to.state.weight());

        fn diff(c1: &Vec<Component>, c2: &Vec<Component>) -> Vec<usize> { 
            let (n1, n2) = (c1.len(), c2.len());
            assert!(n1 == n2 + 1 || n1 + 1 == n2);
            c1.iter().enumerate().filter_map(|(i, c)| { 
                if !c2.contains(c) { Some(i) } else { None }
            }).collect_vec()
        }

        let c_from = diff(&from.circles, &to.circles);
        let c_to   = diff(&to.circles, &from.circles);

        let trans = match (c_from.len(), c_to.len()) { 
            (2, 1) => KhCubeEdgeTrans::Merge((c_from[0], c_from[1]), c_to[0]),
            (1, 2) => KhCubeEdgeTrans::Split(c_from[0], (c_to[0], c_to[1])),
            _ => panic!()
        };

        let sign = Self::sign(&from.state, &to.state);
        KhCubeEdge { from: from.state.clone(), to: to.state.clone(), trans, sign }
    }

    fn sign(from: &State, to: &State) -> Sign { 
        use Resolution::Res1;
        debug_assert_eq!(from.len(), to.len());
        debug_assert_eq!(from.weight() + 1, to.weight());

        let n = from.len();
        let i = (0..n).find(|&i| from[i] != to[i]).unwrap();
        let k = (0..i).filter(|&j| from[j] == Res1).count() as u32;

        Sign::from( (-1).pow_mod2(k) )
    }
}

pub struct KhCube<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    str: KhAlgStr<R>,
    dim: usize,
    vertices: HashMap<State, KhCubeVertex>,
    edges: HashMap<State, Vec<KhCubeEdge>>
}

impl<R> KhCube<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(l: &Link, str: KhAlgStr<R>) -> Self { 
        let dim = l.crossing_num() as usize;
        let m = 2.pow(dim) as usize;

        let vertices: HashMap<_, _> = (0..m).map(|i| { 
            let s = State::from_bseq(i, dim);
            let v = KhCubeVertex::new(&l, &s);
            (s, v)
        }).collect();

        let edges: HashMap<_, _> = (0..m).map(|i| { 
            let s = State::from_bseq(i, dim);
            let edges = s.targets().into_iter().map(|t| { 
                let v = &vertices[&s];
                let w = &vertices[&t];
                KhCubeEdge::between(v, w)
            }).collect_vec();
            (s, edges)
        }).collect();

        KhCube { str, dim, vertices, edges }
    }

    pub fn dim(&self) -> usize { 
        self.dim
    }

    pub fn generators(&self, k: usize) -> Vec<&KhEnhState> { 
        self.vertices(k).into_iter().flat_map(|v| v.generators()).collect_vec()
    }

    pub fn differentiate(&self, x: &KhEnhState) -> Vec<(KhEnhState, R)> {
        let edges = self.edges_from(&x.state);
        edges.iter().flat_map(|e| { 
            self.apply(x, e)
        }).collect_vec()
    }

    fn vertex(&self, s: &State) -> &KhCubeVertex { 
        &self.vertices[s]
    }

    fn vertices(&self, k: usize) -> Vec<&KhCubeVertex> { 
        let k = k as u32;
        self.vertices
            .iter()
            .sorted_by(|(s1, _), (s2, _)| Ord::cmp(s1, s2))
            .filter_map(|(s, v)| {
                if s.weight() == k { 
                    Some(v)
                } else {
                    None
                }
            })
            .collect_vec()
    }

    fn edges_from(&self, s: &State) -> &Vec<KhCubeEdge> {
        &self.edges[s]
    }

    fn apply(&self, x: &KhEnhState, e: &KhCubeEdge) -> Vec<(KhEnhState, R)> {
        use KhCubeEdgeTrans::*;
        
        assert_eq!(x.state, e.from);
        let sign = R::from(e.sign);

        match e.trans { 
            Merge((i, j), k) => {
                let (x_i, x_j) = (x.label[i], x.label[j]);
                self.str.prod(x_i, x_j).into_iter().map(|(y_k, a)| { 
                    let t = e.to.clone();
                    let mut label = x.label.clone();
                    label.remove(j);
                    label.remove(i);
                    label.insert(k, y_k);

                    let y = KhEnhState::new(t, label);
                    let r = &sign * &a;
                    (y, r)
                }).collect_vec()
            },
            Split(i, (j, k)) => {
                let x_i = x.label[i];
                self.str.coprod(x_i).into_iter().map(|(y_j, y_k, a)| { 
                    let t = e.to.clone();
                    let mut label = x.label.clone();
                    label.remove(i);
                    label.insert(j, y_j);
                    label.insert(k, y_k);

                    let y = KhEnhState::new(t, label);
                    let r = &sign * &a;
                    (y, r)
                }).collect_vec()
            }
        }
    }
}

#[cfg(test)]
mod tests { 
    use crate::links::links::Resolution::*;
    use super::*;
    
    #[test]
    fn empty() { 
        let l = Link::empty();
        let s = State::empty();
        let v = KhCubeVertex::new(&l, &s);

        assert_eq!(v.state, s);
        assert_eq!(v.circles.len(), 0);
        assert_eq!(v.generators.len(), 1);
    }
    
    #[test]
    fn unknot() { 
        let l = Link::unknot();
        let s = State::empty();
        let v = KhCubeVertex::new(&l, &s);

        assert_eq!(v.state, s);
        assert_eq!(v.circles.len(), 1);
        assert_eq!(v.generators.len(), 2);
    }

    #[test]
    fn unlink_2() {
        let l = Link::from([[0, 0, 1, 1]]).resolve_at(0, Res0);
        let s = State::empty();
        let v = KhCubeVertex::new(&l, &s);

        assert_eq!(v.state, s);
        assert_eq!(v.circles.len(), 2);
        assert_eq!(v.generators.len(), 4);
    }

    #[test]
    fn edge_merge() { 
        let l = Link::from([[0, 0, 1, 1]]);
        let s = State::from(vec![0]);
        let t = State::from(vec![1]);
        let v = KhCubeVertex::new(&l, &s);
        let w = KhCubeVertex::new(&l, &t);
        let e = KhCubeEdge::between(&v, &w);

        assert_eq!(e.to, t);
        assert!(e.sign.is_positive());

        let KhCubeEdgeTrans::Merge(from, to) = e.trans else { 
            panic!()
        };
        assert_eq!(from, (0, 1));
        assert_eq!(to, 0);
    }

    #[test]
    fn edge_split() { 
        let l = Link::from([[0, 1, 1, 0]]);
        let s = State::from(vec![0]);
        let t = State::from(vec![1]);
        let v = KhCubeVertex::new(&l, &s);
        let w = KhCubeVertex::new(&l, &t);
        let e = KhCubeEdge::between(&v, &w);

        assert_eq!(e.to, t);
        assert!(e.sign.is_positive());

        let KhCubeEdgeTrans::Split(from, to) = e.trans else { 
            panic!()
        };
        assert_eq!(from, 0);
        assert_eq!(to, (0, 1));
    }

    #[test]
    fn edge_sign() { 
        let s = State::from(vec![0, 0, 0]);
        let t = State::from(vec![1, 0, 0]);
        let e = KhCubeEdge::sign(&s, &t);
        assert!(e.is_positive());

        let s = State::from(vec![1, 0, 0]);
        let t = State::from(vec![1, 1, 0]);
        let e = KhCubeEdge::sign(&s, &t);
        assert!(e.is_negative());

        let s = State::from(vec![1, 1, 0]);
        let t = State::from(vec![1, 1, 1]);
        let e = KhCubeEdge::sign(&s, &t);
        assert!(e.is_positive());

        let s = State::from(vec![0, 1, 0]);
        let t = State::from(vec![0, 1, 1]);
        let e = KhCubeEdge::sign(&s, &t);
        assert!(e.is_negative());
    }

    #[test]
    fn cube_empty() { 
        let l = Link::empty();
        let str = KhAlgStr::new(0, 0);
        let cube = KhCube::new(&l, str);

        assert_eq!(cube.dim, 0);
        assert_eq!(cube.vertices.len(), 1);
        assert_eq!(cube.edges.len(), 1);

        let s = State::empty();
        let v = cube.vertex(&s);

        assert_eq!(v.circles.len(), 0);
        assert_eq!(v.generators.len(), 1);

        assert!(cube.edges_from(&s).is_empty());
    }

    #[test]
    fn cube_unknot() { 
        let l = Link::unknot();
        let str = KhAlgStr::new(0, 0);
        let cube = KhCube::new(&l, str);

        assert_eq!(cube.dim, 0);
        assert_eq!(cube.vertices.len(), 1);

        let s = State::empty();
        let v = cube.vertex(&s);

        assert_eq!(v.circles.len(), 1);
        assert_eq!(v.generators.len(), 2);

        assert!(cube.edges_from(&s).is_empty());
    }

    #[test]
    fn cube_twist_unknot() { 
        let l = Link::from([[0, 0, 1, 1]]);
        let str = KhAlgStr::new(0, 0);
        let cube = KhCube::new(&l, str);

        assert_eq!(cube.dim, 1);
        assert_eq!(cube.vertices.len(), 2);

        let s0 = State::from(vec![0]);
        let s1 = State::from(vec![1]);

        let v0 = cube.vertex(&s0);
        let v1 = cube.vertex(&s1);

        assert_eq!(v0.circles.len(), 2);
        assert_eq!(v0.generators.len(), 4);

        assert_eq!(v1.circles.len(), 1);
        assert_eq!(v1.generators.len(), 2);

        assert_eq!(cube.edges_from(&s0).len(), 1);
        assert_eq!(cube.edges_from(&s1).len(), 0);
    }

    #[test]
    fn cube_hopf_link() { 
        let l = Link::hopf_link();
        let str = KhAlgStr::new(0, 0);
        let cube = KhCube::new(&l, str);

        assert_eq!(cube.dim, 2);
        assert_eq!(cube.vertices.len(), 4);
   }

   #[test]
   fn cube_trefoil() { 
       let l = Link::trefoil();
       let str = KhAlgStr::new(0, 0);
       let cube = KhCube::new(&l, str);

       assert_eq!(cube.dim, 3);
       assert_eq!(cube.vertices.len(), 8);
  }
}