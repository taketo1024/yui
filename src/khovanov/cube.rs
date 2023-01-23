use std::collections::HashMap;
use std::ops::RangeInclusive;
use itertools::Itertools;
use num_traits::Pow;
use crate::links::links::{Link, State, Component, Resolution, Edge};
use crate::math::traits::{Ring, RingOps, PowMod2};
use crate::math::sign::Sign;
use super::algebra::{KhAlgStr, KhEnhState};

#[derive(Debug)]
pub struct KhCubeVertex { 
    state: State,
    circles: Vec<Component>
}

impl KhCubeVertex { 
    pub fn new(l: &Link, state: State) -> Self {
        let circles = l.resolved_by(&state).components();
        KhCubeVertex { state, circles }
    }

    pub fn generators(&self) -> Vec<KhEnhState> { 
        self.collect_generators(None)
    }

    pub fn reduced_generators(&self, red_e: &Edge) -> Vec<KhEnhState> { 
        self.collect_generators(Some(red_e))
    }

    fn collect_generators(&self, red_e: Option<&Edge>) -> Vec<KhEnhState> { 
        use super::algebra::KhAlgGen::{X, I};
        let s = &self.state;
        let r = self.circles.len();

        let red_i = red_e.map(|red_e| 
            self.circles.iter().find_position(|c| 
                c.edges().contains(red_e)
            ).unwrap().0 // must exist
        );

        return (0..2.pow(r)).filter_map(|b| { 
            if let Some(red_i) = red_i { 
                if (b >> red_i) & 1 == 1 { 
                    return None
                }
            }
            
            let state = s.clone();
            let labels = (0..r).map(|i| {
                if (b >> i) & 1 == 1 { I } else { X }
            }).collect();
            let x = KhEnhState::new( state, labels );

            Some(x)
        }).collect();
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KhCubeEdgeTrans { 
    Merge((usize, usize), usize),
    Split(usize, (usize, usize))
}

#[derive(Debug, Clone)]
pub struct KhCubeEdge { 
    trans: KhCubeEdgeTrans,
    sign: Sign
}

impl KhCubeEdge { 
    pub fn sign(&self) -> Sign { 
        self.sign
    }

    pub fn trans(&self) -> &KhCubeEdgeTrans { 
        &self.trans
    }
    
    fn edge_between(from: &KhCubeVertex, to: &KhCubeVertex) -> Self { 
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

        let sign = Self::sign_between(&from.state, &to.state);
        KhCubeEdge { trans, sign }
    }

    fn sign_between(from: &State, to: &State) -> Sign { 
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
    edges: HashMap<State, Vec<(State, KhCubeEdge)>>
}

impl<R> KhCube<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(l: &Link) -> Self { 
        Self::new_ht(l, R::zero(), R::zero())
    }

    pub fn new_ht(l: &Link, h: R, t: R) -> Self { 
        let str = KhAlgStr::new(h, t);
        Self::_new(l, str)
    }
    
    fn _new(l: &Link, str: KhAlgStr<R>) -> Self { 
        let dim = l.crossing_num() as usize;
        let m = 2.pow(dim) as usize;
        let vertices: HashMap<_, _> = (0..m).map(|i| { 
            let s = State::from_bseq(i, dim);
            let v = KhCubeVertex::new(&l, s.clone());
            (s, v)
        }).collect();

        let edges: HashMap<_, _> = (0..m).map(|i| { 
            let s = State::from_bseq(i, dim);
            let edges = s.targets().into_iter().map(|t| { 
                let v = &vertices[&s];
                let w = &vertices[&t];
                (t, KhCubeEdge::edge_between(v, w))
            }).collect_vec();
            (s, edges)
        }).collect();

        KhCube { str, dim, vertices, edges }
    }

    pub fn structure(&self) -> &KhAlgStr<R> {
        &self.str
    }

    pub fn dim(&self) -> usize { 
        self.dim
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        0 ..= (self.dim as isize)
    }

    pub fn q_range(&self) -> RangeInclusive<isize> { 
        let n = self.dim;

        let s0 = State::from_bseq(0, n);
        let v0 = self.vertex(&s0);
        let q0 = -(v0.circles.len() as isize); // tensor factors are all X

        let s1 = State::from_bseq((2.pow(n) - 1) as usize, n);
        let v1 = self.vertex(&s1);
        let q1 = (n + v1.circles.len()) as isize; // tensor factors are all 1

        q0 ..= q1
    }

    pub fn generators(&self, i: isize) -> Vec<KhEnhState> { 
        self.collect_generators(i, None)
    }

    pub fn reduced_generators(&self, i: isize, red_e: &Edge) -> Vec<KhEnhState> { 
        self.collect_generators(i, Some(red_e))
    }

    fn collect_generators(&self, i: isize, red_e: Option<&Edge>) -> Vec<KhEnhState> { 
        if self.h_range().contains(&i) { 
            let i = i as usize;
            self.vertices_of_weight(i).into_iter().flat_map(|v| 
                if let Some(red_e) = red_e { 
                    v.reduced_generators(red_e)
                } else { 
                    v.generators() 
                }
            ).collect()
        } else {
            vec![]
        }
    }

    pub fn differentiate(&self, x: &KhEnhState) -> Vec<(KhEnhState, R)> {
        let edges = self.edges_from(x.state());
        edges.iter().flat_map(|(t, e)| { 
            self.apply(x, t, e)
        }).collect_vec()
    }

    pub fn vertex(&self, s: &State) -> &KhCubeVertex { 
        &self.vertices[s]
    }

    pub fn edge(&self, from: &State, to: &State) -> Option<&KhCubeEdge> { 
        self.edges[from].iter().find(|(t, _)| t == to).map(|(_, e)| e)
    }

    fn vertices_of_weight(&self, k: usize) -> Vec<&KhCubeVertex> { 
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

    fn edges_from(&self, s: &State) -> &Vec<(State, KhCubeEdge)> {
        &self.edges[s]
    }

    fn apply(&self, x: &KhEnhState, to: &State, e: &KhCubeEdge) -> Vec<(KhEnhState, R)> {
        use KhCubeEdgeTrans::*;
        
        let sign = R::from(e.sign);

        match e.trans { 
            Merge((i, j), k) => {
                let (x_i, x_j) = (x.label_at(i), x.label_at(j));
                self.str.prod(x_i, x_j).into_iter().map(|(y_k, a)| { 
                    let mut label = x.label().clone();
                    label.remove(j);
                    label.remove(i);
                    label.insert(k, y_k);

                    let t = to.clone();
                    let y = KhEnhState::new(t, label);
                    let r = &sign * &a;
                    (y, r)
                }).collect_vec()
            },
            Split(i, (j, k)) => {
                let x_i = x.label_at(i);
                self.str.coprod(x_i).into_iter().map(|(y_j, y_k, a)| { 
                    let mut label = x.label().clone();
                    label.remove(i);
                    label.insert(j, y_j);
                    label.insert(k, y_k);

                    let t = to.clone();
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
        let v = KhCubeVertex::new(&l, s.clone());

        assert_eq!(v.state, s);
        assert_eq!(v.circles.len(), 0);
        assert_eq!(v.generators().len(), 1);
    }
    
    #[test]
    fn unknot() { 
        let l = Link::unknot();
        let s = State::empty();
        let v = KhCubeVertex::new(&l, s.clone());

        assert_eq!(v.state, s);
        assert_eq!(v.circles.len(), 1);
        assert_eq!(v.generators().len(), 2);
    }

    #[test]
    fn unlink_2() {
        let l = Link::from(&[[0, 0, 1, 1]]).resolved_at(0, Res0);
        let s = State::empty();
        let v = KhCubeVertex::new(&l, s.clone());

        assert_eq!(v.state, s);
        assert_eq!(v.circles.len(), 2);
        assert_eq!(v.generators().len(), 4);
    }

    #[test]
    fn edge_merge() { 
        let l = Link::from(&[[0, 0, 1, 1]]);
        let s = State::from([0]);
        let t = State::from([1]);
        let v = KhCubeVertex::new(&l, s);
        let w = KhCubeVertex::new(&l, t);
        let e = KhCubeEdge::edge_between(&v, &w);

        assert!(e.sign.is_positive());

        let KhCubeEdgeTrans::Merge(from, to) = e.trans else { 
            panic!()
        };
        assert_eq!(from, (0, 1));
        assert_eq!(to, 0);
    }

    #[test]
    fn edge_split() { 
        let l = Link::from(&[[0, 1, 1, 0]]);
        let s = State::from([0]);
        let t = State::from([1]);
        let v = KhCubeVertex::new(&l, s);
        let w = KhCubeVertex::new(&l, t);
        let e = KhCubeEdge::edge_between(&v, &w);

        assert!(e.sign.is_positive());

        let KhCubeEdgeTrans::Split(from, to) = e.trans else { 
            panic!()
        };
        assert_eq!(from, 0);
        assert_eq!(to, (0, 1));
    }

    #[test]
    fn edge_sign() { 
        let s = State::from([0, 0, 0]);
        let t = State::from([1, 0, 0]);
        let e = KhCubeEdge::sign_between(&s, &t);
        assert!(e.is_positive());

        let s = State::from([1, 0, 0]);
        let t = State::from([1, 1, 0]);
        let e = KhCubeEdge::sign_between(&s, &t);
        assert!(e.is_negative());

        let s = State::from([1, 1, 0]);
        let t = State::from([1, 1, 1]);
        let e = KhCubeEdge::sign_between(&s, &t);
        assert!(e.is_positive());

        let s = State::from([0, 1, 0]);
        let t = State::from([0, 1, 1]);
        let e = KhCubeEdge::sign_between(&s, &t);
        assert!(e.is_negative());
    }

    #[test]
    fn cube_empty() { 
        let l = Link::empty();
        let cube = KhCube::<i32>::new(&l);

        assert_eq!(cube.dim, 0);
        assert_eq!(cube.vertices.len(), 1);
        assert_eq!(cube.edges.len(), 1);

        let s = State::empty();
        let v = cube.vertex(&s);

        assert_eq!(v.circles.len(), 0);
        assert_eq!(v.generators().len(), 1);

        assert!(cube.edges_from(&s).is_empty());
    }

    #[test]
    fn cube_unknot() { 
        let l = Link::unknot();
        let cube = KhCube::<i32>::new(&l);

        assert_eq!(cube.dim, 0);
        assert_eq!(cube.vertices.len(), 1);

        let s = State::empty();
        let v = cube.vertex(&s);

        assert_eq!(v.circles.len(), 1);
        assert_eq!(v.generators().len(), 2);

        assert!(cube.edges_from(&s).is_empty());
    }

    #[test]
    fn cube_twist_unknot() { 
        let l = Link::from(&[[0, 0, 1, 1]]);
        let cube = KhCube::<i32>::new(&l);

        assert_eq!(cube.dim, 1);
        assert_eq!(cube.vertices.len(), 2);

        let s0 = State::from([0]);
        let s1 = State::from([1]);

        let v0 = cube.vertex(&s0);
        let v1 = cube.vertex(&s1);

        assert_eq!(v0.circles.len(), 2);
        assert_eq!(v0.generators().len(), 4);

        assert_eq!(v1.circles.len(), 1);
        assert_eq!(v1.generators().len(), 2);

        let Some(e) = cube.edge(&s0, &s1) else { panic!() };

        assert_eq!(e.sign(), Sign::Pos);
        assert_eq!(e.trans(), &KhCubeEdgeTrans::Merge((0, 1), 0));
    }

    #[test]
    fn cube_hopf_link() { 
        let l = Link::hopf_link();
        let cube = KhCube::<i32>::new(&l);

        assert_eq!(cube.dim, 2);
        assert_eq!(cube.vertices.len(), 4);
   }

   #[test]
   fn cube_trefoil() { 
       let l = Link::trefoil();
       let cube = KhCube::<i32>::new(&l);

       assert_eq!(cube.dim, 3);
       assert_eq!(cube.vertices.len(), 8);
  }
}