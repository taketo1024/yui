use std::collections::HashMap;
use std::ops::RangeInclusive;
use itertools::Itertools;
use yui_core::{Ring, RingOps, PowMod2, Sign, GetSign};
use yui_homology::XChainComplex;
use yui_link::{Link, State, LinkComp, Edge};

use crate::{KhAlgStr, KhLabel, KhEnhState};

#[derive(Debug)]
pub struct KhCubeVertex { 
    state: State,
    circles: Vec<LinkComp>,
    gens: Vec<KhEnhState>
}

impl KhCubeVertex { 
    pub fn new(l: &Link, state: State) -> Self {
        let circles = l.resolved_by(&state).components();
        let r = circles.len();
        let gens = KhLabel::generate(r).into_iter().map(|label| { 
            KhEnhState::new( state, label )
        }).collect();
        KhCubeVertex { state, circles, gens }
    }

    pub fn generators(&self) -> Vec<&KhEnhState> { 
        self.gens.iter().collect()
    }

    pub fn reduced_generators(&self, red_e: &Edge) -> Vec<&KhEnhState> { 
        let red_i = self.circles.iter().position(|c| 
            c.edges().contains(red_e)
        ).unwrap(); // must exist

        self.generators().into_iter().filter(|x| { 
            x.label[red_i].is_X()
        }).collect()
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

        fn diff(c1: &Vec<LinkComp>, c2: &Vec<LinkComp>) -> Vec<usize> { 
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
        debug_assert_eq!(from.len(), to.len());
        debug_assert_eq!(from.weight() + 1, to.weight());

        let n = from.len();
        let i = (0..n).find(|&i| from[i] != to[i]).unwrap();
        let k = (0..i).filter(|&j| from[j].is_one()).count() as u32;

        (-1).pow_mod2(k).sign()
    }
}

pub struct KhCube<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    str: KhAlgStr<R>,
    dim: usize,
    vertices: HashMap<State, KhCubeVertex>,
    edges: HashMap<State, Vec<(State, KhCubeEdge)>>,
    base_pt: Option<Edge>
}

impl<R> KhCube<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(l: &Link, h: &R, t: &R) -> Self { 
        let str = KhAlgStr::new(h, t);

        let n = l.crossing_num();
        let vertices: HashMap<_, _> = State::generate(n).into_iter().map(|s| { 
            let v = KhCubeVertex::new(&l, s.clone());
            (s, v)
        }).collect();

        let edges: HashMap<_, _> = vertices.keys().map(|s| { 
            let edges = Self::targets(s).into_iter().map(|t| { 
                let v = &vertices[s];
                let w = &vertices[&t];
                (t, KhCubeEdge::edge_between(v, w))
            }).collect_vec();
            (s.clone(), edges)
        }).collect();

        let base_pt = l.first_edge();

        KhCube { str, dim: n, vertices, edges, base_pt }
    }

    fn targets(from: &State) -> Vec<State> { 
        let n = from.len();
        (0..n).filter(|&i| from[i].is_zero() ).map(|i| { 
            from.edit(|b| b.set_1(i))
        }).collect()
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

        let s0 = State::zeros(n);
        let v0 = self.vertex(&s0);
        let q0 = -(v0.circles.len() as isize); // tensor factors are all X

        let s1 = State::ones(n);
        let v1 = self.vertex(&s1);
        let q1 = (n + v1.circles.len()) as isize; // tensor factors are all 1

        q0 ..= q1
    }

    pub fn generators(&self, i: isize) -> Vec<&KhEnhState> { 
        self.collect_generators(i, None)
    }

    pub fn reduced_generators(&self, i: isize, red_e: &Edge) -> Vec<&KhEnhState> { 
        self.collect_generators(i, Some(red_e))
    }

    fn collect_generators(&self, i: isize, red_e: Option<&Edge>) -> Vec<&KhEnhState> { 
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
        let edges = self.edges_from(&x.state);
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
        
        let sign = R::from_sign(e.sign);

        match e.trans { 
            Merge((i, j), k) => {
                let (x_i, x_j) = (x.label[i], x.label[j]);
                self.str.prod(x_i, x_j).into_iter().map(|(y_k, a)| { 
                    let mut label = x.label.clone();
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
                let x_i = x.label[i];
                self.str.coprod(x_i).into_iter().map(|(y_j, y_k, a)| { 
                    let mut label = x.label.clone();
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

    pub fn as_complex(self, i0: isize, reduced: bool) -> XChainComplex<KhEnhState, R> {
        let range = self.h_range();
        let range = (range.start() + i0) ..= (range.end() + i0);
        
        XChainComplex::new(range, 1, 
            |i| {
                let i = i - i0;
                let gens = if reduced {
                    let e = self.base_pt.unwrap();
                    self.reduced_generators(i, &e)
                } else { 
                    self.generators(i) 
                };
                gens.into_iter().cloned().collect()
            },
            |_i, x| { 
                self.differentiate(x)
            }
        )
    }   
}

#[cfg(test)]
mod tests { 
    use yui_utils::bitseq::Bit;
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
        let l = Link::from_pd_code([[0, 0, 1, 1]]).resolved_at(0, Bit::Bit0);
        let s = State::empty();
        let v = KhCubeVertex::new(&l, s.clone());

        assert_eq!(v.state, s);
        assert_eq!(v.circles.len(), 2);
        assert_eq!(v.generators().len(), 4);
    }

    #[test]
    fn edge_merge() { 
        let l = Link::from_pd_code([[0, 0, 1, 1]]);
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
        let l = Link::from_pd_code([[0, 1, 1, 0]]);
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
        let cube = KhCube::<i32>::new(&l, &0, &0);

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
        let cube = KhCube::<i32>::new(&l, &0, &0);

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
        let l = Link::from_pd_code([[0, 0, 1, 1]]);
        let cube = KhCube::<i32>::new(&l, &0, &0);

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
        let cube = KhCube::<i32>::new(&l, &0, &0);

        assert_eq!(cube.dim, 2);
        assert_eq!(cube.vertices.len(), 4);
   }

   #[test]
   fn cube_trefoil() { 
       let l = Link::trefoil();
       let cube = KhCube::<i32>::new(&l, &0, &0);

       assert_eq!(cube.dim, 3);
       assert_eq!(cube.vertices.len(), 8);
    }
}