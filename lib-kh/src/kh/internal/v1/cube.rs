use std::collections::HashMap;
use std::ops::RangeInclusive;
use itertools::Itertools;
use num_traits::Zero;
use yui_core::{Ring, RingOps, PowMod2, Sign, GetSign};
use yui_homology::{ChainComplex, Grid, Summand};
use yui_link::{Link, State, Path, Edge};

use crate::kh::{KhAlg, KhChain, KhState, KhTensor};

#[derive(Clone, Debug)]
pub struct KhCubeVertex { 
    state: State,
    circles: Vec<Path>,
    gens: Vec<KhState>
}

impl KhCubeVertex { 
    pub fn new(l: &Link, state: State, red_e: Option<Edge>, deg_shift: (isize, isize)) -> Self {
        let mut circles = l.resolve_by(&state).comps();
        circles.sort_by_key(|c| c.min_edge());
        
        let r = circles.len();

        let red_i = red_e.and_then(|e| { 
            circles.iter().position(|c| 
                c.edges().contains(&e)
            )
        });

        let gens = KhTensor::generate(r).filter_map(|label| { 
            let ok = if let Some(red_i) = red_i { 
                label[red_i].is_X()
            } else { 
                true
            };
            ok.then(|| 
                KhState::new(state, label, deg_shift)
            )
        }).collect();

        KhCubeVertex { state, circles, gens }
    }

    pub fn generators(&self) -> Vec<&KhState> { 
        self.gens.iter().collect()
    }

    pub fn circles(&self) -> &[Path] { 
        &self.circles
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

        fn diff(c1: &Vec<Path>, c2: &Vec<Path>) -> Vec<usize> { 
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

#[derive(Clone)]
pub struct KhCube<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    str: KhAlg<R>,
    dim: usize,
    vertices: HashMap<State, KhCubeVertex>,
    edges: HashMap<State, Vec<(State, KhCubeEdge)>>,
    deg_shift: (isize, isize)
}

impl<R> KhCube<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(l: &Link, h: &R, t: &R, reduce_e: Option<Edge>, deg_shift: (isize, isize)) -> Self { 
        assert!(reduce_e.is_none() || t.is_zero());

        let n = l.n_crossings();
        let str = KhAlg::new(h, t);

        let vertices: HashMap<_, _> = State::generate(n).map(|s| { 
            let v = KhCubeVertex::new(l, s, reduce_e, deg_shift);
            (s, v)
        }).collect();

        let edges: HashMap<_, _> = vertices.keys().map(|s| { 
            let edges = Self::generate_targets(s).map(|t| { 
                let v = &vertices[s];
                let w = &vertices[&t];
                (t, KhCubeEdge::edge_between(v, w))
            }).collect_vec();
            (*s, edges)
        }).collect();

        KhCube { str, dim: n, vertices, edges, deg_shift }
    }

    fn generate_targets(from: &State) -> impl Iterator<Item = State> + '_ { 
        let n = from.len();
        (0..n).filter(|&i| from[i].is_zero() ).map(move |i| { 
            from.edit(|b| b.set_1(i))
        })
    }

    pub fn str(&self) -> &KhAlg<R> {
        &self.str
    }

    pub fn dim(&self) -> usize { 
        self.dim
    }

    pub fn h_range(&self) -> RangeInclusive<isize> { 
        let i0 = self.deg_shift.0;
        let i1 = i0 + (self.dim as isize);
        i0 ..= i1
    }

    pub fn q_range(&self) -> RangeInclusive<isize> { 
        let j0 = self.deg_shift.1;
        let n = self.dim;

        let s0 = State::zeros(n);
        let v0 = self.vertex(&s0);
        let q0 = j0 - (v0.circles.len() as isize); // tensor factors are all X

        let s1 = State::ones(n);
        let v1 = self.vertex(&s1);
        let q1 = j0 + (n + v1.circles.len()) as isize; // tensor factors are all 1

        q0 ..= q1
    }

    pub fn vertex(&self, s: &State) -> &KhCubeVertex { 
        &self.vertices[s]
    }

    pub fn edge(&self, from: &State, to: &State) -> Option<&KhCubeEdge> { 
        self.edges.get(from)?.iter().find(|(t, _)| t == to).map(|(_, e)| e)
    }

    pub fn targets_from(&self, s: &State) -> impl Iterator<Item = &State> + use<'_, R> {
        self.edges[s].iter().map(|e| &e.0)
    }

    pub fn d_to(&self, x: &KhState, target: &State, signed: bool) -> KhChain<R> {
        use KhCubeEdgeTrans::*;
        
        let Some(e) = self.edge(&x.state, target) else { 
            return KhChain::zero();
        };

        let res = match e.trans { 
            Merge(ij, k) => self.merge(x, ij, k, *target),
            Split(i, jk) => self.split(x, i, jk, *target)
        };

        if signed { 
            let sign = R::from_sign(e.sign());
            res * sign
        } else { 
            res
        }
    }

    pub fn rev_d_to(&self, x: &KhState, target: &State, signed: bool) -> KhChain<R> {
        use KhCubeEdgeTrans::*;
        
        let Some(e) = self.edge(target, &x.state) else { 
            return KhChain::zero();
        };

        let res = match e.trans { 
            Merge(ij, k) => self.split(x, k, ij, *target),
            Split(i, jk) => self.merge(x, jk, i, *target)
        };

        if signed { 
            let sign = R::from_sign(e.sign());
            res * sign
        } else { 
            res
        }
    }

    fn merge(&self, x: &KhState, in_indices: (usize, usize), out_index: usize, target: State) -> KhChain<R> { 
        self.str.mul_tensor(&x.tensor, in_indices, out_index).map_keys(|y| { 
            KhState::new(target, y, x.deg_shift)
        })
    }

    fn split(&self, x: &KhState, in_index: usize, out_indices: (usize, usize), target: State) -> KhChain<R> { 
        self.str.comul_tensor(&x.tensor, in_index, out_indices).map_keys(|y| { 
            KhState::new(target, y, x.deg_shift)
        })
    }

    pub fn d(&self, x: &KhState) -> KhChain<R> {
        self.targets_from(&x.state).flat_map(|t| { 
            self.d_to(x, t, true)
        }).collect()
    }

    pub fn generators(&self, i: isize) -> Vec<&KhState> { 
        let i0 = self.deg_shift.0;
        if self.h_range().contains(&i) { 
            let i = (i - i0) as usize;
            self.states_of_weight(i).flat_map(|s| 
                self.vertex(s).generators() 
            ).sorted_by_key(|x| 
                -x.q_deg()
            ).collect()
        } else {
            vec![]
        }
    }

    fn states_of_weight(&self, k: usize) -> impl Iterator<Item = &State> { 
        self.vertices
            .iter()
            .filter_map(|(s, _)| {
                if s.weight() == k { 
                    Some(s)
                } else {
                    None
                }
            })
            .sorted_by(|s1, s2| Ord::cmp(s1, s2))
    }

    pub fn into_complex(self) -> ChainComplex<KhState, R> {
        let summands = Grid::generate(self.h_range(), |i| { 
            let gens = self.generators(i);
            Summand::from_raw_generators(gens.into_iter().cloned())
        });

        ChainComplex::new(summands, 1, move |_, z| { 
            z.apply(|x| self.d(x))
        })
    }   
}

#[cfg(test)]
mod tests { 
    use yui_core::bitseq::Bit;
    use super::*;
    
    #[test]
    fn empty() { 
        let l = Link::empty();
        let s = State::empty();
        let v = KhCubeVertex::new(&l, s, None, (0, 0));

        assert_eq!(v.state, s);
        assert_eq!(v.circles.len(), 0);
        assert_eq!(v.generators().len(), 1);
    }
    
    #[test]
    fn unknot() { 
        let l = Link::unknot();
        let s = State::empty();
        let v = KhCubeVertex::new(&l, s, None, (0, 0));

        assert_eq!(v.state, s);
        assert_eq!(v.circles.len(), 1);
        assert_eq!(v.generators().len(), 2);
    }

    #[test]
    fn unknot_red() { 
        let l = Link::unknot();
        let s = State::empty();
        let v = KhCubeVertex::new(&l, s, Some(1), (0, 0));

        assert_eq!(v.state, s);
        assert_eq!(v.circles.len(), 1);
        assert_eq!(v.generators().len(), 1);
    }

    #[test]
    fn unlink_2() {
        let l = Link::test_data("unknot_l_twist").resolve_at(0, Bit::Bit0);
        let s = State::empty();
        let v = KhCubeVertex::new(&l, s, None, (0, 0));

        assert_eq!(v.state, s);
        assert_eq!(v.circles.len(), 2);
        assert_eq!(v.generators().len(), 4);
    }

    #[test]
    fn edge_merge() { 
        let l = Link::test_data("unknot_l_twist");
        let s = State::from([0]);
        let t = State::from([1]);
        let v = KhCubeVertex::new(&l, s, None, (0, 0));
        let w = KhCubeVertex::new(&l, t, None, (0, 0));
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
        let l = Link::test_data("unknot_r_twist");
        let s = State::from([0]);
        let t = State::from([1]);
        let v = KhCubeVertex::new(&l, s, None, (0, 0));
        let w = KhCubeVertex::new(&l, t, None, (0, 0));
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
        let cube = KhCube::<i32>::new(&l, &0, &0, None, (0, 0));

        assert_eq!(cube.dim, 0);
        assert_eq!(cube.vertices.len(), 1);
        assert_eq!(cube.edges.len(), 1);

        let s = State::empty();
        let v = cube.vertex(&s);

        assert_eq!(v.circles.len(), 0);
        assert_eq!(v.generators().len(), 1);

        assert_eq!(cube.targets_from(&s).next(), None);
    }

    #[test]
    fn cube_unknot() { 
        let l = Link::unknot();
        let cube = KhCube::<i32>::new(&l, &0, &0, None, (0, 0));

        assert_eq!(cube.dim, 0);
        assert_eq!(cube.vertices.len(), 1);

        let s = State::empty();
        let v = cube.vertex(&s);

        assert_eq!(v.circles.len(), 1);
        assert_eq!(v.generators().len(), 2);

        assert_eq!(cube.targets_from(&s).next(), None);
    }

    #[test]
    fn cube_twist_unknot() { 
        let l = Link::test_data("unknot_l_twist");
        let cube = KhCube::<i32>::new(&l, &0, &0, None, (0, 0));

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
        let l = Link::test_data("L2a1");
        let cube = KhCube::<i32>::new(&l, &0, &0, None, (0, 0));

        assert_eq!(cube.dim, 2);
        assert_eq!(cube.vertices.len(), 4);
   }

   #[test]
   fn cube_trefoil() { 
       let l = Link::test_data("3_1");
       let cube = KhCube::<i32>::new(&l, &0, &0, None, (0, 0));

       assert_eq!(cube.dim, 3);
       assert_eq!(cube.vertices.len(), 8);
    }
}