use itertools::Itertools;
use yui::bitseq::Bit;
use yui::lc::Lc;
use yui::{CloneAnd, Ring, RingOps, Sign};
use yui_homology::ChainMap;
use yui_link::{Link, Path, State};
use num_traits::Zero;

use crate::kh::internal::v1::cube::KhCube;
use crate::kh::{KhAlg, KhChain, KhChainGen, KhComplex, KhGen, KhTensor};

type KhChainMap<R> = ChainMap<isize, KhChainGen, KhChainGen, R>;

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn cc_pair(l: &Link, h: &R, t: &R, reduced: bool, i: usize) -> (KhComplex<R>, KhComplex<R>) {
        assert!(l.node(i).is_crossing());
        
        let l2 = l.crossing_change(i);
        let c1 = KhComplex::new_no_simplify(l, h, t, reduced);
        let c2 = KhComplex::new_no_simplify(&l2, h, t, reduced);
        (c1, c2)
    }

    pub fn cc_map0(c1: &KhComplex<R>, c2: &KhComplex<R>, i: usize) -> KhChainMap<R> {
        let deg = c2.deg_shift().0 - c1.deg_shift().0 + 1;
        let c2_deg_shift = c2.deg_shift();

        ChainMap::new(c1.inner(), c2.inner(), deg, move |_, z| { 
            z.apply(|x| {
                if !x.state[i].is_zero() { 
                    return KhChain::zero();
                }

                let e = Sign::from_parity( count_1s(&x.state, i) );
                let y = x.clone_and(|y| { 
                    y.state.set_1(i);
                    y.deg_shift = c2_deg_shift;
                });

                KhChain::from(y) * R::from_sign(e)
            })
        })
    }

    pub fn cc_map1(c1: &KhComplex<R>, c2: &KhComplex<R>, l: &Link, i: usize) -> KhChainMap<R> {
        assert!(l.node(i).is_crossing());

        let deg = c2.deg_shift().0 - c1.deg_shift().0 - 1;
        let c2_deg_shift = c2.deg_shift();

        let alg = c1.str().clone();
        let (a0, a1) = l.node(i).resolved(Bit::Bit0).arcs();

        // TODO We don't want to reproduce the cube. 
        let (h, t) = c1.str().ht();
        let red_e = if c1.is_reduced() { l.first_edge() } else { None };
        let cube = KhCube::new(l, h, t, red_e, c1.deg_shift());

        ChainMap::new(c1.inner(), c2.inner(), deg, move |_, z| { 
            z.apply(|x| {
                if !x.state[i].is_one() { 
                    return KhChain::zero();
                }

                let circles = cube.vertex(&x.state).circles();
                let (k0, k1) = (circle_index(circles, &a0), circle_index(circles, &a1));

                if k0 == k1 { 
                    return Lc::zero();
                }

                let mut s = x.state;
                s.set_0(i);

                let e = Sign::from_parity( count_1s(&x.state, i) );
                let t = apply_f1(&alg, &x.tensor, k0, k1) * R::from_sign(e);
                
                t.into_map_gens(|y| { 
                    KhChainGen::new(s, y, c2_deg_shift)
                })
            })
        })
    }
}

fn count_1s(s: &State, i: usize) -> u32 { 
    s.iter().enumerate().filter(|(j, b)| j > &i && b.is_one()).count() as u32
}

fn circle_index(circles: &[Path], arc: &Path) -> usize { 
    circles.iter().find_position(|c| c.edges().contains(&arc.min_edge())).unwrap().0
}

fn apply_f1<R>(alg: &KhAlg<R>, x: &KhTensor, i0: usize, i1: usize) -> Lc<KhTensor, R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    use KhGen::X;

    let w0 = x.apply_at(i0, |x0| 
        alg.mul(*x0, X)  // multiply X at i0
    );
    let w1 = x.apply_at(i1, |x1| 
        alg.mul(*x1, X)  // multiply X at i0
    );

    w0 - w1 
}

#[cfg(test)]
mod tests {
    use yui::poly::Poly2;
    use yui::{EucRing, EucRingOps};
    use yui_homology::DisplaySeq;
    use yui_link::Link;

    use crate::kh::ext::cc::KhChainMap;
    use crate::kh::KhComplex;
 
    #[test]
    fn test_cc0_pos_to_neg() { 
        let i = 0;
        let l = Link::load("5_1").unwrap().crossing_change(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, false, i);
        let f = KhComplex::cc_map0(&c1, &c2, i);

        assert_eq!(f.deg(), 0);

        f.check_all(c1.inner(), c2.inner());
    }

    #[test]
    fn test_cc0_neg_to_pos() { 
        let i = 0;
        let l = Link::load("5_1").unwrap().mirror().crossing_change(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, false, i);
        let f = KhComplex::cc_map0(&c1, &c2, i);

        assert_eq!(f.deg(), 2);

        f.check_all(c1.inner(), c2.inner());
    }

    #[test]
    fn test_cc1_pos_to_neg() { 
        let i = 0;
        let l = Link::load("5_1").unwrap().crossing_change(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, false, i);
        let f = KhComplex::cc_map1(&c1, &c2, &l, i);

        assert_eq!(f.deg(), -2);

        f.check_all(c1.inner(), c2.inner());
    }

    #[test]
    fn test_cc1_neg_to_pos() { 
        let i = 0;
        let l = Link::load("5_1").unwrap().mirror().crossing_change(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, false, i);
        let f = KhComplex::cc_map1(&c1, &c2, &l, i);

        assert_eq!(f.deg(), 0);

        f.check_all(c1.inner(), c2.inner());
    }

    #[test]
    fn test_red_cc0_pos_to_neg() { 
        let i = 0;
        let l = Link::load("5_1").unwrap().crossing_change(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, true, i);
        let f = KhComplex::cc_map0(&c1, &c2, i);

        assert_eq!(f.deg(), 0);

        f.check_all(c1.inner(), c2.inner());
    }

    #[test]
    fn test_red_cc0_neg_to_pos() { 
        let i = 0;
        let l = Link::load("5_1").unwrap().mirror().crossing_change(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, true, i);
        let f = KhComplex::cc_map0(&c1, &c2, i);

        assert_eq!(f.deg(), 2);

        f.check_all(c1.inner(), c2.inner());
    }

    #[test]
    fn test_red_cc1_pos_to_neg() { 
        let i = 0;
        let l = Link::load("5_1").unwrap().crossing_change(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, true, i);
        let f = KhComplex::cc_map1(&c1, &c2, &l, i);

        assert_eq!(f.deg(), -2);

        f.check_all(c1.inner(), c2.inner());
    }

    #[test]
    fn test_red_cc1_neg_to_pos() { 
        let i = 0;
        let l = Link::load("5_1").unwrap().mirror().crossing_change(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, true, i);
        let f = KhComplex::cc_map1(&c1, &c2, &l, i);

        assert_eq!(f.deg(), 0);

        f.check_all(c1.inner(), c2.inner());
    }

    #[test]
    fn test_cc1_pos_to_neg_ht() { 
        type P = Poly2<'h', 't', i64>;

        let i = 0;
        let l = Link::load("5_1").unwrap().crossing_change(i);
        let (h, t) = (P::variable(0), P::variable(1));

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, false, i);
        let f = KhComplex::cc_map1(&c1, &c2, &l, i);

        assert_eq!(f.deg(), -2);

        f.check_all(c1.inner(), c2.inner());
    }

    #[test]
    fn test_cc1_neg_to_pos_ht() { 
        type P = Poly2<'h', 't', i64>;

        let i = 0;
        let l = Link::load("5_1").unwrap().mirror().crossing_change(i);
        let (h, t) = (P::variable(0), P::variable(1));

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, false, i);
        let f = KhComplex::cc_map1(&c1, &c2, &l, i);

        assert_eq!(f.deg(), 0);

        f.check_all(c1.inner(), c2.inner());
    }

    #[allow(unused)]
    fn print_h_map<R>(c1: &KhComplex<R>, c2: &KhComplex<R>, f: &KhChainMap<R>) 
    where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
        let h1 = c1.homology();
        let h2 = c2.homology();

        println!("L1");
        h1.print_seq("i");

        println!("L2");
        h2.print_seq("i");

        println!("f: deg {}\n", f.deg());

        for i in h1.h_range() { 
            let j = i + f.deg();
            println!("({i}) {} -> ({j}) {}", h1[i], h2[j]);
            for z in h1[i].gens() { 
                let w = f.apply(i, &z);
                let x = h1[i].vectorize_euc(&z).into_vec();
                let y = h2[j].vectorize_euc(&w).into_vec();
                println!("\t{:?} -> {:?}", x, y);
            }
            println!();
        }
    }
}