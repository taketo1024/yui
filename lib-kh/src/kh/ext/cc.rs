use itertools::Itertools;
use yui_core::bitseq::Bit;
use yui_core::lc::Lc;
use yui_core::{Ring, RingOps, Sign};
use yui_homology::ChainMap;
use yui_link::{Link, Path, State};
use num_traits::Zero;

use crate::kh::{KhAlg, KhCube, KhChain, KhGen, KhComplex, KhAlgGen, KhTensor};

pub type KhChainMap<'a, 'c, R> = ChainMap<'a, 'c, isize, KhGen, KhGen, R>;

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn cc_pair(l: &Link, h: &R, t: &R, reduced: bool, i: usize) -> (KhComplex<R>, KhComplex<R>) {
        assert!(l.node(i).is_crossing());
        
        let l2 = l.cc_at(i);
        let c1 = KhComplex::new_no_simplify(l, h, t, reduced);
        let c2 = KhComplex::new_no_simplify(&l2, h, t, reduced);
        (c1, c2)
    }

    pub fn cc_map0<'a>(c1: &'a KhComplex<R>, c2: &'a KhComplex<R>, i: usize) -> KhChainMap<'a, 'static, R> {
        let deg = c2.deg_shift().0 - c1.deg_shift().0 + 1;

        ChainMap::new(c1.inner(), c2.inner(), deg, move |_, z| {
            z.apply(|x: &KhGen| {
                if !x.state()[i].is_zero() {
                    return KhChain::zero();
                }

                let e = Sign::from_parity( count_1s(x.state(), i) );
                let t = x.state().edit(|s| s.set_1(i));
                let y = KhGen::new(t, *x.tensor());

                KhChain::from(y) * R::from_sign(e)
            })
        })
    }

    pub fn cc_map1<'a>(c1: &'a KhComplex<R>, c2: &'a KhComplex<R>, l: &Link, i: usize) -> KhChainMap<'a, 'static, R> {
        assert!(l.node(i).is_crossing());
        assert!(!c1.is_reduced() || l.base_pt().is_some());

        let deg = c2.deg_shift().0 - c1.deg_shift().0 - 1;

        let alg = c1.alg().clone();
        let (a0, a1) = l.node(i).resolve(Bit::Bit0).arcs();

        // TODO We don't want to reproduce the cube.
        let (h, t) = c1.alg().ht();
        let base_pt = if c1.is_reduced() { l.base_pt() } else { None };
        let cube = KhCube::new(l, h, t, base_pt, c1.deg_shift());

        ChainMap::new(c1.inner(), c2.inner(), deg, move |_, z| {
            z.apply(|x: &KhGen| {
                if !x.state()[i].is_one() {
                    return KhChain::zero();
                }

                let circles = cube.vertex(x.state()).circles();
                let (k0, k1) = (circle_index(circles, &a0), circle_index(circles, &a1));

                if k0 == k1 {
                    return Lc::zero();
                }

                let s = x.state().edit(|s| s.set_0(i));

                let e = Sign::from_parity( count_1s(x.state(), i) );
                let t = apply_f1(&alg, x.tensor(), k0, k1) * R::from_sign(e);

                t.map_keys(|y| {
                    KhGen::new(s, y)
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
    use KhAlgGen::X;

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
    use yui_core::poly::Poly2;
    use yui_core::{EucRing, EucRingOps};
    use yui_homology::ToSeqString;
    use yui_link::Link;

    use crate::kh::ext::cc::KhChainMap;
    use crate::kh::KhComplex;
 
    #[test]
    fn test_cc0_pos_to_neg() { 
        let i = 0;
        let l = Link::test_data("5_1").cc_at(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, false, i);
        let f = KhComplex::cc_map0(&c1, &c2, i);

        assert_eq!(f.deg(), 0);

        f.check_all();
    }

    #[test]
    fn test_cc0_neg_to_pos() { 
        let i = 0;
        let l = Link::test_data("5_1").mirror().cc_at(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, false, i);
        let f = KhComplex::cc_map0(&c1, &c2, i);

        assert_eq!(f.deg(), 2);

        f.check_all();
    }

    #[test]
    fn test_cc1_pos_to_neg() { 
        let i = 0;
        let l = Link::test_data("5_1").cc_at(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, false, i);
        let f = KhComplex::cc_map1(&c1, &c2, &l, i);

        assert_eq!(f.deg(), -2);

        f.check_all();
    }

    #[test]
    fn test_cc1_neg_to_pos() { 
        let i = 0;
        let l = Link::test_data("5_1").mirror().cc_at(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, false, i);
        let f = KhComplex::cc_map1(&c1, &c2, &l, i);

        assert_eq!(f.deg(), 0);

        f.check_all();
    }

    #[test]
    fn test_red_cc0_pos_to_neg() { 
        let i = 0;
        let l = Link::test_data("5_1").cc_at(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, true, i);
        let f = KhComplex::cc_map0(&c1, &c2, i);

        assert_eq!(f.deg(), 0);

        f.check_all();
    }

    #[test]
    fn test_red_cc0_neg_to_pos() { 
        let i = 0;
        let l = Link::test_data("5_1").mirror().cc_at(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, true, i);
        let f = KhComplex::cc_map0(&c1, &c2, i);

        assert_eq!(f.deg(), 2);

        f.check_all();
    }

    #[test]
    fn test_red_cc1_pos_to_neg() { 
        let i = 0;
        let l = Link::test_data("5_1").cc_at(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, true, i);
        let f = KhComplex::cc_map1(&c1, &c2, &l, i);

        assert_eq!(f.deg(), -2);

        f.check_all();
    }

    #[test]
    fn test_red_cc1_neg_to_pos() { 
        let i = 0;
        let l = Link::test_data("5_1").mirror().cc_at(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, true, i);
        let f = KhComplex::cc_map1(&c1, &c2, &l, i);

        assert_eq!(f.deg(), 0);

        f.check_all();
    }

    #[test]
    fn test_cc1_pos_to_neg_ht() { 
        type P = Poly2<'h', 't', i64>;

        let i = 0;
        let l = Link::test_data("5_1").cc_at(i);
        let (h, t) = (P::variable(0), P::variable(1));

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, false, i);
        let f = KhComplex::cc_map1(&c1, &c2, &l, i);

        assert_eq!(f.deg(), -2);

        f.check_all();
    }

    #[test]
    fn test_cc1_neg_to_pos_ht() { 
        type P = Poly2<'h', 't', i64>;

        let i = 0;
        let l = Link::test_data("5_1").mirror().cc_at(i);
        let (h, t) = (P::variable(0), P::variable(1));

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, false, i);
        let f = KhComplex::cc_map1(&c1, &c2, &l, i);

        assert_eq!(f.deg(), 0);

        f.check_all();
    }

    #[allow(unused)]
    fn print_h_map<R>(c1: &KhComplex<R>, c2: &KhComplex<R>, f: &KhChainMap<'_, '_, R>)
    where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
        let h1 = c1.homology();
        let h2 = c2.homology();

        println!("L1");
        println!("{}", h1.to_seq_string());

        println!("L2");
        println!("{}", h2.to_seq_string());

        println!("f: deg {}\n", f.deg());

        for i in h1.h_range() { 
            let j = i + f.deg();
            println!("({i}) {} -> ({j}) {}", h1[i], h2[j]);
            for z in h1[i].generators() {
                let w = f.apply(i, &z);
                let x = h1[i].vectorize_euc(&z).into_dense();
                let y = h2[j].vectorize_euc(&w).into_dense();
                println!("\t{:?} -> {:?}", x, y);
            }
            println!();
        }
    }
}