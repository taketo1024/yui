use yui::bitseq::Bit;
use yui::{Ring, RingOps, Sign};
use yui_homology::ChainMap;
use yui_link::Link;
use num_traits::Zero;

use crate::kh::{KhChain, KhChainGen, KhComplex};

type KhChainMap<R> = ChainMap<isize, KhChainGen, KhChainGen, R>;

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn cc_pair(l: &Link, i: usize, h: &R, t: &R, reduced: bool) -> (KhComplex<R>, KhComplex<R>) {
        let l2 = l.crossing_changed_at(i);
        let c1 = KhComplex::new_no_simplify(l, h, t, reduced);
        let c2 = KhComplex::new_no_simplify(&l2, h, t, reduced);
        (c1, c2)
    }

    pub fn cc_map0(c1: &KhComplex<R>, c2: &KhComplex<R>, i: usize) -> KhChainMap<R> {
        let deg = c2.deg_shift().0 - c1.deg_shift().0 + 1;
        let c2_deg_shift = c2.deg_shift();
        ChainMap::new(c1.inner(), c2.inner(), deg, move |_, z| { 
            z.apply(|x| {
                if x.state[i].is_one() { 
                    return KhChain::zero();
                }

                let e = Sign::from_parity(
                    x.state.iter().enumerate().filter(|(j, b)| j > &i && b.is_one()).count() as i32
                );
                let mut y = x.clone();
                y.state.set_1(i);
                y.deg_shift = c2_deg_shift;

                KhChain::from(y) * R::from_sign(e)
            })
        })
    }
}

#[cfg(test)]
mod tests {
    use yui::{EucRing, EucRingOps};
    use yui_homology::DisplaySeq;
    use yui_link::Link;

    use crate::kh::ext::cc::KhChainMap;
    use crate::kh::KhComplex;
 
    #[test]
    fn test_cc0_pos_to_neg() { 
        let i = 0;
        let l = Link::load("5_1").unwrap().crossing_changed_at(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, i, &h, &t, false);
        let f = KhComplex::cc_map0(&c1, &c2, i);

        assert_eq!(f.deg(), 0);

        f.check_all(c1.inner(), c2.inner());
    }

    #[test]
    fn test_cc0_neg_to_pos() { 
        let i = 0;
        let l = Link::load("5_1").unwrap().mirror().crossing_changed_at(i);
        let (h, t) = (0, 0);

        let (c1, c2) = KhComplex::cc_pair(&l, i, &h, &t, false);
        let f = KhComplex::cc_map0(&c1, &c2, i);

        assert_eq!(f.deg(), 2);

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