use yui_core::{Ring, RingOps};
use yui_link::{Edge, Link, State};

use crate::ext::LinkExt;
use crate::kh::{KhGen, KhChain, KhComplex, KhChainGen, KhTensor};

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub(crate) fn make_canon_cycles(l: &Link, base: Edge, a: &R, b: &R, reduced: bool, deg_shift: (isize, isize)) -> Vec<KhChain<R>> {
        if reduced { 
            vec![
                Self::make_canon_cycle(l, base, a, b, deg_shift),
            ]
        } else { 
            vec![
                Self::make_canon_cycle(l, base, a, b, deg_shift),
                Self::make_canon_cycle(l, base, b, a, deg_shift),
            ]
        }
    }

    fn make_canon_cycle(l: &Link, base: Edge, a: &R, b: &R, deg_shift: (isize, isize)) -> KhChain<R> {
        let s = l.seifert_state();
        let x_a = Self::color_factor(a); // X - a
        let x_b = Self::color_factor(b); // X - b

        let colors = l.colored_seifert_circles(base);
        let xs = colors.into_iter().map(|(_, c)| if c.is_a() {
            &x_a
        } else { 
            &x_b
        });

        let init = KhChain::from(
            KhChainGen::new(s, KhTensor::empty(), deg_shift)
        );

        xs.fold(init, |res, next| { 
            res.combine(next, |a, b| 
                KhChainGen::new(s, a.tensor + b.tensor, deg_shift)
            )
        })
    }

    fn color_factor(a: &R) -> KhChain<R> // a -> X - a
    where R: Ring, for<'x> &'x R: RingOps<R> { 
        use KhGen::{I, X};
    
        fn init(x: KhGen) -> KhChainGen { 
            KhChainGen::new(
                State::empty(),
                KhTensor::from(x),
                (0, 0)
            )
        }
    
        KhChain::from_iter([
            (init(X), R::one()), 
            (init(I), -a)
        ])
    }    
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use crate::kh::KhComplex;

    use super::*;
 
    #[test]
    fn trefoil() { 
        let l = Link::trefoil().mirror();
        let p = l.min_edge().unwrap();
        let r = false;
        let c = KhComplex::new_no_simplify(&l, &1, &0, r);
        let zs = KhComplex::make_canon_cycles(&l, p, &0, &1, r, c.deg_shift());

        assert_eq!(zs.len(), 2);
        assert_ne!(zs[0], zs[1]);

        for z in zs { 
            assert!(z.gens().all(|x| x.h_deg() == 0));
            assert!(!z.is_zero());
            
            let dz = c.d(0, &z);
            assert!(dz.is_zero());
        }
    }
 
    #[test]
    fn figure8() { 
        let l = Link::figure8();
        let p = l.min_edge().unwrap();
        let r = false;
        let c = KhComplex::new_no_simplify(&l, &1, &0, r);
        let zs = KhComplex::make_canon_cycles(&l, p, &0, &1, r, c.deg_shift());
        
        assert_eq!(zs.len(), 2);
        assert_ne!(zs[0], zs[1]);

        for z in zs { 
            assert!(z.gens().all(|x| x.h_deg() == 0));
            assert!(!z.is_zero());

            let dz = c.d(0, &z);
            assert!(dz.is_zero());
        }
    }

    #[test]
    fn trefoil_red() { 
        let l = Link::trefoil().mirror();
        let p = l.min_edge().unwrap();
        let r = true;
        let c = KhComplex::new_no_simplify(&l, &1, &0, r);
        let zs = KhComplex::make_canon_cycles(&l, p, &0, &1, r, c.deg_shift());

        assert_eq!(zs.len(), 1);

        for z in zs { 
            assert!(z.gens().all(|x| x.h_deg() == 0));
            assert!(!z.is_zero());
            
            let dz = c.d(0, &z);
            assert!(dz.is_zero());
        }
    }
 
    #[test]
    fn figure8_red() { 
        let l = Link::figure8();
        let p = l.min_edge().unwrap();
        let r = true;
        let c = KhComplex::new_no_simplify(&l, &1, &0, r);
        let zs = KhComplex::make_canon_cycles(&l, p, &0, &1, r, c.deg_shift());
        
        assert_eq!(zs.len(), 1);

        for z in zs { 
            assert!(z.gens().all(|x| x.h_deg() == 0));
            assert!(!z.is_zero());

            let dz = c.d(0, &z);
            assert!(dz.is_zero());
        }
    }
}