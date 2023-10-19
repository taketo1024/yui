use yui_core::{Ring, RingOps};
use yui_link::{Link, State};

use crate::ext::LinkExt;
use crate::{KhAlgGen, KhLabel, KhEnhState, KhChain};

pub trait CanonCycles<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    fn canon_cycle(l: &Link, a: &R, b: &R, ori: bool) -> KhChain<R>;
}

impl<R> CanonCycles<R> for KhChain<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn canon_cycle(l: &Link, a: &R, b: &R, ori: bool) -> KhChain<R> {
        let s = l.ori_pres_state();
        let colors = l.colored_seifert_circles(ori);

        let mut z = KhChain::from_gen(
            KhEnhState::new(s, KhLabel::empty())
        );

        let x_a = make_factor(a); // X - a
        let x_b = make_factor(b); // X - b
    
        for (_, c) in colors {
            z *= if c.is_a() {
                &x_a
            } else { 
                &x_b
            }
        }
    
        z
    }
}

fn make_factor<R>(a: &R) -> KhChain<R> // a -> X - a
where R: Ring, for<'x> &'x R: RingOps<R> { 
    use KhAlgGen::{I, X};

    fn init(x: KhAlgGen) -> KhEnhState { 
        KhEnhState::new(
            State::empty(),
            KhLabel::from(x)
        )
    }

    KhChain::from_iter([
        (init(X), R::one()), 
        (init(I), -a)
    ])
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use crate::KhComplex;

    use super::*;
 
    #[test]
    fn trefoil() { 
        let l = Link::trefoil().mirror();
        let c = KhComplex::new_v1(&l, &1, &0, false);

        let zs = [true, false].map(|ori| 
            KhChain::canon_cycle(&l, &0, &1, ori)
        );

        for z in zs { 
            let dz = c.differentiate(0, &z);
            assert_eq!(z.is_zero(), false);
            assert_eq!(dz.is_zero(), true);
        }
    }
 
    #[test]
    fn figure8() { 
        let l = Link::figure8();
        let c = KhComplex::new_v1(&l, &1, &0, false);
        
        let zs = [true, false].map(|ori| 
            KhChain::canon_cycle(&l, &0, &1, ori)
        );

        for z in zs { 
            let dz = c.differentiate(0, &z);
            assert_eq!(z.is_zero(), false);
            assert_eq!(dz.is_zero(), true);
        }
    }
}