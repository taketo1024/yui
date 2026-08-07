//! The canonical (Lee/Bar-Natan) cycles of a link, built on [`KhComplex`]. Their
//! `H`-divisibility is what the `s` and `ss` invariants read off.

use yui_core::abst::{Ring, RingOps};
use yui_link::{Link, State};

use crate::ext::LinkExt;
use crate::kh::{KhAlgGen, KhChain, KhComplex, KhGen, KhTensor};

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub(crate) fn make_canon_cycles(l: &Link, a: &R, b: &R, reduced: bool) -> Vec<KhChain<R>> {
        if reduced {
            vec![
                Self::make_canon_cycle(l, a, b),
            ]
        } else {
            vec![
                Self::make_canon_cycle(l, a, b),
                Self::make_canon_cycle(l, b, a),
            ]
        }
    }

    fn make_canon_cycle(l: &Link, a: &R, b: &R) -> KhChain<R> {
        let s = l.seifert_state();
        let x_a = Self::color_factor(a); // X - a
        let x_b = Self::color_factor(b); // X - b

        let colors = l.colored_seifert_circles();
        let xs = colors.into_iter().map(|(_, c)| if c.is_a() {
            &x_a
        } else {
            &x_b
        });

        let init = KhChain::from(
            KhGen::new(s, KhTensor::empty())
        );

        xs.fold(init, |res, next| {
            res.apply_bilin(next, |a, b|
                KhGen::new(s, *a.tensor() + *b.tensor())
            )
        })
    }

    fn color_factor(a: &R) -> KhChain<R> // a -> X - a
    where R: Ring, for<'x> &'x R: RingOps<R> {
        use KhAlgGen::{I, X};

        fn init(x: KhAlgGen) -> KhGen {
            KhGen::new(
                State::empty(),
                KhTensor::from(x)
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
        let l = Link::test_data("3_1").mirror();
        let r = false;
        let c = KhComplex::new_no_simplify(&l, &1, &0, r);
        let zs = KhComplex::make_canon_cycles(&l, &0, &1, r);

        assert_eq!(zs.len(), 2);
        assert_ne!(zs[0], zs[1]);

        for z in zs {
            assert!(z.homogeneous_value(|x| c.h_deg_of(x)) == Some(0));
            assert!(!z.is_zero());

            let dz = c.d(0, &z);
            assert!(dz.is_zero());
        }
    }

    #[test]
    fn figure8() {
        let l = Link::test_data("4_1");
        let r = false;
        let c = KhComplex::new_no_simplify(&l, &1, &0, r);
        let zs = KhComplex::make_canon_cycles(&l, &0, &1, r);

        assert_eq!(zs.len(), 2);
        assert_ne!(zs[0], zs[1]);

        for z in zs {
            assert!(z.homogeneous_value(|x| c.h_deg_of(x)) == Some(0));
            assert!(!z.is_zero());

            let dz = c.d(0, &z);
            assert!(dz.is_zero());
        }
    }

    #[test]
    fn trefoil_red() {
        let l = Link::test_data("3_1").mirror();
        let r = true;
        let c = KhComplex::new_no_simplify(&l, &1, &0, r);
        let zs = KhComplex::make_canon_cycles(&l, &0, &1, r);

        assert_eq!(zs.len(), 1);

        for z in zs {
            assert!(z.homogeneous_value(|x| c.h_deg_of(x)) == Some(0));
            assert!(!z.is_zero());

            let dz = c.d(0, &z);
            assert!(dz.is_zero());
        }
    }

    #[test]
    fn figure8_red() {
        let l = Link::test_data("4_1");
        let r = true;
        let c = KhComplex::new_no_simplify(&l, &1, &0, r);
        let zs = KhComplex::make_canon_cycles(&l, &0, &1, r);

        assert_eq!(zs.len(), 1);

        for z in zs {
            assert!(z.homogeneous_value(|x| c.h_deg_of(x)) == Some(0));
            assert!(!z.is_zero());

            let dz = c.d(0, &z);
            assert!(dz.is_zero());
        }
    }
}