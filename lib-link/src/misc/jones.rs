use num_traits::Pow;
use yui_core::{Sign, Ring, AddMon};
use yui_core::poly::LPoly;
use crate::{Link, State};

pub fn jones_polynomial(l: &Link) -> LPoly<'q', i32> {
    type P = LPoly<'q', i32>;

    // the writhe normalization needs the orientation; unoriented signs would silently read as 0.
    assert!(l.is_oriented(), "jones_polynomial requires an oriented link");

    let n = l.n_crossings();
    let n_signed = l.n_signed_crossings();
    let (n_pos, n_neg) = (n_signed.0 as i32, n_signed.1 as i32);

    let e = P::from_sign( Sign::from_parity(n_neg) );
    let q = P::variable();
    let a = e * q.pow(n_pos - 2 * n_neg); // a = (-1)^{n^-} q^{n^+ - 2n^-}

    let q0: P = &q + q.pow(-1);
    let body = P::sum(State::generate(n).map(|s| { 
        let w = s.weight();
        let l_s = l.resolve_by(&s);
        let r = l_s.n_comps();

        (-&q).pow(w) * q0.pow(r) // (-q)^w (q + q^{-1})^r
    }));

    a * body
}

#[cfg(test)]
mod tests { 
    use super::*;
    use num_traits::One;
    use yui_core::bitseq::Bit;

    type P = LPoly<'q', i32>;

    #[test]
    fn empty() {
        let l = Link::empty();
        let p = jones_polynomial(&l);
        assert_eq!(p, P::one());
    }

    #[test]
    fn unknot() {
        let l = Link::unknot();
        let p = jones_polynomial(&l);
        let q = P::mono;
        assert_eq!(p, P::from_iter([(q(-1), 1), (q(1), 1)]));
    }

    #[test]
    fn unlink_2() {
        let l = Link::test_data("unknot_r_twist").resolve_at(0, Bit::Bit1);
        let p = jones_polynomial(&l);
        let q = P::mono;
        assert_eq!(p, P::from_iter([(q(-2), 1), (q(0), 2), (q(2), 1)]));
    }

    #[test]
    fn trefoil() {
        let l = Link::test_data("3_1");
        let p = jones_polynomial(&l);
        let q = P::mono;
        assert_eq!(p, P::from_iter([(q(-9), -1), (q(-5), 1), (q(-3), 1), (q(-1), 1)]));
    }

    #[test]
    fn figure8() {
        let l = Link::test_data("4_1");
        let p = jones_polynomial(&l);
        let q = P::mono;
        assert_eq!(p, P::from_iter([(q(-5), 1), (q(5), 1)]));
    }

    #[test]
    fn hopf_link() {
        let l = Link::test_data("L2a1");
        let p = jones_polynomial(&l);
        let q = P::mono;
        assert_eq!(p, P::from_iter([(q(-6), 1), (q(-4), 1), (q(-2), 1), (q(0), 1)]));
    }

    #[test]
    fn unlink_2_loops() {
        // J(unlink(2)) = (q + q^{-1})^2 = q^{-2} + 2 + q^2.
        let l = Link::unlink(2);
        let p = jones_polynomial(&l);
        let q = P::mono;
        assert_eq!(p, P::from_iter([(q(-2), 1), (q(0), 2), (q(2), 1)]));
    }

    #[test]
    fn unlink_3_loops() {
        // J(unlink(3)) = (q + q^{-1})^3 = q^{-3} + 3 q^{-1} + 3 q + q^3.
        let l = Link::unlink(3);
        let p = jones_polynomial(&l);
        let q = P::mono;
        assert_eq!(p, P::from_iter([(q(-3), 1), (q(-1), 3), (q(1), 3), (q(3), 1)]));
    }

    #[test]
    fn trefoil_with_free_loop() {
        // J(L ⊔ unknot) = J(L) · (q + q^{-1}).
        let trefoil = Link::test_data("3_1");
        let p_trefoil = jones_polynomial(&trefoil);

        let l = Link::new(trefoil.nodes().cloned(), [7]);
        let p = jones_polynomial(&l);

        let q = P::mono;
        let q_factor = P::from_iter([(q(-1), 1), (q(1), 1)]);
        assert_eq!(p, p_trefoil * q_factor);
    }
}