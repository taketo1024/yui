use num_traits::Pow;
use yui_core::PowMod2;
use yui_polynomial::LPoly;
use super::links::{Link, State};

pub fn jones_polynomial(l: &Link) -> LPoly<'q', i32> {
    type P = LPoly<'q', i32>;

    let n = l.crossing_num() as usize;
    let n_signed = l.signed_crossing_nums();
    let (n_pos, n_neg) = (n_signed.0 as i32, n_signed.1 as i32);

    let e = P::from( (-1).pow_mod2(n_neg) );
    let q = P::variable();
    let a = e * q.pow(n_pos - 2 * n_neg); // a = (-1)^{n^-} q^{n^+ - 2n^-}

    let q0: P = &q + q.pow(-1);
    let body: P = State::generate(n).into_iter().map(|s| { 
        let w = s.weight();
        let l_s = l.resolved_by(&s);
        let r = l_s.components().len();

        (-&q).pow(w) * q0.pow(r) // (-q)^w (q + q^{-1})^r
    }).sum();

    a * body
}

#[cfg(test)]
mod tests { 
    use super::*;
    use num_traits::One;
    use crate::links::Resolution;

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
        assert_eq!(p, P::from_deg([(-1, 1), (1, 1)]));
    }

    #[test]
    fn unlink_2() {
        let l = Link::from(&[[0, 1, 1, 0]]).resolved_at(0, Resolution::Res1);
        let p = jones_polynomial(&l);
        assert_eq!(p, P::from_deg([(-2, 1), (0, 2), (2, 1)]));
    }

    #[test]
    fn trefoil() {
        let l = Link::trefoil();
        let p = jones_polynomial(&l);
        assert_eq!(p, P::from_deg([(-9, -1), (-5, 1), (-3, 1), (-1, 1)]));
    }

    #[test]
    fn figure8() {
        let l = Link::figure8();
        let p = jones_polynomial(&l);
        assert_eq!(p, P::from_deg([(-5, 1), (5, 1)]));
    }

    #[test]
    fn hopf_link() { 
        let l = Link::hopf_link();
        let p = jones_polynomial(&l);
        assert_eq!(p, P::from_deg([(-6, 1), (-4, 1), (-2, 1), (0, 1)]));
    }
}