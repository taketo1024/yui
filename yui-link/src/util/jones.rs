use num_traits::Pow;
use yui::{PowMod2, GetSign, Ring};
use yui::poly::LPoly;
use crate::{Link, State};

pub fn jones_polynomial(l: &Link) -> LPoly<'q', i32> {
    type P = LPoly<'q', i32>;

    let n = l.crossing_num();
    let n_signed = l.signed_crossing_nums();
    let (n_pos, n_neg) = (n_signed.0 as i32, n_signed.1 as i32);

    let e = P::from_sign( (-1).pow_mod2(n_neg).sign() );
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
    use yui::bitseq::Bit;

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
        let l = Link::from_pd_code([[0, 1, 1, 0]]).resolved_at(0, Bit::Bit1);
        let p = jones_polynomial(&l);
        let q = P::mono;
        assert_eq!(p, P::from_iter([(q(-2), 1), (q(0), 2), (q(2), 1)]));
    }

    #[test]
    fn trefoil() {
        let l = Link::trefoil();
        let p = jones_polynomial(&l);
        let q = P::mono;
        assert_eq!(p, P::from_iter([(q(-9), -1), (q(-5), 1), (q(-3), 1), (q(-1), 1)]));
    }

    #[test]
    fn figure8() {
        let l = Link::figure8();
        let p = jones_polynomial(&l);
        let q = P::mono;
        assert_eq!(p, P::from_iter([(q(-5), 1), (q(5), 1)]));
    }

    #[test]
    fn hopf_link() { 
        let l = Link::hopf_link();
        let p = jones_polynomial(&l);
        let q = P::mono;
        assert_eq!(p, P::from_iter([(q(-6), 1), (q(-4), 1), (q(-2), 1), (q(0), 1)]));
    }
}