use num_traits::{Pow, Zero};
use crate::math::laurent_polynomial::LaurentPolynomial;
use crate::math::traits::PowMod2;
use crate::links::links::{Link, State};

pub fn jones_polynomial(l: &Link) -> LaurentPolynomial<i32> {
    type P = LaurentPolynomial<i32>;

    let n = l.crossing_num() as usize;
    let n_signed = l.signed_crossing_nums();
    let (n_pos, n_neg) = (n_signed.0 as i32, n_signed.1 as i32);

    let q    = P::variable();
    let qinv = P::new(vec![1], -1);
    
    let e = P::constant( (-1).pow_mod2(n_neg) );
    let s = n_pos - 2 * n_neg;
    let a = match s >= 0 {
        true => q.pow(s),
        false => qinv.pow(-s)
    };

    let m = 2.pow(n) as usize;
    let body = (0..m).fold(P::zero(), |res, i| { 
        let s = State::from_bseq(i, n);
        let l_s = l.clone().resolve(&s);

        let w = s.weight() as i32;
        let r = l_s.components().len() as i32;
        let p = (-&q).pow(w) * (&q + &qinv).pow(r);

        res + p
    });

    e * a * body
}

#[cfg(test)]
mod tests { 
    use num_traits::One;

    use crate::math::laurent_polynomial::LaurentPolynomial;
    use crate::links::links::{Link, Resolution::Res1};
    use super::jones_polynomial;

    type P = LaurentPolynomial<i32>;

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
        assert_eq!(p, P::new(vec![1, 0, 1], -1));
    }

    #[test]
    fn unlink_2() {
        let l = Link::from([[0, 1, 1, 0]]).resolve_at(0, Res1);
        let p = jones_polynomial(&l);
        assert_eq!(p, P::new(vec![1, 0, 2, 0, 1], -2));
    }

    #[test]
    fn trefoil() {
        let l = Link::trefoil();
        let p = jones_polynomial(&l);
        assert_eq!(p, P::new(vec![-1, 0, 0, 0, 1, 0, 1, 0, 1], -9));
    }

    #[test]
    fn figure8() {
        let l = Link::figure8();
        let p = jones_polynomial(&l);
        assert_eq!(p, P::new(vec![1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], -5));
    }

    #[test]
    fn hopf_link() { 
        let l = Link::hopf_link();
        let p = jones_polynomial(&l);
        assert_eq!(p, P::new(vec![1, 0, 1, 0, 1, 0, 1], -6));
    }
}