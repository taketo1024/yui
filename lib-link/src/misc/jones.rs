use num_traits::Pow;
use num_integer::Integer;
use yui_core::{CloneAnd, Sign, Ring, AddMon};
use yui_core::poly::{LPoly, Mono};
use yui_core::num::GaussInt;
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

// Determinant det(L) = |Δ_L(-1)| = |V_L(-1)|. `jones_polynomial` is Khovanov-unnormalized
// (J(unknot)=q+q⁻¹), so reduced V = J/(q+q⁻¹) and det = |V(q=-i)| (since t=-1 ⟹ q²=-1). At q=-i both
// J and q+q⁻¹ vanish, so by L'Hôpital det = |J'(-i)|/2 = |Σₖ k·cₖ·(-i)ᵏ| / 2 (writhe/shift-invariant).
pub fn det(l: &Link) -> i32 {
    type Z = GaussInt<i64>;
    let unit = |d: isize| match d.rem_euclid(4) {   // (-i)ᵏ as a Gaussian integer
        0 => Z::new(1, 0),
        1 => Z::new(0, -1),
        2 => Z::new(-1, 0),
        _ => Z::new(0, 1),
    };
    // det is orientation-independent (the formula is shift-invariant), so orient freely if needed
    // (`reorient` with every port claimable always succeeds on a valid diagram).
    let jones = if l.is_oriented() {
        jones_polynomial(l)
    } else {
        jones_polynomial(&l.clone_and(|l| {
            l.reorient(|_, _| true);
        }))
    };

    let s = jones.iter().fold(Z::from(0), |s, (x, c)| {
        let d = x.deg();
        s + Z::from(d as i64 * *c as i64) * unit(d)
    });

    // s is purely real or imaginary (a real value × the unit (-i)^writhe); which one depends on the
    // writhe, not just #components, since this Jones is unnormalized.
    let v = match s.pair_into() {
        (0, m) | (m, 0) => m.abs() as i32,
        _ => unreachable!("Jones at q=-i must be purely real or imaginary"),
    };
    debug_assert!(v.is_even(), "determinant numerator must be even");
    v / 2
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
        let l = Link::test_data("3_1").mirror();
        let p = jones_polynomial(&l);
        let q = P::mono;
        assert_eq!(p, P::from_iter([(q(-9), -1), (q(-5), 1), (q(-3), 1), (q(-1), 1)]));
    }

    #[test]
    fn determinant_values() {
        assert_eq!(det(&Link::unknot()), 1);
        assert_eq!(det(&Link::test_data("3_1")), 3);
        assert_eq!(det(&Link::test_data("4_1")), 5);
        assert_eq!(det(&Link::test_data("5_2")), 7);
        assert_eq!(det(&Link::test_data("7_1")), 7);   // (2,7) torus
        assert_eq!(det(&Link::test_data("L4a1")), 4);  // (2,4) torus link: |H₁(L(4,1))| = 4
        assert_eq!(det(&Link::test_data("unlink2")), 0);  // split (and unoriented: det orients freely)
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