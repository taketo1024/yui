//! The equivariant Rasmussen invariant `(s̲, s̄)` for a strongly invertible knot,
//! obtained from the `c`-divisibilities of the two equivariant Lee classes in
//! `KhI` (§3 of the reference).
//!
//! Reference:
//! - T. Sano, "Involutive Khovanov homology and equivariant knots",
//!   Algebr. Geom. Topol. 25 (2025), 5059–5111.
//!   <https://doi.org/10.2140/agt.2025.25.5059>, <https://arxiv.org/abs/2404.08568>

use itertools::Itertools;
use num_traits::Zero;
use log::info;

use yui_core::{EucRing, EucRingOps};
use yui_link::InvLink;

use crate::util::calc::div_vec;
use crate::khi::KhIHomology;

pub fn ssi_invariants<R>(l: &InvLink, c: &R, reduced: bool) -> (i32, i32)
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    assert!(!c.is_zero());
    assert!(!c.is_unit());
    assert!(l.is_knot());

    info!("compute ssi, c = {c} over {}.", R::math_symbol());

    let w = l.writhe();
    let r = l.seifert_circles().len() as i32;
    let (d0, d1) = div(l, c, reduced);

    let ss0 = 2 * d0 + w - r + 1;
    let ss1 = 2 * d1 + w - r + 1;

    info!("w = {w}, r = {r}, d0 = {d0}, d1 = {d1}.");
    info!("ssi = ({ss0}, {ss1}).");

    (ss0, ss1)
}

fn div<R>(l: &InvLink, c: &R, reduced: bool) -> (i32, i32)
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    let r = if reduced { 1 } else { 2 };
    let t = R::zero(); 

    // bottom..=1: building the cheap low degrees and truncating only at the top is faster than
    // the doubly-truncated `0..=1` slice (which widens to the dense `-1..=2`). Builder clamps the start.
    let kh = KhIHomology::new_partial(l, c, &t, reduced, Some(isize::MIN + 1 ..= 1));

    assert_eq!(kh[0].rank(), r);
    assert_eq!(kh[1].rank(), r);    

    info!("KhI[0]: {}", kh[0]);    
    info!("KhI[1]: {}", kh[1]);    

    let zs = kh.canon_cycles();
    
    assert_eq!(zs.len(), 2 * r);
    for (i, z) in zs.iter().enumerate() {
        let expected = if i < r { 0 } else { 1 };
        assert!(!z.is_zero());
        assert_eq!(z.homogeneous_value(|x| kh.h_deg_of(x)), Some(expected));
    }

    let ds = zs.iter().enumerate().map(|(i, z)| {
        let h = kh.h_deg_of_chain(z);
        let v = kh[h].vectorize_euc(z);
        info!("a[{i}] in Kh[{h}]: ({})", v.clone().into_dense().iter().join(","));
        v
    }).map(|v| 
        div_vec(&v.subvec(0..r), c).expect("invalid divisibility.")
    ).collect_vec();

    let (d0, d1) = if reduced { 
        (ds[0], ds[1])
    } else { 
        assert_eq!(ds[0], ds[1]);
        assert_eq!(ds[2], ds[3]);
        (ds[0], ds[2])
    };

    assert!(d0 <= d1);

    (d0, d1)
}

#[cfg(test)]
mod tests {
    use yui_core::poly::Poly;
    use yui_core::num::FF2;

    use super::*;

    type R = FF2;
    type P = Poly<'H', R>;

    #[test]
    fn test_unknot_pos_twist() {
        let l = InvLink::test_data("unknot_r_twist");
        let c = P::variable();

        let ssi = ssi_invariants(&l, &c, false);
        assert_eq!(ssi.0, 0);
        assert_eq!(ssi.1, 0);
    }

    #[test]
    fn test_unknot_neg_twist() {
        let l = InvLink::test_data("unknot_l_twist");
        let c = P::variable();

        let ssi = ssi_invariants(&l, &c, false);
        assert_eq!(ssi.0, 0);
        assert_eq!(ssi.1, 0);
    }

    #[test]
    fn test_unknot_neg_twist2() {
        let l = InvLink::test_data("unknot_l_twist2");
        let c = P::variable();

        let ssi = ssi_invariants(&l, &c, false);
        assert_eq!(ssi.0, 0);
        assert_eq!(ssi.1, 0);
    }

    #[test]
    fn test_3_1() { 
        let l = InvLink::test_data("3_1");
        let c = P::variable();

        let ssi = ssi_invariants(&l, &c, false);
        assert_eq!(ssi.0, 2);
        assert_eq!(ssi.1, 2);
    }

    #[test]
    fn test_3_1_m() { 
        let l = InvLink::test_data("3_1").mirror();
        let c = P::variable();

        let ssi = ssi_invariants(&l, &c, false);
        assert_eq!(ssi.0, -2);
        assert_eq!(ssi.1, -2);
    }

    #[test]
    fn test_3_1_red() { 
        let l = InvLink::test_data("3_1");
        let c = P::variable();

        let ssi = ssi_invariants(&l, &c, true);
        assert_eq!(ssi.0, 2);
        assert_eq!(ssi.1, 2);
    }

    macro_rules! test {
        ($(#[$m:meta])* $test:ident, $name:literal, $expected:expr) => {
            $(#[$m])* 
            #[test]
            fn $test() -> Result<(), Box<dyn std::error::Error>> { 
                type R = FF2;
                type P = Poly<'H', R>;
                let c = P::variable();
    
                let l = InvLink::load($name)?;
                let ssi = ssi_invariants(&l, &c, false);
                assert_eq!(ssi, $expected);
    
                Ok(())
            }
        }
    }
    
    test!(k3_1, "3_1", (2, 2));
    test!(k4_1, "4_1", (0, 0));
    test!(k5_1, "5_1", (4, 4));
    test!(k5_2a, "5_2a", (2, 2));
    test!(k5_2b, "5_2b", (2, 2));
    test!(k6_1a, "6_1a", (0, 0));
    test!(k6_1b, "6_1b", (0, 0));
    test!(k6_2a, "6_2a", (2, 2));
    test!(k6_2b, "6_2b", (2, 2));
    test!(k6_3, "6_3", (0, 0));
    test!(k7_1, "7_1", (6, 6));
    test!(k7_2a, "7_2a", (2, 2));
    test!(k7_2b, "7_2b", (2, 2));
    test!(k7_3a, "7_3a", (4, 4));
    test!(k7_3b, "7_3b", (4, 4));
    test!(k7_4a, "7_4a", (2, 2));
    test!(k7_4b, "7_4b", (2, 2));
    test!(k7_5a, "7_5a", (4, 4));
    test!(k7_5b, "7_5b", (4, 4));
    test!(k7_6a, "7_6a", (-2, -2));
    test!(k7_6b, "7_6b", (-2, -2));
    test!(k7_7a, "7_7a", (0, 0));
    test!(k7_7b, "7_7b", (0, 0));

    #[test]
    fn k9_46() { 
        let l = InvLink::from_symmetric_pd_code(
            [[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]
        );

        let c = P::variable();
        let ssi = ssi_invariants(&l, &c, false);

        assert_eq!(ssi, (0, 2));
    }
}