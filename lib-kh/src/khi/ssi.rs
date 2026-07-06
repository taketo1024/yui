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
use yui_homology::algo::{ChainReducer, HomologyCalc};
use yui_matrix::MatTrait;
use yui_matrix::sparse::SpMat;
use yui_link::{InvLink, Link};

use crate::kh::KhComplex;
use crate::tng::builder::SymBuildConfig;
use crate::util::calc::div_vec;
use crate::khi::{KhIComplex, KhIHomology};

// injected into the reducer on the heavy cone path: bounds each Schur round's `a⁻¹b` transient.
pub(crate) const MAX_PIVOTS_PER_ROUND: usize = 32_768;

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

/// `ssi` via the cobordism-level cone (`ConeBuilder`), without any bigraded structure: the canon
/// classes are transported as coordinate vectors through a trans-free capped reduction, and the
/// basis-change is computed only at the reduced scale (two small SNFs). Memory-safe on huge diagrams.
pub fn ssi_invariants_via_cone<R>(l: &InvLink, c: &R, reduced: bool, config: SymBuildConfig) -> (i32, i32)
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    assert!(!c.is_zero());
    assert!(!c.is_unit());
    assert!(l.is_knot());

    info!("compute ssi via cone, c = {c} over {}.", R::math_symbol());

    let w = l.writhe();
    let r = l.seifert_circles().len() as i32;
    let (d0, d1) = div_via_cone(l, c, reduced, config);

    let ss0 = 2 * d0 + w - r + 1;
    let ss1 = 2 * d1 + w - r + 1;

    info!("w = {w}, r = {r}, d0 = {d0}, d1 = {d1}.");
    info!("ssi = ({ss0}, {ss1}).");

    (ss0, ss1)
}

fn div_via_cone<R>(l: &InvLink, c: &R, reduced: bool, config: SymBuildConfig) -> (i32, i32)
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    let r = if reduced { 1 } else { 2 };
    let t = R::zero();

    // Respect a requested `h_range` (must cover 0 and 1, the ssi degrees); else the full support up to 1.
    // Built one degree wider on both ends for the boundary maps.
    let requested = config.h_range.clone().unwrap_or(-(Link::MAX_CROSSING as isize) ..= 1);
    let range = KhComplex::<R>::clamp_h_range(l.inner(), reduced, requested);
    let (a, b) = (*range.start(), *range.end());
    assert!(a <= 0 && b >= 1, "ssi h-range must include 0 and 1, got {a}..={b}");
    let config = SymBuildConfig { h_range: Some((a - 1)..=(b + 1)), ..config };
    let kc = KhIComplex::from_cone(l, c, &t, reduced, config);

    let zs = kc.canon_cycles(); // sorted by h-degree: B classes at 0, then Q classes at 1
    assert_eq!(zs.len(), 2 * r);
    for (i, z) in zs.iter().enumerate() {
        let expected = if i < r { 0 } else { 1 };
        assert!(!z.is_zero());
        assert_eq!(z.homogeneous_value(|x| kc.h_deg_of(x)), Some(expected));
    }

    let mut red = ChainReducer::from_complex(kc.inner(), false);
    red.set_max_pivots(MAX_PIVOTS_PER_ROUND);

    for z in zs.iter() {
        let h = kc.h_deg_of_chain(z);
        red.add_vec(h, kc.inner()[h].vectorize(z));
    }

    red.reduce_all(false);
    red.reduce_all(true);

    // basis change only at the reduced scale: homology with trans at the two canon degrees.
    let ds = [0, 1].map(|h| {
        let d1 = red.matrix(h).expect("d[h] must be set").clone();
        // below the clamped window bottom there is no incoming differential.
        let d0 = red.matrix(h - 1).cloned().unwrap_or_else(|| SpMat::zero((d1.n_cols(), 0)));
        let (rank, _, tr) = HomologyCalc::calculate(d0, d1, true);
        let tr = tr.unwrap();

        assert_eq!(rank, r);
        info!("KhI[{h}]: rank {rank}");

        red.vecs(h).expect("transported vecs at canon degree").iter().map(|v| {
            let w = tr.forward(v);
            info!("a in KhI[{h}]: ({})", w.clone().into_dense().iter().join(","));
            div_vec(&w.subvec(0..r), c).expect("invalid divisibility.")
        }).collect_vec()
    });

    let (d0, d1) = if reduced {
        (ds[0][0], ds[1][0])
    } else {
        assert_eq!(ds[0][0], ds[0][1]);
        assert_eq!(ds[1][0], ds[1][1]);
        (ds[0][0], ds[1][0])
    };

    assert!(d0 <= d1);

    (d0, d1)
}

fn div<R>(l: &InvLink, c: &R, reduced: bool) -> (i32, i32)
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    let r = if reduced { 1 } else { 2 };
    let t = R::zero(); 

    // bottom..=1: building the cheap low degrees and truncating only at the top is faster than
    // the doubly-truncated `0..=1` slice (which widens to the dense `-1..=2`). Builder clamps the start.
    let kh = KhIHomology::new_partial(l, c, &t, reduced, Some(-(Link::MAX_CROSSING as isize) ..= 1));

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

    macro_rules! test_cone {
        ($(#[$m:meta])* $test:ident, $name:literal, $expected:expr) => {
            $(#[$m])*
            #[test]
            fn $test() -> Result<(), Box<dyn std::error::Error>> {
                let c = P::variable();
                let l = InvLink::load($name)?;
                let ssi = ssi_invariants_via_cone(&l, &c, false, SymBuildConfig::default());
                assert_eq!(ssi, $expected);

                Ok(())
            }
        }
    }

    // all-negative diagram: the Kh h-range tops at 0, so the KhI top degree 1 must survive the
    // cone's window pruning (regression: prune_isolated_top dropped it, killing the Q classes).
    #[test]
    fn k3_1_m_cone() {
        let l = InvLink::test_data("3_1").mirror();
        let c = P::variable();

        let ssi = ssi_invariants_via_cone(&l, &c, false, SymBuildConfig::default());
        assert_eq!(ssi, (-2, -2));
    }

    // direct symmetry-broken emission (Sano2026 Prop 4.6) with the SDR canon transport.
    macro_rules! test_cone_direct {
        ($test:ident, $name:literal, $expected:expr) => {
            #[test]
            fn $test() -> Result<(), Box<dyn std::error::Error>> {
                let c = P::variable();
                let l = InvLink::load($name)?;
                let config = SymBuildConfig { ..Default::default() };
                let ssi = ssi_invariants_via_cone(&l, &c, false, config);
                assert_eq!(ssi, $expected);

                Ok(())
            }
        }
    }

    test_cone_direct!(k3_1_cone_direct, "3_1", (2, 2));
    test_cone_direct!(k4_1_cone_direct, "4_1", (0, 0));
    test_cone_direct!(k6_2a_cone_direct, "6_2a", (2, 2));
    test_cone_direct!(k7_6a_cone_direct, "7_6a", (-2, -2));

    #[test]
    fn k3_1_m_cone_direct() {
        let l = InvLink::test_data("3_1").mirror();
        let c = P::variable();

        let config = SymBuildConfig { ..Default::default() };
        let ssi = ssi_invariants_via_cone(&l, &c, false, config);
        assert_eq!(ssi, (-2, -2));
    }

    #[test]
    fn k9_46_cone_direct() {
        let l = InvLink::from_symmetric_pd_code(
            [[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]
        );

        let c = P::variable();
        let config = SymBuildConfig { ..Default::default() };
        let ssi = ssi_invariants_via_cone(&l, &c, false, config);

        assert_eq!(ssi, (0, 2));
    }

    test_cone!(k3_1_cone, "3_1", (2, 2));
    test_cone!(k4_1_cone, "4_1", (0, 0));
    test_cone!(k5_1_cone, "5_1", (4, 4));
    test_cone!(k6_2a_cone, "6_2a", (2, 2));
    test_cone!(k6_3_cone, "6_3", (0, 0));
    test_cone!(k7_6a_cone, "7_6a", (-2, -2));
    test_cone!(k7_7a_cone, "7_7a", (0, 0));

    #[test]
    fn k3_1_cone_red() {
        let l = InvLink::test_data("3_1");
        let c = P::variable();

        let ssi = ssi_invariants_via_cone(&l, &c, true, SymBuildConfig::default());
        assert_eq!(ssi, (2, 2));
    }

    #[test]
    fn k9_46_cone() {
        let l = InvLink::from_symmetric_pd_code(
            [[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]
        );

        let c = P::variable();
        let ssi = ssi_invariants_via_cone(&l, &c, false, SymBuildConfig::default());

        assert_eq!(ssi, (0, 2));
    }

    #[test]
    fn k15n_103488_cone() {
        let l = InvLink::from_symmetric_pd_code(
            [[1,11,2,10],[2,20,3,19],[5,17,6,16],[6,25,7,26],[9,22,10,23],[12,30,13,29],[14,8,15,7],[15,27,16,26],[18,4,19,3],[20,11,21,12],[21,1,22,30],[23,4,24,5],[24,18,25,17],[27,8,28,9],[28,14,29,13]]
        );

        let c = P::variable();
        let ssi = ssi_invariants_via_cone(&l, &c, false, SymBuildConfig::default());

        assert_eq!(ssi, (0, 2));
    }

    #[test]
    fn k15n_103488() { 
        let l = InvLink::from_symmetric_pd_code(
            [[1,11,2,10],[2,20,3,19],[5,17,6,16],[6,25,7,26],[9,22,10,23],[12,30,13,29],[14,8,15,7],[15,27,16,26],[18,4,19,3],[20,11,21,12],[21,1,22,30],[23,4,24,5],[24,18,25,17],[27,8,28,9],[28,14,29,13]]
        );

        let c = P::variable();
        let ssi = ssi_invariants(&l, &c, false);

        assert_eq!(ssi, (0, 2));
    }

    #[test]
    #[ignore]
    fn k17nh_73() {
        let l = InvLink::from_symmetric_pd_code(
            [[1,27,2,26],[19,2,20,3],[3,13,4,12],[4,31,5,32],[30,5,31,6],[13,7,14,6],[8,27,9,28],[9,1,10,34],[10,18,11,17],[24,11,25,12],[14,21,15,22],[28,16,29,15],[33,16,34,17],[18,26,19,25],[20,8,21,7],[29,23,30,22],[23,33,24,32]]
        );

        let c = P::variable();
        let ssi = ssi_invariants(&l, &c, false);

        assert_eq!(ssi, (0, 2));
    }
}