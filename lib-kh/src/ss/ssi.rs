//! The equivariant Rasmussen invariant `(s̲, s̄)` for a strongly invertible knot,
//! obtained from the `H`-divisibilities of the two equivariant Lee classes in
//! `KhI` over `𝔽₂[H]` (§3 of the reference).
//!
//! Reference:
//! - T. Sano, "Involutive Khovanov homology and equivariant knots",
//!   Algebr. Geom. Topol. 25 (2025), 5059–5111.
//!   <https://doi.org/10.2140/agt.2025.25.5059>, <https://arxiv.org/abs/2404.08568>
use itertools::Itertools;
use num_traits::Zero;
use log::info;

use yui_core::MathType;
use yui_core::num::FF2;
use yui_link::{InvLink, Link};

use crate::kh::KhComplex;
use crate::khi::{KhIChain, KhIComplex, KhIHomology};
use crate::tng::builder::SymBuildConfig;
use crate::util::FastPoly;
use super::util::div_vec;
use super::util::{assert_homogeneous, canon_q_deg, max_solvable_level, SsVersion};

type F = FF2;
type P = FastPoly<'H', F>;

/// `ssi` over `𝔽₂[H]`, via the current default pipeline ([`SsVersion::V2`]).
pub fn ssi_invariant(l: &InvLink, reduced: bool) -> (i32, i32) {
    ssi_invariant_with(l, reduced, SymBuildConfig::default(), None, SsVersion::V2)
}

/// The same, with the build configuration, the guessed s-value and the pipeline given explicitly.
pub fn ssi_invariant_with(l: &InvLink, reduced: bool, config: SymBuildConfig, expected: Option<isize>, ver: SsVersion) -> (i32, i32) {
    assert!(l.is_knot());
    assert_h_range(&config);

    let bot = -(l.n_signed_crossings().1 as isize);
    let config = SymBuildConfig {
        h_range: Some(config.h_range.clone().unwrap_or(bot ..= 1)),
        ..config
    };

    info!("compute ssi ({ver:?}) over {}.", P::math_symbol());

    let w = l.writhe();
    let r = l.seifert_circles().len() as i32;
    let (d0, d1) = match ver {
        SsVersion::V1 => ssi_divisibility_v1(l, reduced, config),
        SsVersion::V2 => ssi_divisibility_v2(l, reduced, config, expected),
    };

    let ss0 = 2 * d0 + w - r + 1;
    let ss1 = 2 * d1 + w - r + 1;

    info!("w = {w}, r = {r}, d0 = {d0}, d1 = {d1}.");
    info!("ssi = ({ss0}, {ss1}).");

    (ss0, ss1)
}

// The two equivariant Lee classes live at h = 0 and h = 1, so a window missing either computes
// nothing. `None` = the default (bottom..=1).
fn assert_h_range(config: &SymBuildConfig) {
    assert!(
        config.h_range.as_ref().is_none_or(|r| r.contains(&0) && r.contains(&1)),
        "the build h-range must contain 0 and 1, got {:?}", config.h_range
    );
}

fn ssi_divisibility_v1(l: &InvLink, reduced: bool, config: SymBuildConfig) -> (i32, i32) {
    let r = if reduced { 1 } else { 2 };
    let c = P::variable();
    let t = P::zero();

    let kh = KhIHomology::new_with_config(l, &c, &t, reduced, config);

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
        div_vec(&v.subvec(0..r), &c).expect("invalid divisibility.")
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

// initial high-q build cut at `q = expected − 1` — only the high end is cut, so `d_H` (mod torsion)
// is unchanged. A too-tight window (the reduction lifts the canon representative above it) is
// widened by 2 and retried. `None` = full (un-windowed) build.
fn ssi_divisibility_v2(l: &InvLink, reduced: bool, config: SymBuildConfig, expected: Option<isize>) -> (i32, i32) {
    let q0 = canon_q_deg(l.writhe(), l.seifert_circles().len(), reduced);

    let Some(s) = expected else {
        return divisibility_in_window(l, reduced, config, q0)
            .expect("the full (un-windowed) build must not truncate the canon cycle");
    };

    let q_hi0 = s - 1; // the top divisibility generator for the guessed s sits at q = s − 1.
    let mut q_hi = q_hi0;
    loop {
        info!("q-window: q0 = {q0}, cut above {q_hi}.");
        let cfg = SymBuildConfig { q_range: Some((isize::MIN + 1) ..= q_hi), ..config.clone() };

        if let Some(ds) = divisibility_in_window(l, reduced, cfg, q0) {
            return ds;
        }
        q_hi += 2; // one divisibility level = q-degree 2 (deg H = −2); raise the cut minimally.
        info!("canon cycle above window; raising cut to {q_hi}.");
        assert!(q_hi <= q_hi0 + 128, "q-window widening runaway");
    }
}

fn divisibility_in_window(l: &InvLink, reduced: bool, config: SymBuildConfig, q0: isize) -> Option<(i32, i32)> {
    let r = if reduced { 1 } else { 2 };
    let h = P::variable();
    let t = P::zero();

    // Build over `(a−1)..=1`: the solves only involve `C[−1] → C[0] → C[1]`, so nothing above 1
    // is built; the cheap low degrees stay (truncating the bottom makes the build frontier dense).
    // Only KhI degrees `−1..=1` are converted to the raw complex — the rest is never expanded.
    let requested = config.h_range.clone().unwrap_or(-(Link::MAX_CROSSING as isize) ..= 1);
    let range = KhComplex::<P>::clamp_h_range(l.inner(), reduced, requested);
    let (a, b) = (*range.start(), *range.end());
    let config = SymBuildConfig { h_range: Some((a - 1)..=b), ..config };
    let kc = KhIComplex::<P>::new_windowed(l, &h, &t, reduced, config, -1..=1);

    let zs = kc.canon_cycles(); // sorted by h-degree: B classes at 0, then Q classes at 1
    assert_eq!(zs.len(), 2 * r);

    // index gives the h-degree (0 for the first `r`, else 1) — robust to a cycle truncated to zero.
    let mut ds = [vec![], vec![]];
    for (i, z) in zs.iter().enumerate() {
        let i0 = if i < r { 0 } else { 1 };
        ds[i0 as usize].push(solvable_level(&kc, z, i0, q0)?);
    }

    let (d0, d1) = if reduced {
        (ds[0][0], ds[1][0])
    } else {
        assert_eq!(ds[0][0], ds[0][1]);
        assert_eq!(ds[1][0], ds[1][1]);
        (ds[0][0], ds[1][0])
    };

    assert!(d0 <= d1); // s̲ ≤ s̄.

    Some((d0, d1))
}

// Hand the shared level-climb the pieces it needs out of the KhI complex.
fn solvable_level(kc: &KhIComplex<P>, z: &KhIChain<P>, at: isize, q0: isize) -> Option<i32> {
    // an empty cycle means the window cut the whole representative away — signal widen.
    if z.is_zero() {
        return None;
    }
    assert_homogeneous(z, |x| kc.q_deg_of(x), q0);

    let q_degs = |i: isize| kc.inner()[i].raw_generators().iter().map(|x| kc.q_deg_of(x)).collect_vec();
    let d = kc.inner().d_matrix(at - 1);
    let v = kc.inner()[at].vectorize(z);

    max_solvable_level(&d, &v, &q_degs(at - 1), &q_degs(at), q0, at)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ssi_invariant_uses_the_default_config_and_v2() {
        let l = InvLink::test_data("3_1");
        assert_eq!(
            ssi_invariant(&l, false),
            ssi_invariant_with(&l, false, SymBuildConfig::default(), None, SsVersion::V2)
        );
    }

    #[test]
    fn test_unknot_pos_twist() {
        let l = InvLink::test_data("unknot_r_twist");
        let ssi = ssi_invariant_with(&l, false, SymBuildConfig::default(), None, SsVersion::V1);
        assert_eq!(ssi.0, 0);
        assert_eq!(ssi.1, 0);
    }

    #[test]
    fn test_unknot_neg_twist() {
        let l = InvLink::test_data("unknot_l_twist");
        let ssi = ssi_invariant_with(&l, false, SymBuildConfig::default(), None, SsVersion::V1);
        assert_eq!(ssi.0, 0);
        assert_eq!(ssi.1, 0);
    }

    #[test]
    fn test_unknot_neg_twist2() {
        let l = InvLink::test_data("unknot_l_twist2");
        let ssi = ssi_invariant_with(&l, false, SymBuildConfig::default(), None, SsVersion::V1);
        assert_eq!(ssi.0, 0);
        assert_eq!(ssi.1, 0);
    }

    // diagrams come from `test_data`, not the data dir, so the tests need no external resources.
    fn test(name: &str, ver: SsVersion, reduced: bool, expected: (i32, i32)) {
        let l = InvLink::test_data(name);
        let ssi = ssi_invariant_with(&l, reduced, SymBuildConfig::default(), None, ver);
        assert_eq!(ssi, expected);
    }

    macro_rules! test {
        ($test:ident, $name:literal, $expected:expr) => {
            mod $test {
                use super::*;

                #[test]
                fn v1() {
                    test($name, SsVersion::V1, false, $expected);
                }

                #[test]
                fn v2() {
                    test($name, SsVersion::V2, false, $expected);
                }

                #[test]
                fn v1_red() {
                    test($name, SsVersion::V1, true, $expected);
                }

                #[test]
                fn v2_red() {
                    test($name, SsVersion::V2, true, $expected);
                }
            }
        }
    }
    
    test!(k3_1, "3_1", (2, 2));
    test!(k4_1a, "4_1a", (0, 0));
    test!(k4_1b, "4_1b", (0, 0));
    test!(k5_1, "5_1", (4, 4));
    test!(k5_2a, "5_2a", (2, 2));
    test!(k5_2b, "5_2b", (2, 2));
    test!(k6_1a, "6_1a", (0, 0));
    test!(k6_1b, "6_1b", (0, 0));
    test!(k6_2a, "6_2a", (2, 2));
    test!(k6_2b, "6_2b", (2, 2));
    test!(k6_3a, "6_3a", (0, 0));
    test!(k6_3b, "6_3b", (0, 0));
    test!(k7_1, "7_1", (6, 6));
    test!(k7_2a, "7_2a", (2, 2));
    test!(k7_2b, "7_2b", (2, 2));
    test!(k7_3a, "7_3a", (4, 4));
    test!(k7_3b, "7_3b", (4, 4));
    test!(k7_4a, "7_4a", (2, 2));
    test!(k7_4b, "7_4b", (2, 2));
    test!(k7_5a, "7_5a", (4, 4));
    test!(k7_5b, "7_5b", (4, 4));
    test!(k7_6a, "7_6a", (2, 2));
    test!(k7_6b, "7_6b", (2, 2));
    test!(k7_7a, "7_7a", (0, 0));
    test!(k7_7b, "7_7b", (0, 0));

    // The (s̲, s̄) ≠ (s, s) cases of [Sano, InvKh II, Prop. 1.3] — everything above has s̲ = s̄, so
    // these are what keep a bug that collapses ssi to (s, s) from passing. 8_21's two inversion
    // classes differ, which is ssi seeing the involution and not just the knot.
    test!(k8_21a, "8_21a", (2, 2));
    test!(k8_21b, "8_21b", (2, 4));
    test!(k9_46a, "9_46a", (-2, 0));
    test!(k9_46b, "9_46b", (-2, 0));
}
