//! The Rasmussen-type slice-torus invariants `ss̃_c(K) = 2·d_c(D) + w(D) − r(D) + 1`,
//! where `d_c` is the `c`-divisibility of the (unreduced or reduced) Lee class —
//! selected by the `reduced` argument to [`ss_invariant`].
//!
//! Reference:
//! - T. Sano and K. Sato, "A family of slice-torus invariants from the divisibility of Lee classes",
//!   Topol. Appl. 357 (2024), 109059.
//!   <https://doi.org/10.1016/j.topol.2024.109059>, <https://arxiv.org/abs/2211.02494>
//!
//! Theorem 2 there identifies `ss̃_H` over `F[H]` with the Rasmussen invariant `s` over the
//! field `F` — the specialization [`s_invariant`] computes.


use itertools::Itertools;
use log::info;
use num_traits::Zero;
use yui_link::Link;
use yui_core::{EucRing, EucRingOps};

use yui_core::{Field, FieldOps, MathType};

use crate::kh::{KhChain, KhComplex, KhHomology};
use crate::tng::builder::BuildConfig;
use crate::util::FastPoly;
use super::util::{assert_homogeneous, canon_q_deg, div_vec, max_solvable_level, SsVersion};

type P<F> = FastPoly<'H', F>;

pub fn ss_invariant<R>(l: &Link, c: &R, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    ss_invariant_with(l, c, reduced, BuildConfig::default())
}

/// The same, with the build configuration given explicitly.
pub fn ss_invariant_with<R>(l: &Link, c: &R, reduced: bool, config: BuildConfig) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    assert!(!c.is_zero());
    assert!(!c.is_unit());
    assert!(l.is_knot());
    assert_h_range(&config);

    let bot = -(l.n_signed_crossings().1 as isize);
    let config = BuildConfig {
        h_range: Some(config.h_range.clone().unwrap_or(bot ..= 0)),
        ..config
    };

    info!("compute ss, c = {c} ({}).", std::any::type_name::<R>());

    let w = l.writhe();
    let r = l.seifert_circles().len() as i32;
    let d = ss_divisibility(l, c, reduced, config);
    let ss = 2 * d + w - r + 1;

    info!("d = {d}, w = {w}, r = {r}.");
    info!("ss = {ss} (c = {c}, {}).", if reduced { "reduced" } else { "unreduced" } );

    ss
}

// The Lee class lives at h = 0, so a build window that excludes it computes nothing.
fn assert_h_range(config: &BuildConfig) {
    assert!(
        config.h_range.as_ref().is_none_or(|r| r.contains(&0)),
        "the build h-range must contain 0, got {:?}", config.h_range
    );
}

fn ss_divisibility<R>(l: &Link, c: &R, reduced: bool, config: BuildConfig) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    let r = if reduced { 1 } else { 2 };

    let kh = KhHomology::new_with_config(l, c, &R::zero(), reduced, config);

    assert_eq!(kh[0].rank(), r);
    info!("Kh[0]: {}", kh[0]);
    
    let zs = kh.canon_cycles();

    assert_eq!(zs.len(), r);
    for z in zs.iter() {
        assert!(!z.is_zero());
        assert_eq!(z.homogeneous_value(|x| kh.h_deg_of(x)), Some(0));
    }

    let ds = zs.iter().enumerate().map(|(i, z)| {
        let v = kh[0].vectorize_euc(z).subvec(0..r);
        info!("a[{i}] in Kh[0]: ({})", v.clone().into_dense().iter().join(", "));
        div_vec(&v, c).expect("invalid divisibility.")
    }).collect_vec();

    assert!(ds.iter().all_equal());

    ds[0]
}

/// The Rasmussen invariant `s` over a field `F`: the specialization of `ss̃_c` to `R = F[H]`,
/// `c = H` (Theorem 2 of the reference), via the `H = 1` pipeline.
pub fn s_invariant<F>(l: &Link, reduced: bool) -> i32
where F: Field, for<'x> &'x F: FieldOps<F> {
    s_invariant_with::<F>(l, reduced, BuildConfig::default(), SsVersion::default())
}

/// The same, with the build configuration and the pipeline given explicitly.
pub fn s_invariant_with<F>(l: &Link, reduced: bool, config: BuildConfig, ver: SsVersion) -> i32
where F: Field, for<'x> &'x F: FieldOps<F> {
    // V1 is exactly `ss̃_H` over `F[H]`; only V2 needs the specialized route.
    if ver == SsVersion::V1 {
        let c = P::<F>::variable();
        return ss_invariant_with::<P<F>>(l, &c, reduced, config);
    }
    
    assert!(l.is_knot());
    assert_h_range(&config);

    info!("compute s ({ver:?}) over {}.", P::<F>::math_symbol());

    let w = l.writhe();
    let r = l.seifert_circles().len() as i32;
    let d = s_divisibility_v2::<F>(l, reduced, config);
    let s = 2 * d + w - r + 1;

    info!("d = {d}, w = {w}, r = {r}.");
    info!("s = {s} ({}).", if reduced { "reduced" } else { "unreduced" });

    s
}

// The `H`-divisibility of the Lee class by `H = 1` solves — no homology, no basis tracking.
fn s_divisibility_v2<F>(l: &Link, reduced: bool, config: BuildConfig) -> i32
where F: Field, for<'x> &'x F: FieldOps<F> {
    let n = if reduced { 1 } else { 2 };
    let q0 = canon_q_deg(l.writhe(), l.seifert_circles().len(), reduced);
    let (h, t) = (P::<F>::variable(), P::<F>::zero());

    // the solve involves `C[-1] -> C[0]` only. Build one degree wider on each side: the extra
    // degree below feeds the differential, and the one above is needed to close the circles at 0
    // (the canon cycle is not evaluable if the build is cut there).
    let requested = config.h_range.clone().unwrap_or(-(Link::MAX_CROSSING as isize) ..= 0);
    let range = KhComplex::<P<F>>::clamp_h_range(l, reduced, requested);
    let (a, b) = (*range.start(), *range.end());
    let config = BuildConfig { h_range: Some((a - 1)..=(b + 1)), ..config };
    let kc = KhComplex::<P<F>>::new_with_config(l, &h, &t, reduced, config);

    let zs = kc.canon_cycles();
    assert_eq!(zs.len(), n);

    let ds = zs.iter().map(|z|
        solvable_level(&kc, z, 0, q0).expect("the full build must not truncate the canon cycle")
    ).collect_vec();

    assert!(ds.iter().all_equal());

    ds[0]
}

// Hand the shared level-climb the pieces it needs out of the Khovanov complex.
fn solvable_level<F>(kc: &KhComplex<P<F>>, z: &KhChain<P<F>>, at: isize, q0: isize) -> Option<i32>
where F: Field, for<'x> &'x F: FieldOps<F> {
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
    use yui_core::num::{FF, FF2, Ratio};
    use yui_link::Link;
    use super::*;

    // `s` over F[H] must agree with `ss̃_H` computed the V1 way (Theorem 2 of the reference), and
    // with the classical values pinned by the `ss` tests below.
    #[test]
    fn s_matches_ss_over_f2() {
        type F = FF2;
        for (name, expected) in [("3_1", 2), ("4_1", 0), ("5_1", 4), ("5_2", 2), ("6_2", 2), ("6_3", 0)] {
            let l = Link::test_data(name);
            for reduced in [false, true] {
                let v1 = s_invariant_with::<F>(&l, reduced, BuildConfig::default(), SsVersion::V1);
                let v2 = s_invariant::<F>(&l, reduced);
                assert_eq!(v1, expected, "{name} V1, reduced = {reduced}");
                assert_eq!(v2, expected, "{name} V2, reduced = {reduced}");
            }
        }
    }

    #[test]
    fn s_over_q_and_f3() {
        // the invariant is a priori field-dependent; for these knots the three fields agree.
        for (name, expected) in [("3_1", 2), ("4_1", 0), ("8_19", 6)] {
            let l = Link::test_data(name);
            assert_eq!(s_invariant::<Ratio<i64>>(&l, false), expected, "{name} over Q");
            assert_eq!(s_invariant::<FF<3>>(&l, false), expected, "{name} over F3");
        }
    }

    #[test]
    #[should_panic(expected = "h-range must contain 0")]
    fn s_rejects_an_h_range_missing_zero() {
        let l = Link::test_data("3_1");
        let config = BuildConfig { h_range: Some(1..=2), ..Default::default() };
        let _ = s_invariant_with::<FF2>(&l, false, config, SsVersion::V2);
    }

    #[test]
    fn s_of_the_mirror_negates() {
        type F = FF2;
        for name in ["3_1", "5_2", "8_19"] {
            let l = Link::test_data(name);
            assert_eq!(s_invariant::<F>(&l.mirror(), false), -s_invariant::<F>(&l, false), "{name}");
        }
    }

    // Every case below: reduced and unreduced agree, and the mirror negates.
    fn check<R>(l: &Link, c: &R, expected: i32)
    where R: EucRing, for<'x> &'x R: EucRingOps<R> {
        assert_eq!(ss_invariant(l, c, false), expected);
        assert_eq!(ss_invariant(l, c, true ), expected);
        assert_eq!(ss_invariant(&l.mirror(), c, false), -expected);
        assert_eq!(ss_invariant(&l.mirror(), c, true ), -expected);
    }

    // `c = 2` over Z, the classical Rasmussen invariant.
    macro_rules! test_c2 {
        ($test:ident, $name:literal, $expected:expr) => {
            #[test]
            fn $test() {
                check(&Link::test_data($name), &2, $expected);
            }
        };
    }

    #[test]
    fn unknot() {
        check(&Link::unknot(), &2, 0);
    }

    test_c2!(unknot_rm1,     "unknot_l_twist", 0);
    test_c2!(unknot_rm1_neg, "unknot_r_twist", 0);
    test_c2!(k3_1,  "3_1",  2);
    test_c2!(k4_1,  "4_1",  0);
    test_c2!(k5_1,  "5_1",  4);
    test_c2!(k5_2,  "5_2",  2);
    test_c2!(k6_1,  "6_1",  0);
    test_c2!(k6_2,  "6_2",  2);
    test_c2!(k6_3,  "6_3",  0);
    test_c2!(k7_1,  "7_1",  6);
    test_c2!(k7_2,  "7_2",  2);
    test_c2!(k7_3,  "7_3",  4);
    test_c2!(k8_19, "8_19", 6);

    // The invariant genuinely depends on the coefficients, so one computation per choice is the
    // point here — the mirror/reduced invariance is already pinned by the knots above.
    //
    // 14n_19265: `ss` over Z differs for c = 2 and c = 3; `s` over F2 differs from Q and F3.
    #[test]
    fn k14_ring_dependence() {
        let l = Link::test_data("14n_19265");

        assert_eq!(ss_invariant(&l, &2_i64, false), -2, "c = 2");
        assert_eq!(ss_invariant(&l, &3_i64, false),  0, "c = 3");

        assert_eq!(s_invariant::<Ratio<i64>>(&l, false),  0, "Q");
        assert_eq!(s_invariant::<FF2>(&l, false),        -2, "F2");
        assert_eq!(s_invariant::<FF<3>>(&l, false),       0, "F3");
    }

    // The F3 counterpart of 14n_19265: here Q and F2 agree and F3 differs. Values from the
    // computations behind [Sano-Sato, Topol. Appl. 357 (2024)].
    #[test]
    #[ignore = "slow in debug (~15s): two 18-crossing knots over three fields"]
    fn k18_ring_dependence() {
        for (name, s_q, s_f3) in [("18nh_05566876", 2, 0), ("18nh_37144251", -2, 0)] {
            let l = Link::test_data(name);
            assert_eq!(s_invariant::<Ratio<i64>>(&l, false), s_q,  "{name} over Q");
            assert_eq!(s_invariant::<FF2>(&l, false),        s_q,  "{name} over F2");
            assert_eq!(s_invariant::<FF<3>>(&l, false),      s_f3, "{name} over F3");
        }
    }
}
