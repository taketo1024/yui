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
use yui_core::poly::Poly;
use yui_link::{InvLink, Link};

use crate::tng::builder::SymBuildConfig;
use crate::util::calc::div_vec;
use crate::khi::KhIHomology;
use super::ssi_h1::ssi_divisibility_v2;

type P = Poly<'H', FF2>;

/// The `ssi` computation pipeline.
/// - `V1`: full bigraded `KhIHomology`; simple, memory-heavy (`config`/`expected` are ignored).
/// - `V2`: `H = 1` specialization — per-level q-truncated `𝔽₂` solves, no homology or basis
///   tracking. `expected` (the guessed s-value) seeds the high-q build cut; `None` = full build.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SsiVersion {
    V1,
    V2
}

/// `ssi` over `𝔽₂[H]`, via the current default pipeline ([`SsiVersion::V2`]).
pub fn ssi_invariant(l: &InvLink, reduced: bool, config: SymBuildConfig, expected: Option<isize>) -> (i32, i32) {
    ssi_invariant_ver(l, reduced, config, expected, SsiVersion::V2)
}

pub fn ssi_invariant_ver(l: &InvLink, reduced: bool, config: SymBuildConfig, expected: Option<isize>, ver: SsiVersion) -> (i32, i32) {
    assert!(l.is_knot());

    info!("compute ssi ({ver:?}) over {}.", P::math_symbol());

    let w = l.writhe();
    let r = l.seifert_circles().len() as i32;
    let (d0, d1) = match ver {
        SsiVersion::V1 => ssi_divisibility_v1(l, reduced),
        SsiVersion::V2 => ssi_divisibility_v2(l, reduced, config, expected),
    };

    let ss0 = 2 * d0 + w - r + 1;
    let ss1 = 2 * d1 + w - r + 1;

    info!("w = {w}, r = {r}, d0 = {d0}, d1 = {d1}.");
    info!("ssi = ({ss0}, {ss1}).");

    (ss0, ss1)
}

fn ssi_divisibility_v1(l: &InvLink, reduced: bool) -> (i32, i32) {
    let r = if reduced { 1 } else { 2 };
    let c = P::variable();
    let t = P::zero();

    // bottom..=1: building the cheap low degrees and truncating only at the top is faster than
    // the doubly-truncated `0..=1` slice (which widens to the dense `-1..=2`). Builder clamps the start.
    let kh = KhIHomology::new_partial(l, &c, &t, reduced, Some(-(Link::MAX_CROSSING as isize) ..= 1));

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unknot_pos_twist() {
        let l = InvLink::test_data("unknot_r_twist");
        let ssi = ssi_invariant_ver(&l, false, SymBuildConfig::default(), None, SsiVersion::V1);
        assert_eq!(ssi.0, 0);
        assert_eq!(ssi.1, 0);
    }

    #[test]
    fn test_unknot_neg_twist() {
        let l = InvLink::test_data("unknot_l_twist");
        let ssi = ssi_invariant_ver(&l, false, SymBuildConfig::default(), None, SsiVersion::V1);
        assert_eq!(ssi.0, 0);
        assert_eq!(ssi.1, 0);
    }

    #[test]
    fn test_unknot_neg_twist2() {
        let l = InvLink::test_data("unknot_l_twist2");
        let ssi = ssi_invariant_ver(&l, false, SymBuildConfig::default(), None, SsiVersion::V1);
        assert_eq!(ssi.0, 0);
        assert_eq!(ssi.1, 0);
    }

    fn test(name: &str, ver: SsiVersion, reduced: bool, expected: (i32, i32)) -> Result<(), Box<dyn std::error::Error>> {
        let l = InvLink::load(name)?;

        let ssi = ssi_invariant_ver(&l, reduced, SymBuildConfig::default(), None, ver);
        assert_eq!(ssi, expected);

        Ok(())
    }

    macro_rules! test {
        ($test:ident, $name:literal, $expected:expr) => {
            mod $test {
                use super::*;

                #[test]
                fn v1() -> Result<(), Box<dyn std::error::Error>> {
                    test($name, SsiVersion::V1, false, $expected)
                }

                #[test]
                fn v2() -> Result<(), Box<dyn std::error::Error>> {
                    test($name, SsiVersion::V2, false, $expected)
                }

                #[test]
                fn v1_red() -> Result<(), Box<dyn std::error::Error>> {
                    test($name, SsiVersion::V1, true, $expected)
                }

                #[test]
                fn v2_red() -> Result<(), Box<dyn std::error::Error>> {
                    test($name, SsiVersion::V2, true, $expected)
                }
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

}