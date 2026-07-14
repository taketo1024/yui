//! `ssi` via specialization at `H = 1`: the divisibility `d_H` of each equivariant Lee class
//! is read off from a family of linear systems `d(w) = z` over `F`, posed in q-truncations of
//! the Bar-Natan specialization of `CKhI(D)`. No homology, SNF or basis tracking is involved.
//!
//! Reference: `research/h1-divisibility.tex`. Solvability at level `m` (window `q < q₀ + 2m`)
//! is equivalent to `d_H([z]) ≥ m`, the `H`-divisibility modulo torsion — exactly the quantity
//! defining `ssi`.

use itertools::Itertools;
use num_traits::Zero;
use log::info;

use yui_core::{Field, FieldOps, MathType};
use yui_core::poly::Poly;
use yui_matrix::MatTrait;
use yui_matrix::sparse::{SpMat, SpVec};
use yui_matrix::sparse::pluq::solve_pluq;
use yui_link::{InvLink, Link};

use crate::kh::KhComplex;
use crate::tng::builder::SymBuildConfig;
use crate::khi::{KhIChain, KhIComplex};

type PolyH<F> = Poly<'H', F>;

/// `ssi` from the `H = 1` linear systems. The complex is built over `F[H]` (q-homogeneous)
/// exactly as in `ssi_invariant`; only the solves happen over `F`. Char-2 only (the cone).
pub fn ssi_invariant_h1<F>(l: &InvLink, reduced: bool, config: SymBuildConfig) -> (i32, i32)
where F: Field, for<'x> &'x F: FieldOps<F> {
    assert!(l.is_knot());

    info!("compute ssi via H=1 solves over {}.", PolyH::<F>::math_symbol());

    let w = l.writhe();
    let r = l.seifert_circles().len() as i32;
    let (d0, d1) = ssi_divisibility_h1::<F>(l, reduced, config);

    let ss0 = 2 * d0 + w - r + 1;
    let ss1 = 2 * d1 + w - r + 1;

    info!("w = {w}, r = {r}, d0 = {d0}, d1 = {d1}.");
    info!("ssi = ({ss0}, {ss1}).");

    (ss0, ss1)
}

/// `ssi_invariant_h1` with the build's high-q cut set from a target divisibility `max_div`: only
/// generators of q-degree `≤ q₀ + 2·(max_div + 2)` are built (`q₀` = the canon cycle's q-degree).
/// Since only the high end is cut, `d_H` (mod torsion) is unchanged, so the result is exact provided
/// the true divisibility is `≤ max_div + 1` (otherwise the solve reaches the built ceiling — raise
/// `max_div`). This is the memory lever for large diagrams (Wh-doubles).
pub fn ssi_invariant_h1_windowed<F>(l: &InvLink, reduced: bool, config: SymBuildConfig, max_div: i32) -> (i32, i32)
where F: Field, for<'x> &'x F: FieldOps<F> {
    let q0 = canon_q_deg::<F>(l, reduced);
    let q_hi = q0 + 2 * (max_div as isize + 2);
    info!("q-window: q0 = {q0}, cut above {q_hi} (max_div = {max_div}).");

    // only the high-q cut; keep everything below (the full quotient needs the low generators for
    // the torsion corrections that make this `d_H` mod torsion, not honest divisibility).
    let config = SymBuildConfig { q_range: Some((isize::MIN + 1) ..= q_hi), ..config };
    ssi_invariant_h1::<F>(l, reduced, config)
}

// The canon cycle's homogeneous q-degree `q₀`, from the Seifert-state cycle alone (no complex build).
fn canon_q_deg<F>(l: &InvLink, reduced: bool) -> isize
where F: Field, for<'x> &'x F: FieldOps<F> {
    let deg_shift = KhComplex::<PolyH<F>>::deg_shift_for(l.inner(), reduced);
    let (a, b) = (PolyH::<F>::zero(), PolyH::<F>::variable());
    let zs = KhComplex::<PolyH<F>>::make_canon_cycles(l.inner(), &a, &b, reduced);

    zs.iter().flat_map(|z|
        z.iter().map(|(x, c)| deg_shift.1 + x.rel_q_deg() - 2 * (c.lead_deg() as isize))
    ).min().expect("empty canon cycle")
}

fn ssi_divisibility_h1<F>(l: &InvLink, reduced: bool, config: SymBuildConfig) -> (i32, i32)
where F: Field, for<'x> &'x F: FieldOps<F> {
    let r = if reduced { 1 } else { 2 };
    let h = PolyH::<F>::variable();
    let t = PolyH::<F>::zero();

    // same windowing as `ssi_divisibility`: cover the canon degrees 0, 1, one wider on both ends.
    let requested = config.h_range.clone().unwrap_or(-(Link::MAX_CROSSING as isize) ..= 1);
    let range = KhComplex::<PolyH<F>>::clamp_h_range(l.inner(), reduced, requested);
    let (a, b) = (*range.start(), *range.end());
    assert!(a <= 0 && b >= 1, "ssi h-range must include 0 and 1, got {a}..={b}");
    let config = SymBuildConfig { h_range: Some((a - 1)..=(b + 1)), ..config };
    let kc = KhIComplex::<PolyH<F>>::new_with_config(l, &h, &t, reduced, config);

    let zs = kc.canon_cycles(); // sorted by h-degree: B classes at 0, then Q classes at 1
    assert_eq!(zs.len(), 2 * r);
    for (i, z) in zs.iter().enumerate() {
        let expected = if i < r { 0 } else { 1 };
        assert!(!z.is_zero());
        assert_eq!(z.homogeneous_value(|x| kc.h_deg_of(x)), Some(expected));
    }

    let ds = [0, 1].map(|i0|
        zs.iter()
            .filter(|z| kc.h_deg_of_chain(z) == i0)
            .map(|z| max_solvable_level::<F>(&kc, z, i0))
            .collect_vec()
    );

    let (d0, d1) = if reduced {
        (ds[0][0], ds[1][0])
    } else {
        assert_eq!(ds[0][0], ds[0][1]);
        assert_eq!(ds[1][0], ds[1][1]);
        (ds[0][0], ds[1][0])
    };

    assert!(d0 <= d1); // s̲ ≤ s̄.

    (d0, d1)
}

// Climb `m` while the level-`m` system `d(w) = z` in `C̄ / F_{q₀+2m}` solves; the first failure
// gives `d_H = m − 1`. Level 0 is trivially solvable (`z̄` truncated to zero), so the climb starts at 1.
fn max_solvable_level<F>(kc: &KhIComplex<PolyH<F>>, z: &KhIChain<PolyH<F>>, i0: isize) -> i32
where F: Field, for<'x> &'x F: FieldOps<F> {
    // the homogeneous degree of `z`: coefficient `H^k` shifts a term's degree by `−2k`.
    // (`q_deg_of_chain` ignores coefficients — wrong here when the transported `z` is `H`-divisible.)
    let q0 = z.iter().map(|(x, a)| kc.q_deg_of(x) - 2 * (a.lead_deg() as isize)).min()
        .expect("canon cycle is zero");
    assert_homogeneous(kc, z, q0);

    let src_q = q_degs(kc, i0 - 1);
    let tgt_q = q_degs(kc, i0);
    let d = kc.inner().d_matrix(i0 - 1);
    let v = kc.inner()[i0].vectorize(z);

    // beyond this the window contains everything and the system stops changing.
    let q_top = Iterator::chain(src_q.iter(), tgt_q.iter()).max().copied().unwrap_or(q0);

    let mut m = 0;
    loop {
        let q_hi = q0 + 2 * (m + 1);
        assert!(q_hi <= q_top + 2, "untruncated z̄ solved as a boundary — the Lee class cannot vanish");

        let (a, y) = truncated_system(&d, &v, &src_q, &tgt_q, q_hi);
        info!("h = {i0}, level {}: solve system of size {:?}.", m + 1, a.shape());

        if solve_pluq(&a, &y).is_some() {
            m += 1;
        } else {
            info!("h = {i0}: d = {m}.");
            return m as i32;
        }
    }
}

// The level system: rows/cols of `d^{i0-1}` restricted to the q-window `q < q_hi`
// (the full quotient `C̄ / F_{q_hi}`), entries evaluated at `H = 1`.
fn truncated_system<F>(
    d: &SpMat<PolyH<F>>, v: &SpVec<PolyH<F>>,
    src_q: &[isize], tgt_q: &[isize], q_hi: isize,
) -> (SpMat<F>, SpVec<F>)
where F: Field, for<'x> &'x F: FieldOps<F> {
    let in_win = |q: isize| q < q_hi;
    let (n_rows, row_map) = window_index(tgt_q, &in_win);
    let (n_cols, col_map) = window_index(src_q, &in_win);

    let a = SpMat::from_entries((n_rows, n_cols), d.iter_nz().filter_map(|(i, j, p)| {
        let i1 = row_map[i]?;
        let j1 = col_map[j]?;
        debug_assert!(
            p.nterms() == 1 && 2 * (p.lead_deg() as isize) == tgt_q[i] - src_q[j],
            "the differential is not q-homogeneous"
        );
        let e = eval_at_one(p);
        (!e.is_zero()).then_some((i1, j1, e))
    }));

    let y = SpVec::from_entries(n_rows, v.iter_nz().filter_map(|(i, p)| {
        let i1 = row_map[i]?;
        let e = eval_at_one(p);
        (!e.is_zero()).then_some((i1, e))
    }));

    (a, y)
}

fn q_degs<F>(kc: &KhIComplex<PolyH<F>>, i: isize) -> Vec<isize>
where F: Field, for<'x> &'x F: FieldOps<F> {
    kc.inner()[i].raw_generators().iter().map(|x| kc.q_deg_of(x)).collect()
}

fn window_index(qs: &[isize], in_win: impl Fn(isize) -> bool) -> (usize, Vec<Option<usize>>) {
    let mut n = 0;
    let map = qs.iter().map(|&q|
        in_win(q).then(|| {
            let i = n;
            n += 1;
            i
        })
    ).collect();
    (n, map)
}

fn eval_at_one<F>(p: &PolyH<F>) -> F
where F: Field, for<'x> &'x F: FieldOps<F> {
    F::sum(p.iter().map(|(_, a)| a.clone()))
}

// every coefficient of `z` must be a monomial `λ·H^k` with `q(x) − 2k = q₀`.
fn assert_homogeneous<F>(kc: &KhIComplex<PolyH<F>>, z: &KhIChain<PolyH<F>>, q0: isize)
where F: Field, for<'x> &'x F: FieldOps<F> {
    let bad = z.iter().filter(|(x, a)|
        a.nterms() != 1 || kc.q_deg_of(*x) - 2 * (a.lead_deg() as isize) != q0
    ).map(|(x, a)|
        format!("  q(x) = {}, coeff = {} (nterms {}, lead_deg {})", kc.q_deg_of(x), a, a.nterms(), a.lead_deg())
    ).collect_vec();

    assert!(
        bad.is_empty(),
        "canon cycle is not q-homogeneous at q0 = {q0}:\n{}", bad.iter().join("\n")
    );
}

#[cfg(test)]
mod tests {
    use yui_core::num::FF2;
    use super::*;

    type F = FF2;

    fn ssi_h1(l: &InvLink, reduced: bool) -> (i32, i32) {
        ssi_invariant_h1::<F>(l, reduced, SymBuildConfig::default())
    }

    macro_rules! test {
        ($test:ident, $name:literal, $expected:expr) => {
            #[test]
            fn $test() -> Result<(), Box<dyn std::error::Error>> {
                let l = InvLink::load($name)?;
                assert_eq!(ssi_h1(&l, false), $expected);
                Ok(())
            }
        }
    }

    test!(k3_1, "3_1", (2, 2));
    test!(k4_1, "4_1", (0, 0));
    test!(k5_1, "5_1", (4, 4));
    test!(k6_2a, "6_2a", (2, 2));
    test!(k6_3, "6_3", (0, 0));
    test!(k7_6a, "7_6a", (-2, -2));

    // the q-windowed build (high-q cut at the true divisibility) reproduces the full ssi.
    #[test]
    fn windowed_matches_full() -> Result<(), Box<dyn std::error::Error>> {
        for name in ["3_1", "4_1", "5_1", "5_2a", "6_2a", "6_3", "7_6a"] {
            let l = InvLink::load(name)?;
            let full = ssi_invariant_h1::<F>(&l, false, SymBuildConfig::default());

            // tight window at the true divisibility d = (ss − w + r − 1) / 2 (max over the two classes).
            let (w, r) = (l.writhe(), l.seifert_circles().len() as i32);
            let d_max = [(full.0 - w + r - 1) / 2, (full.1 - w + r - 1) / 2].into_iter().max().unwrap();
            let win = ssi_invariant_h1_windowed::<F>(&l, false, SymBuildConfig::default(), d_max);

            assert_eq!(win, full, "{name} (d_max = {d_max})");
        }
        Ok(())
    }

    #[test]
    fn k3_1_m() {
        let l = InvLink::test_data("3_1").mirror();
        assert_eq!(ssi_h1(&l, false), (-2, -2));
    }

    #[test]
    fn k3_1_red() {
        let l = InvLink::test_data("3_1");
        assert_eq!(ssi_h1(&l, true), (2, 2));
    }

    #[test]
    fn k9_46() {
        let l = InvLink::from_symmetric_pd_code(
            [[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]
        );
        assert_eq!(ssi_h1(&l, false), (0, 2));
    }

    // h1 must agree with the divisibility pipeline (`ssi_invariant`), incl. chunked builds.
    #[test]
    fn matches_ssi_invariant_chunked() -> Result<(), Box<dyn std::error::Error>> {
        use yui_core::poly::Poly;
        use crate::tng::builder::CutOption;
        use crate::khi::ssi_invariant;

        type P = Poly<'H', F>;
        let c = P::variable();

        for name in ["3_1", "4_1", "6_1a"] {
            let l = InvLink::load(name)?;
            let config = SymBuildConfig { cut: CutOption::Auto(3), ..Default::default() };
            let expected = ssi_invariant(&l, &c, false, config.clone());
            let h1 = ssi_invariant_h1::<F>(&l, false, config);
            assert_eq!(h1, expected, "{name}");
        }
        Ok(())
    }
}
