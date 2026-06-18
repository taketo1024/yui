//! The Rasmussen-type slice-torus invariants `ss̃_c(K) = 2·d_c(D) + w(D) − r(D) + 1`,
//! where `d_c` is the `c`-divisibility of the (unreduced or reduced) Lee class —
//! selected by the `reduced` argument to [`ss_invariant`].
//!
//! Reference:
//! - T. Sano and K. Sato, "A family of slice-torus invariants from the divisibility of Lee classes",
//!   Topol. Appl. 357 (2024), 109059.
//!   <https://doi.org/10.1016/j.topol.2024.109059>, <https://arxiv.org/abs/2211.02494>

use itertools::Itertools;
use log::info;
use num_traits::Zero;
use yui_link::Link;
use yui_core::{EucRing, EucRingOps};

use crate::util::calc::div_vec;
use crate::kh::KhHomology;

pub fn ss_invariant<R>(l: &Link, c: &R, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    assert!(!c.is_zero());
    assert!(!c.is_unit());
    assert!(l.is_knot());

    info!("compute ss, c = {c} ({}).", std::any::type_name::<R>());

    let w = l.writhe();
    let r = l.seifert_circles().len() as i32;
    let d = compute_div(l, c, reduced);
    let ss = 2 * d + w - r + 1;

    info!("d = {d}, w = {w}, r = {r}.");
    info!("ss = {ss} (c = {c}, {}).", if reduced { "reduced" } else { "unreduced" } );

    ss
}

fn compute_div<R>(l: &Link, c: &R, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    let r = if reduced { 1 } else { 2 };

    // bottom..=0: building the cheap low degrees and truncating only at the top is faster than
    // the doubly-truncated `0..=0` slice (which widens to the dense `-1..=1`). Builder clamps the start.
    let kh = KhHomology::new_partial(l, c, &R::zero(), reduced, Some(isize::MIN + 1 ..= 0));

    assert_eq!(kh[0].rank(), r);
    info!("Kh[0]: {}", kh[0]);
    
    let zs = kh.canon_cycles();

    assert_eq!(zs.len(), r);
    for z in zs.iter() {
        assert!(!z.is_zero());
        assert_eq!(z.homogeneous_value(|x| kh.h_deg_of(x)), Some(0));
    }

    let ds = zs.iter().enumerate().map(|(i, z)| {
        let v = kh[0].vectorize_euc(z);
        info!("a[{i}] in Kh[0]: ({})", v.clone().into_dense().iter().join(","));
        v
    }).map(|v| 
        div_vec(&v.subvec(0..r), c).expect("invalid divisibility.")
    ).collect_vec();

    assert!(ds.iter().all_equal());

    ds[0]
}

#[cfg(test)]
mod tests {
    use yui_link::Link;
    use super::*;

    #[test]
    fn test_unknot() { 
        let l = Link::unknot();
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 0);
        assert_eq!(ss_invariant(&l, &c, true ), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 0);
    }

    #[test]
    fn test_unknot_rm1() { 
        let l = Link::test_data("unknot_l_twist");
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 0);
        assert_eq!(ss_invariant(&l, &c, true ), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 0);
    }

    #[test]
    fn test_unknot_rm1_neg() { 
        let l = Link::test_data("unknot_r_twist");
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 0);
        assert_eq!(ss_invariant(&l, &c, true ), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 0);
    }

    #[test]
    fn test_3_1() { 
        let l = Link::test_data("3_1");
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), -2);
        assert_eq!(ss_invariant(&l, &c, true ), -2);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 2);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 2);
    }

    #[test]
    fn test_4_1() { 
        let l = Link::test_data("4_1");
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 0);
        assert_eq!(ss_invariant(&l, &c, true ), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 0);
    }

    #[test]
    fn test_5_1() { 
        let l = Link::test_data("5_1");
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), -4);
        assert_eq!(ss_invariant(&l, &c, true ), -4);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 4);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 4);
    }

    #[test]
    fn test_5_2() { 
        let l = Link::test_data("5_2");
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), -2);
        assert_eq!(ss_invariant(&l, &c, true ), -2);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 2);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 2);
    }

    #[test]
    fn test_6_1() { 
        let l = Link::test_data("6_1");
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 0);
        assert_eq!(ss_invariant(&l, &c, true ), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 0);
    }

    #[test]
    fn test_6_2() { 
        let l = Link::test_data("6_2");
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), -2);
        assert_eq!(ss_invariant(&l, &c, true ), -2);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 2);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 2);
    }

    #[test]
    fn test_6_3() { 
        let l = Link::test_data("6_3");
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 0);
        assert_eq!(ss_invariant(&l, &c, true ), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 0);
    }

    #[test]
    fn test_7_1() { 
        let l = Link::test_data("7_1");
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), -6);
        assert_eq!(ss_invariant(&l, &c, true ), -6);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 6);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 6);
    }

    #[test]
    fn test_7_2() { 
        let l = Link::test_data("7_2");
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), -2);
        assert_eq!(ss_invariant(&l, &c, true ), -2);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 2);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 2);
    }

    #[test]
    fn test_7_3() { 
        let l = Link::test_data("7_3");
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 4);
        assert_eq!(ss_invariant(&l, &c, true ), 4);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), -4);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), -4);
    }

    #[test]
    fn test_8_19() { 
        let l = Link::test_data("8_19");
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 6);
        assert_eq!(ss_invariant(&l, &c, true ), 6);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), -6);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), -6);
    }

    #[test]
    #[ignore]
    fn test_k14_c2() { 
        let l = Link::test_data("14n_19265");
        let c = 2_i64;
        
        assert_eq!(ss_invariant(&l, &c, false), -2);
        assert_eq!(ss_invariant(&l, &c, true ), -2);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 2);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 2);
    }

    #[test]
    #[ignore]
    fn test_k14_c3() { 
        let l = Link::test_data("14n_19265");
        let c = 3_i64;
        
        assert_eq!(ss_invariant(&l, &c, false), 0);
        assert_eq!(ss_invariant(&l, &c, true ), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 0);
    }
}