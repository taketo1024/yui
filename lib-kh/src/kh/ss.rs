// "A family of slice-torus invariants from the divisibility of reduced Lee classes"
// T. Sano, K. Sato
// https://arxiv.org/abs/2211.02494

use itertools::Itertools;
use log::info;
use num_traits::Zero;
use yui_homology::SummandTrait;
use yui_link::Link;
use yui::{EucRing, EucRingOps};

use crate::misc::div_vec;
use crate::kh::{KhChainExt, KhHomology};

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

    let kh = KhHomology::new(l, c, &R::zero(), reduced).truncated(0..=0);
    // let kh = KhHomology::new_partial(l, c, &R::zero(), reduced, Some(0..=0));

    assert_eq!(kh[0].rank(), r);
    info!("Kh[0]: {}", kh[0]);
    
    let zs = kh.canon_cycles();

    assert_eq!(zs.len(), r);
    assert!(zs.iter().all(|z| !z.is_zero()));
    assert!(zs.iter().all(|z| z.h_deg() == 0));    

    let ds = zs.iter().enumerate().map(|(i, z)| {
        let v = kh[0].vectorize_euc(z);
        info!("a[{i}] in Kh[0]: ({})", v.to_dense().iter().join(","));
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
        let l = Link::from_pd_code([[0,0,1,1]]);
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 0);
        assert_eq!(ss_invariant(&l, &c, true ), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 0);
    }

    #[test]
    fn test_unknot_rm1_neg() { 
        let l = Link::from_pd_code([[0,1,1,0]]);
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 0);
        assert_eq!(ss_invariant(&l, &c, true ), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 0);
    }

    #[test]
    fn test_3_1() { 
        let l = Link::from_pd_code([[1,4,2,5],[3,6,4,1],[5,2,6,3]]);
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), -2);
        assert_eq!(ss_invariant(&l, &c, true ), -2);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 2);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 2);
    }

    #[test]
    fn test_4_1() { 
        let l = Link::from_pd_code([[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]]);
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 0);
        assert_eq!(ss_invariant(&l, &c, true ), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 0);
    }

    #[test]
    fn test_5_1() { 
        let l = Link::from_pd_code([[1,6,2,7],[3,8,4,9],[5,10,6,1],[7,2,8,3],[9,4,10,5]]);
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), -4);
        assert_eq!(ss_invariant(&l, &c, true ), -4);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 4);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 4);
    }

    #[test]
    fn test_5_2() { 
        let l = Link::from_pd_code([[1,4,2,5],[3,8,4,9],[5,10,6,1],[9,6,10,7],[7,2,8,3]]);
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), -2);
        assert_eq!(ss_invariant(&l, &c, true ), -2);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 2);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 2);
    }

    #[test]
    fn test_6_1() { 
        let l = Link::from_pd_code([[1,4,2,5],[7,10,8,11],[3,9,4,8],[9,3,10,2],[5,12,6,1],[11,6,12,7]]);
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 0);
        assert_eq!(ss_invariant(&l, &c, true ), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 0);
    }

    #[test]
    fn test_6_2() { 
        let l = Link::from_pd_code([[1,4,2,5],[5,10,6,11],[3,9,4,8],[9,3,10,2],[7,12,8,1],[11,6,12,7]]);
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), -2);
        assert_eq!(ss_invariant(&l, &c, true ), -2);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 2);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 2);
    }

    #[test]
    fn test_6_3() { 
        let l = Link::from_pd_code([[4,2,5,1],[8,4,9,3],[12,9,1,10],[10,5,11,6],[6,11,7,12],[2,8,3,7]]);
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 0);
        assert_eq!(ss_invariant(&l, &c, true ), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 0);
    }

    #[test]
    fn test_7_1() { 
        let l = Link::from_pd_code([[1,8,2,9],[3,10,4,11],[5,12,6,13],[7,14,8,1],[9,2,10,3],[11,4,12,5],[13,6,14,7]]);
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), -6);
        assert_eq!(ss_invariant(&l, &c, true ), -6);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 6);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 6);
    }

    #[test]
    fn test_7_2() { 
        let l = Link::from_pd_code([[1,4,2,5],[3,10,4,11],[5,14,6,1],[7,12,8,13],[11,8,12,9],[13,6,14,7],[9,2,10,3]]);
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), -2);
        assert_eq!(ss_invariant(&l, &c, true ), -2);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 2);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 2);
    }

    #[test]
    fn test_7_3() { 
        let l = Link::from_pd_code([[6,2,7,1],[10,4,11,3],[14,8,1,7],[8,14,9,13],[12,6,13,5],[2,10,3,9],[4,12,5,11]]);
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 4);
        assert_eq!(ss_invariant(&l, &c, true ), 4);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), -4);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), -4);
    }

    #[test]
    fn test_8_19() { 
        let l = Link::from_pd_code([[4,2,5,1],[8,4,9,3],[9,15,10,14],[5,13,6,12],[13,7,14,6],[11,1,12,16],[15,11,16,10],[2,8,3,7]]);
        let c = 2;
        
        assert_eq!(ss_invariant(&l, &c, false), 6);
        assert_eq!(ss_invariant(&l, &c, true ), 6);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), -6);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), -6);
    }

    #[test]
    #[ignore]
    fn test_k14_c2() { 
        let l = Link::from_pd_code([[1,19,2,18],[19,1,20,28],[20,13,21,14],[12,17,13,18],[16,21,17,22],[5,15,6,14],[15,5,16,4],[6,27,7,28],[2,7,3,8],[26,3,27,4],[25,23,26,22],[11,9,12,8],[23,10,24,11],[9,24,10,25]]);
        let c = 2_i64;
        
        assert_eq!(ss_invariant(&l, &c, false), -2);
        assert_eq!(ss_invariant(&l, &c, true ), -2);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 2);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 2);
    }

    #[test]
    #[ignore]
    fn test_k14_c3() { 
        let l = Link::from_pd_code([[1,19,2,18],[19,1,20,28],[20,13,21,14],[12,17,13,18],[16,21,17,22],[5,15,6,14],[15,5,16,4],[6,27,7,28],[2,7,3,8],[26,3,27,4],[25,23,26,22],[11,9,12,8],[23,10,24,11],[9,24,10,25]]);
        let c = 3_i64;
        
        assert_eq!(ss_invariant(&l, &c, false), 0);
        assert_eq!(ss_invariant(&l, &c, true ), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, false), 0);
        assert_eq!(ss_invariant(&l.mirror(), &c, true ), 0);
    }
}