// "A family of slice-torus invariants from the divisibility of reduced Lee classes"
// T. Sano, K. Sato
// https://arxiv.org/abs/2211.02494

use core::panic;

use itertools::Itertools;
use log::info;
use yui_homology::{ChainComplexTrait, RModStr};
use yui_link::Link;
use yui_homology::utils::HomologyCalc;
use yui_core::{EucRing, EucRingOps};
use yui_matrix::sparse::*;

use crate::KhComplex;

pub fn ss_invariant<R>(l: &Link, c: &R, n: usize, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    ss_invariant_v2(l, c, n, reduced)
}

pub fn ss_invariant_v1<R>(l: &Link, c: &R, n: usize, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    compute_ss(l, c, n, reduced, 1)
}

pub fn ss_invariant_v2<R>(l: &Link, c: &R, n: usize, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    compute_ss(l, c, n, reduced, 2)
}

fn compute_ss<R>(l: &Link, c: &R, n: usize, reduced: bool, ver: usize) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    assert!(!c.is_zero());
    assert!(!c.is_unit());
    assert!(n >= 1);
    assert!(ver == 2);

    info!("compute ss, c = {c} ({}), n = {n}.", std::any::type_name::<R>());

    let h = (0..n).fold(R::one(), |p, _| p * c); // h = c^n
    let (a0, a1, vs) = match ver { 
        2 => prepare_v2(l, &h, reduced),
        _ => panic!()
    };

    let kh = HomologyCalc::calculate(&a0, &a1, true);

    info!("homology: {}", kh.math_symbol());
    
    let r = kh.rank();
    let t = kh.trans().unwrap();

    let vs = vs.into_iter().enumerate().map(|(i, v)| { 
        let v = t.forward(&v).subvec(0..r);
        info!("a[{i}] = {v}");
        v
    }).collect_vec();

    let v = &vs[0];
    let Some(d) = div_vec(&v, &c) else { 
        panic!("invalid divisibility for v = {}, c = {}", v, c)
    };

    let w = l.writhe() as i32;
    let r = l.seifert_circles().len() as i32;
    let n = n as i32;

    let ss = 2 * d + n * (w - r + 1);

    info!("d = {d}, w = {w}, r = {r}.");
    info!("ss = {ss} (c = {c}, {}).", if reduced { "reduced" } else { "unreduced" } );

    ss
}

fn prepare_v2<R>(l: &Link, h: &R, reduced: bool) -> (SpMat<R>, SpMat<R>, Vec<SpVec<R>>)
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    let c = KhComplex::new_v2(l, h, &R::zero(), reduced);

    let a0 = c.d_matrix(-1);
    let a1 = c.d_matrix( 0);

    let vs = c.canon_cycles().iter().map(|z| 
        c[0].vectorize(&z)
    ).collect_vec();

    (a0, a1, vs)
}

fn div_vec<R>(v: &SpVec<R>, c: &R) -> Option<i32>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    v.iter().filter_map(|(_, a)|
        div(a, c)
    ).min()
}

fn div<R>(a: &R, c: &R) -> Option<i32>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    if a.is_zero() { return None }

    let mut a = a.clone();
    let mut k = 0;

    while (&a % c).is_zero() { 
        a /= c;
        k += 1;
    }

    Some(k)
}

#[cfg(test)]
mod tests {
    use cartesian::cartesian;
    use yui_link::Link;
    use super::*;

    fn test(l: &Link, c: i32, value: i32) { 
        let ver = 2;
        for &r in cartesian!([true, false].iter()) { 
            assert_eq!(compute_ss(&l, &c, 1, r, ver), value);
            assert_eq!(compute_ss(&l.mirror(), &c, 1, r, ver), -value);
        }
    }

    #[test]
    fn test_unknot() { 
        let l = Link::unknot();
        test(&l, 2, 0);
    }

    #[test]
    fn test_unknot_rm1() { 
        let l = Link::from_pd_code([[0,0,1,1]]);
        test(&l, 2, 0);
    }

    #[test]
    fn test_unknot_rm1_neg() { 
        let l = Link::from_pd_code([[0,1,1,0]]);
        test(&l, 2, 0);
    }

    #[test]
    fn test_3_1() { 
        let l = Link::from_pd_code([[1,4,2,5],[3,6,4,1],[5,2,6,3]]);
        test(&l, 2, -2);
    }

    #[test]
    fn test_4_1() { 
        let l = Link::from_pd_code([[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]]);
        test(&l, 2, 0);
    }

    #[test]
    fn test_5_1() { 
        let l = Link::from_pd_code([[1,6,2,7],[3,8,4,9],[5,10,6,1],[7,2,8,3],[9,4,10,5]]);
        test(&l, 2, -4);
    }

    #[test]
    fn test_5_2() { 
        let l = Link::from_pd_code([[1,4,2,5],[3,8,4,9],[5,10,6,1],[9,6,10,7],[7,2,8,3]]);
        test(&l, 2, -2);
    }

    #[test]
    fn test_6_1() { 
        let l = Link::from_pd_code([[1,4,2,5],[7,10,8,11],[3,9,4,8],[9,3,10,2],[5,12,6,1],[11,6,12,7]]);
        test(&l, 2, 0);
    }

    #[test]
    fn test_6_2() { 
        let l = Link::from_pd_code([[1,4,2,5],[5,10,6,11],[3,9,4,8],[9,3,10,2],[7,12,8,1],[11,6,12,7]]);
        test(&l, 2, -2);
    }

    #[test]
    fn test_6_3() { 
        let l = Link::from_pd_code([[4,2,5,1],[8,4,9,3],[12,9,1,10],[10,5,11,6],[6,11,7,12],[2,8,3,7]]);
        test(&l, 2, 0);
    }

    #[test]
    fn test_7_1() { 
        let l = Link::from_pd_code([[1,8,2,9],[3,10,4,11],[5,12,6,13],[7,14,8,1],[9,2,10,3],[11,4,12,5],[13,6,14,7]]);
        test(&l, 2, -6);
    }

    #[test]
    fn test_7_2() { 
        let l = Link::from_pd_code([[1,4,2,5],[3,10,4,11],[5,14,6,1],[7,12,8,13],[11,8,12,9],[13,6,14,7],[9,2,10,3]]);
        test(&l, 2, -2);
    }

    #[test]
    fn test_7_3() { 
        let l = Link::from_pd_code([[6,2,7,1],[10,4,11,3],[14,8,1,7],[8,14,9,13],[12,6,13,5],[2,10,3,9],[4,12,5,11]]);
        test(&l, 2, 4);
    }

    #[test]
    fn test_8_19() { 
        let l = Link::from_pd_code([[4,2,5,1],[8,4,9,3],[9,15,10,14],[5,13,6,12],[13,7,14,6],[11,1,12,16],[15,11,16,10],[2,8,3,7]]);
        test(&l, 2, 6);
    }

    #[test]
    #[cfg(not(debug_assertions))]
    fn test_k14() { 
        let l = Link::from_pd_code([[1,19,2,18],[19,1,20,28],[20,13,21,14],[12,17,13,18],[16,21,17,22],[5,15,6,14],[15,5,16,4],[6,27,7,28],[2,7,3,8],[26,3,27,4],[25,23,26,22],[11,9,12,8],[23,10,24,11],[9,24,10,25]]);
        assert_eq!(ss_invariant_v2(&l, &2, 1, true), -2);
        assert_eq!(ss_invariant_v2(&l, &3, 1, true), 0);
    }
}