// "A family of slice-torus invariants from the divisibility of reduced Lee classes"
// T. Sano, K. Sato
// https://arxiv.org/abs/2211.02494

use core::panic;

use itertools::Itertools;
use log::info;
use yui_link::Link;
use yui_homology::{RModStr, ChainComplex};
use yui_homology::utils::HomologyCalc;
use yui_core::{EucRing, EucRingOps};
use yui_matrix::sparse::*;
use crate::KhComplex;

pub fn ss_invariant<R>(l: &Link, c: &R, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    assert!(!c.is_zero());
    assert!(!c.is_unit());

    info!("compute ss, n = {}, c = {} ({}).", l.crossing_num(), c, std::any::type_name::<R>());

    let (a0, a1, vs) = compute_d(l, c, reduced);
    let (h, p, _) = HomologyCalc::calculate_with_trans(a0, a1);

    info!("homology: {h}");
    
    let r = h.rank();
    let p = p.submat_rows(0..r).to_owned();

    let vs = vs.into_iter().enumerate().map(|(i, v)| { 
        let v = &p * v;
        info!("a[{i}] = {v}");
        v
    }).collect_vec();

    let v = &vs[0];
    let Some(d) = div_vec(&v, &c) else { 
        panic!("invalid divisibility for v = {}, c = {}", v, c)
    };

    let w = l.writhe();
    let r = l.seifert_circles().len() as i32;
    let ss = 2 * d + w - r + 1;

    info!("d = {d}, w = {w}, r = {r}.");
    info!("ss = {ss} (c = {c}, {}).", if reduced { "reduced" } else { "unreduced" } );

    ss
}

fn compute_d<R>(l: &Link, h: &R, reduced: bool) -> (SpMat<R>, SpMat<R>, Vec<SpVec<R>>)
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    let c = KhComplex::<R>::new(l, h, &R::zero(), reduced);

    let a0 = c.d_matrix(-1);
    let a1 = c.d_matrix(0);

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
    use super::*;
    use yui_link::{Link, Edge};

    fn test<I>(pd_code: I, value: i32)
    where I: IntoIterator<Item = [Edge; 4]> { 
        let l = Link::from_pd_code(pd_code);
        assert_eq!(ss_invariant(&l, &2, true), value);
        assert_eq!(ss_invariant(&l.mirror(), &2, true), -value);
    }

    #[test]
    fn test_3_1() { 
        test([[1,4,2,5],[3,6,4,1],[5,2,6,3]], -2);
    }

    #[test]
    fn test_4_1() { 
        test([[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]], 0);
    }

    #[test]
    fn test_5_1() { 
        test([[1,6,2,7],[3,8,4,9],[5,10,6,1],[7,2,8,3],[9,4,10,5]], -4);
    }

    #[test]
    fn test_5_2() { 
        test([[1,4,2,5],[3,8,4,9],[5,10,6,1],[9,6,10,7],[7,2,8,3]], -2);
    }

    #[test]
    fn test_6_1() { 
        test([[1,4,2,5],[7,10,8,11],[3,9,4,8],[9,3,10,2],[5,12,6,1],[11,6,12,7]], 0);
    }

    #[test]
    fn test_6_2() { 
        test([[1,4,2,5],[5,10,6,11],[3,9,4,8],[9,3,10,2],[7,12,8,1],[11,6,12,7]], -2);
    }

    #[test]
    fn test_6_3() { 
        test([[4,2,5,1],[8,4,9,3],[12,9,1,10],[10,5,11,6],[6,11,7,12],[2,8,3,7]], 0);
    }
}