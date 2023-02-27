// "A family of slice-torus invariants from the divisibility of reduced Lee classes"
// T. Sano, K. Sato
// https://arxiv.org/abs/2211.02494

use core::panic;

use log::info;
use sprs::CsVec;
use crate::khovanov::complex::KhComplex;
use crate::links::Link;
use crate::math::homology::base::{RModStr, GenericRModStr};
use crate::math::homology::complex::ChainComplex;
use crate::math::homology::utils::homology_calc::HomologyCalc;
use crate::math::homology::utils::reducer::ChainReducer;
use crate::math::matrix::sparse::*;
use crate::math::traits::{EucRing, EucRingOps};

pub fn ss_invariant<R>(l: &Link, c: &R, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    assert!(!c.is_zero());
    assert!(!c.is_unit());

    if reduced && l.writhe() < 0 { 
        info!("compute ss (mirror), n = {}, c = {} ({}).", l.crossing_num(), c, std::any::type_name::<R>());
        -(compute_ss(&l.mirror(), c, reduced))
    } else { 
        info!("compute ss, n = {}, c = {} ({}).", l.crossing_num(), c, std::any::type_name::<R>());
        compute_ss(l, c, reduced)
    }
}

fn compute_ss<R>(l: &Link, c: &R, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    let cpx = KhComplex::<R>::new(l.clone(), c.clone(), R::zero(), reduced);
    let z = cpx.canon_cycle();
    let v = cpx[0].vectorize(&z);

    let (_, v) = homology_with(
        cpx.d_matrix(-2), 
        cpx.d_matrix(-1), 
        cpx.d_matrix(0), 
        cpx.d_matrix(1), 
        v
    );

    let Some(d) = div_vec(&v, &c) else { 
        panic!("invalid divisibility for v = {}, c = {}", v.to_dense(), c)
    };

    let w = l.writhe();
    let r = l.seifert_circles().len() as i32;
    let s = 2 * d + w - r + 1;

    info!("d = {d}, w = {s}, r = {r}, s = {s}");

    s
}

fn homology_with<R>(a0: SpMat<R>, a1: SpMat<R>, a2: SpMat<R>, a3: SpMat<R>, v: CsVec<R>) -> (GenericRModStr<R>, CsVec<R>)
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    let zero = |m, n| SpMat::<R>::zero((m, n));

    let vs = vec![v];
    let (n0, n3) = (a0.cols(), a3.rows());

    info!("(k = -2)");
    let ( _, a0, a1, _, _) 
        = ChainReducer::reduce_with(zero(n0, 0), a0, a1, vec![], vec![]);

    info!("(k = -1)");
    let ( _, a1, a2, _, vs) 
        = ChainReducer::reduce_with(a0, a1, a2, vec![], vs);

    info!("(k = 0)");
    let (a1, a2, a3, vs, _) 
        = ChainReducer::reduce_with(a1, a2, a3, vs, vec![]);

    info!("(k = 1)");
    let (a2,  _,  _)
        = ChainReducer::reduce(a2, a3, zero(0, n3));

    let (h, p, _) = HomologyCalc::calculate_with_trans(a1, a2);

    info!("homology: {h}");
    
    let r = h.rank();
    let p = p.submatrix_rows(0..r).to_owned();

    let v = &vs[0];
    let w = &p * v;

    info!("canon-cycle: {}", w.to_dense());

    (h, w)
}

fn div_vec<R>(v: &CsVec<R>, c: &R) -> Option<i32>
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
    use indexmap::IndexMap;

    use crate::links::Link;
    use crate::utils::collections::map;
    use super::*;

    fn test(name: &str) { 
        let values: IndexMap<_, _> = map!{ 
            "3_1" => -2,
            "4_1" =>  0,
            "5_1" => -4,
            "5_2" => -2,
            "6_1" => -0,
            "6_2" => -2,
            "6_3" =>  0,
            "7_1" => -6,
            "7_2" => -2,
            "7_3" =>  4, // mirrored
            "7_4" =>  2, // mirrored
            "7_5" => -4,
            "7_6" => -2,
            "7_7" =>  0
        };

        let l = Link::load(name).unwrap();
        let s = values[name];

        assert_eq!(compute_ss(&l, &2, true), s);
        assert_eq!(compute_ss(&l.mirror(), &2, true), -s);
    }

    #[test]
    fn test_3_1() { 
        test("3_1");
    }

    #[test]
    fn test_4_1() { 
        test("4_1");
    }

    #[test]
    fn test_5_1() { 
        test("5_1");
    }

    #[test]
    fn test_5_2() { 
        test("3_1");
    }

    #[test]
    fn test_6_1() { 
        test("6_1");
    }

    #[test]
    fn test_6_2() { 
        test("6_2");
    }

    #[test]
    fn test_6_3() { 
        test("6_3");
    }

    #[test]
    fn test_7_1() { 
        test("7_1");
    }

    #[test]
    fn test_7_2() { 
        test("7_2");
    }

    #[test]
    fn test_7_3() { 
        test("7_3");
    }

    #[test]
    fn test_7_4() { 
        test("7_4");
    }

    #[test]
    fn test_7_5() { 
        test("7_5");
    }

    #[test]
    fn test_7_6() { 
        test("7_6");
    }

    #[test]
    fn test_7_7() { 
        test("7_7");
    }
}