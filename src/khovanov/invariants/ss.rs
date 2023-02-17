// "A family of slice-torus invariants from the divisibility of reduced Lee classes"
// T. Sano, K. Sato
// https://arxiv.org/abs/2211.02494

use log::info;
use sprs::{CsMat, CsVec};
use crate::khovanov::complex::KhComplex;
use crate::links::Link;
use crate::math::homology::base::{RModStr, GenericRModStr};
use crate::math::homology::complex::ChainComplex;
use crate::math::homology::utils::homology_calc::HomologyCalc;
use crate::math::homology::utils::reducer::ChainReducer;
use crate::math::matrix::sparse::CsMatExt;
use crate::math::traits::{EucRing, EucRingOps};

pub fn ss_invariant<R>(l: &Link, c: &R, reduced: bool) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    assert!(!c.is_zero());
    assert!(!c.is_unit());

    if reduced && l.writhe() < 0 { 
        info!("switch to mirror, w = {}.", l.writhe());
        return -ss_invariant(&l.mirror(), c, reduced)
    }

    info!("compute ss, n = {}, c = {} ({}).", l.crossing_num(), c, std::any::type_name::<R>());

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

    let d = div_vec(&v, &c);
    let w = l.writhe();
    let r = l.seifert_circles().len() as i32;
    let s = 2 * d + w - r + 1;

    info!("d = {d}, w = {s}, r = {r}, s = {s}");

    s
}

fn homology_with<R>(a0: CsMat<R>, a1: CsMat<R>, a2: CsMat<R>, a3: CsMat<R>, v: CsVec<R>) -> (GenericRModStr<R>, CsVec<R>)
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    let vs = vec![v];
    let (n0, n3) = (a0.cols(), a3.rows());

    info!("(k = -2)");
    let ( _, a0, a1, _, _) 
        = ChainReducer::reduce_with(CsMat::zero((n0, 0)), a0, a1, vec![], vec![]);

    info!("(k = -1)");
    let ( _, a1, a2, _, vs) 
        = ChainReducer::reduce_with(a0, a1, a2, vec![], vs);

    info!("(k = 0)");
    let (a1, a2, a3, vs, _) 
        = ChainReducer::reduce_with(a1, a2, a3, vs, vec![]);

    info!("(k = 1)");
    let (a2,  _,  _)
        = ChainReducer::reduce(a2, a3, CsMat::zero((0, n3)));

    let (h, p, _) = HomologyCalc::calculate_with_trans(a1, a2);

    info!("homology: {h}");
    
    let r = h.rank();
    let p = p.submatrix(0..r, 0..p.cols());

    let v = &vs[0];
    let w = &p * v;

    info!("canon-cycle: {}", w.to_dense());

    (h, w)
}

fn div_vec<R>(v: &CsVec<R>, c: &R) -> i32 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    v.iter().map(|(_, a)|
        div(a, c)
    ).min().unwrap_or(i32::MAX)
}

fn div<R>(a: &R, c: &R) -> i32 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    if a.is_zero() { return i32::MAX }

    let mut a = a.clone();
    let mut k = 0;

    while (&a % c).is_zero() { 
        a /= c;
        k += 1;
    }

    k
}