use std::collections::HashMap;
use itertools::Itertools;
use log::*;
use sprs::PermOwned;

use yui_matrix::sparse::*;
use yui_matrix::sparse::pivot::{PivotType, PivotCondition, perms_by_pivots, find_pivots};
use yui_matrix::sparse::schur::Schur;
use yui_matrix::sparse::triang::{solve_triangular_vec, TriangularType};
use yui::{Ring, RingOps};

use crate::generic::GenericChainComplexBase;
use crate::{GridDeg, ChainComplexTrait, GridTrait};

//       a0 = [x]      a1 = [a b]      a2 = [z w]
//            [y]           [c d]     
//  C[0] --------> C[1] ---------> C[2] -------> C[3]
//    |    [1 a⁻¹b] |             | [ 1     ]    |
//    |    [    1 ] |             | [-ca⁻¹ 1]    |
//    |             V             V              |    
//  C[0] --------> C[1] ---------> C[2] -------> C[3]
//    |     [0]     |    [a 0]    |    [0  w]    |
//    |     [y]     |    [0 s]    |              |
//    |             V             V              |
//  C[0] --------> C[1]'---------> C[2]'-------> C[3]
//          [y]           [s]           [w]

pub struct ChainReducer<I, R>
where 
    I: GridDeg,
    R: Ring, for<'x> &'x R: RingOps<R>,
{ 
    support: Vec<I>,
    d_deg: I,
    mats: HashMap<I, SpMat<R>>,
    trans: HashMap<I, Trans<R>>,
    vecs: HashMap<I, Vec<SpVec<R>>>
}

impl<I, R> ChainReducer<I, R>
where 
    I: GridDeg,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    pub fn reduce<C>(complex: &C, with_trans: bool) -> Self
    where C: GridTrait<I> + ChainComplexTrait<I, R = R> {
        let mut r = Self::from(complex, with_trans);
        r.reduce_all(false);
        r.reduce_all(true);
        r
    }

    pub fn from<C>(complex: &C, with_trans: bool) -> Self 
    where C: GridTrait<I> + ChainComplexTrait<I, R = R> {
        let support = complex.support();
        let d_deg = complex.d_deg();

        let mut reducer = Self::new(support, d_deg);

        for i in reducer.support.clone() { 
            for j in [i, i + d_deg] { 
                if !reducer.is_set(j) {
                    let d = complex.d_matrix(j);
                    reducer.set_matrix(j, d, with_trans);
                }
            }
        }

        reducer
    }

    pub fn new<Itr>(support: Itr, d_deg: I) -> Self
    where Itr: Iterator<Item = I> {
        let support = support.collect_vec();
        let mats = HashMap::new();
        let trans = HashMap::new();
        let vecs = HashMap::new();
        Self { support, d_deg, mats, trans, vecs }
    }

    pub fn support(&self) -> &[I] { 
        &self.support
    }

    pub fn matrix(&self, i: I) -> Option<&SpMat<R>> {
        self.mats.get(&i)
    }

    pub fn trans(&self, i: I) -> Option<&Trans<R>> {
        self.trans.get(&i)
    }

    pub fn vecs(&self, i: I) -> Option<&Vec<SpVec<R>>> { 
        self.vecs.get(&i)
    }

    pub fn trans_mut(&mut self, i: I) -> Option<&mut Trans<R>> {
        self.trans.get_mut(&i)
    }

    pub fn rank(&self, i: I) -> Option<usize> { 
        self.matrix(i).map(|d| d.ncols())
    }

    pub fn is_set(&self, i: I) -> bool { 
        self.mats.contains_key(&i)
    }

    pub fn is_done(&self) -> bool { 
        self.support.iter().all(|&i| 
            self.matrix(i).map(|d| d.is_zero()).unwrap_or(false)
        )
    }

    pub fn set_matrix(&mut self, i: I, d: SpMat<R>, with_trans: bool) {
        if with_trans { 
            let n = d.ncols();
            self.trans.insert(i, Trans::id(n));
        }
        self.mats.insert(i, d);
    }

    pub fn add_vec(&mut self, i: I, v: SpVec<R>) { 
        self.vecs.entry(i).or_default().push(v)
    }

    pub fn reduce_all(&mut self, deep: bool) { 
        if self.is_done() { 
            return
        }
        
        if deep { 
            info!("reduce all (deep)");
        } else {
            info!("reduce all (shallow)");
        }

        let support = self.support.clone();

        for &i in support.iter() { 
            self.reduce_at(i, deep);
        }
    }

    pub fn reduce_at(&mut self, i: I, deep: bool) { 
        let mut c = 1;
        loop { 
            let (piv_type, piv_cond) = self.preferred_strategy(i);
            let cont = self.reduce_at_spec(i, piv_type, piv_cond);

            if !deep || !cont { 
                break
            }

            c += 1;
            debug!("next itr: {c}");
        }
    }

    pub fn reduce_at_spec(&mut self, i: I, piv_type: PivotType, piv_cond: PivotCondition) -> bool { 
        let Some(a) = self.matrix(i) else { 
            panic!("not initialized at {i}");
        };

        if a.is_zero() { 
            return false;
        }

        info!("reduce C[{i}]: {:?} ..", a.shape());
        debug!("  nnz: {}", a.nnz());
        debug!("  density: {}", a.density());
        debug!("  mean-weight: {}", a.mean_weight());

        let (p, q, r) = pivots(a, piv_type, piv_cond);

        if r == 0 { 
            info!("  done.");
            return false;
        }

        info!("  found {r} pivots.");
        
        let a = a.permute(p.view(), q.view());

        let t = match piv_type { 
            PivotType::Rows => TriangularType::Upper,
            PivotType::Cols => TriangularType::Lower
        };

        let with_trans = 
            self.trans.contains_key(&i) || 
            self.trans.contains_key(&(i + self.d_deg));

        let sch = Schur::from_partial_triangular(t, &a, r, with_trans);
        let (s, t_src, t_tgt) = sch.disassemble();

        info!("  reduced C[{i}]: {:?} -> {:?}", a.shape(), s.shape());

        self.update_mats(i, &p, &q, r, s);

        if with_trans { 
            let t_src = t_src.unwrap();
            let t_tgt = t_tgt.unwrap();
            self.update_trans(i, &p, &q, t_src, t_tgt);
        }

        self.update_vecs(i, &a, &p, &q, r, t);

        true
    }

    pub fn preferred_strategy(&self, i: I) -> (PivotType, PivotCondition) { 
        let Some(a) = self.matrix(i) else { 
            panic!("not initialized at {i}");
        };

        // TODO improve
        let piv_type = PivotType::Cols;

        let piv_cond = if a.iter().any(|(_, _, r)| r.is_pm_one()) { 
            PivotCondition::One
        } else { 
            PivotCondition::AnyUnit
        };

        (piv_type, piv_cond)
    }

    fn update_trans(&mut self, i: I, p: &PermOwned, q: &PermOwned, t_src: Trans<R>, t_tgt: Trans<R>) {
        let (_, i1, i2) = self.deg_trip(i);
        
        if let Some(t1) = self.trans_mut(i1) { 
            t1.append_perm(q.view());
            t1.merge(t_src);
        }

        if let Some(t2) = self.trans_mut(i2) {
            t2.append_perm(p.view());
            t2.merge(t_tgt);
        }
    }

    fn update_mats(&mut self, i: I, p: &PermOwned, q: &PermOwned, r: usize, s: SpMat<R>) {
        let (m, n) = (p.dim(), q.dim());
        let (i0, i1, i2) = self.deg_trip(i);

        if let Some(a0) = self.matrix(i0) {
            assert_eq!(a0.nrows(), n);
            let a0 = reduce_mat_rows(a0, q, r);
            self.mats.insert(i0, a0);
        }

        self.mats.insert(i1, s);

        if let Some(a2) = self.matrix(i2) { 
            assert_eq!(a2.ncols(), m);
            let a2 = reduce_mat_cols(a2, p, r);
            self.mats.insert(i2, a2);
        }
    }

    fn update_vecs(&mut self, i: I, a: &SpMat<R>, p: &PermOwned, q: &PermOwned, r: usize, t: TriangularType) {
        let (m, n) = a.shape();
        let (_, i1, i2) = self.deg_trip(i);
        
        if let Some(vs) = self.vecs.get_mut(&i1) { 
            for v in vs.iter_mut() { 
                assert_eq!(v.dim(), n);

                let w = v.extract(n - r, |i| {
                    let i = q.at(i);
                    (r..n).contains(&i).then(|| i - r)
                });

                *v = w;
            }
        }

        if let Some(vs) = self.vecs.get_mut(&i2) { 
            debug!("update {} vecs in C[{i2}] ..", vs.len());

            let [a, _, c, _] = a.divide4((r, r));
            
            for v in vs.iter_mut() { 
                assert_eq!(v.dim(), m);

                let (x, y) = v.permute(p.view()).split(r);
                let ainvx = solve_triangular_vec(t, &a, &x);
                let w = y - &c * ainvx;

                *v = w;
            }
        }
    }

    fn deg_trip(&self, i: I) -> (I, I, I) { 
        let deg = self.d_deg;
        (i - deg, i, i + deg)
    }

    pub fn into_complex(self) -> GenericChainComplexBase<I, R> { 
        GenericChainComplexBase::generate(
            self.support.clone(), self.d_deg, |i| self.mats[&i].clone()
        )
    }
}

fn pivots<R>(a: &SpMat<R>, piv_type: PivotType, pivot_cond: PivotCondition) -> (PermOwned, PermOwned, usize) 
where R: Ring, for<'x> &'x R: RingOps<R> {
    let pivs = find_pivots(a, piv_type, pivot_cond);
    let (p, q) = perms_by_pivots(a, &pivs);
    let r = pivs.len();
    (p, q, r)
}

fn reduce_mat_rows<R>(a: &SpMat<R>, p: &PermOwned, r: usize) -> SpMat<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (m, n) = a.shape();
    a.extract((m - r, n), |i, j| { 
        let i = p.at(i);
        (r..m).contains(&i).then(|| (i - r, j))
    })
}

fn reduce_mat_cols<R>(a: &SpMat<R>, p: &PermOwned, r: usize) -> SpMat<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (m, n) = a.shape();
    a.extract((m, n - r), |i, j| { 
        let j = p.at(j);
        (r..n).contains(&j).then(|| (i, j - r))
    })
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use crate::generic::GenericChainComplex;
    use crate::SummandTrait;
    use super::*;

    #[test]
    fn zero() { 
        let c = GenericChainComplex::<i32>::zero();
        let r = ChainReducer::reduce(&c, false).into_complex();

        r.check_d_all();

        assert!(r[0].is_zero());
    }

    #[test]
    fn acyclic() { 
        let c = GenericChainComplex::<i32>::one_one(1);
        let r = ChainReducer::reduce(&c, false).into_complex();

        r.check_d_all();

        assert!(r[0].is_zero());
        assert!(r[1].is_zero());
    }

    #[test]
    fn tor() { 
        let c = GenericChainComplex::<i32>::one_one(2);
        let r = ChainReducer::reduce(&c, false).into_complex();

        r.check_d_all();

        assert_eq!(r[0].rank(), 1);
        assert_eq!(r[1].rank(), 1);
    }
    
    #[test]
    fn d3() {
        let c = GenericChainComplex::<i32>::d3();
        let r = ChainReducer::reduce(&c, false).into_complex();

        r.check_d_all();

        assert_eq!(r[0].rank(), 1);
        assert_eq!(r[1].rank(), 0);
        assert_eq!(r[2].rank(), 0);
        assert_eq!(r[3].rank(), 0);

        assert!(r.d_matrix(0).is_zero());
        assert!(r.d_matrix(1).is_zero());
        assert!(r.d_matrix(2).is_zero());
        assert!(r.d_matrix(3).is_zero());
    }

    #[test]
    fn s2() {
        let c = GenericChainComplex::<i32>::s2();
        let r = ChainReducer::reduce(&c, false).into_complex();

        r.check_d_all();

        assert_eq!(r[0].rank(), 1);
        assert_eq!(r[1].rank(), 0);
        assert_eq!(r[2].rank(), 1);

        assert!(r.d_matrix(0).is_zero());
        assert!(r.d_matrix(1).is_zero());
        assert!(r.d_matrix(2).is_zero());
    }

    #[test]
    fn t2() {
        let c = GenericChainComplex::<i32>::t2();
        let r = ChainReducer::reduce(&c, false).into_complex();

        r.check_d_all();

        assert_eq!(r[0].rank(), 1);
        assert_eq!(r[1].rank(), 2);
        assert_eq!(r[2].rank(), 1);

        assert!(r.d_matrix(0).is_zero());
        assert!(r.d_matrix(1).is_zero());
        assert!(r.d_matrix(2).is_zero());
    }

    #[test]
    fn rp2() {
        let c = GenericChainComplex::<i32>::rp2();
        let r = ChainReducer::reduce(&c, false).into_complex();

        r.check_d_all();

        assert_eq!(r[0].rank(), 1);
        assert_eq!(r[1].rank(), 1);
        assert_eq!(r[2].rank(), 1);

        assert!( r.d_matrix(0).is_zero());
        assert!( r.d_matrix(1).is_zero());
        assert!(!r.d_matrix(2).is_zero());

        let a = r.d_matrix(2).into_dense()[(0, 0)];
        assert!(a == 2 || a == -2);
    }

    #[test]
    fn s2_trans() {
        let c = GenericChainComplex::<i32>::s2();
        let r = ChainReducer::reduce(&c, true);

        let t0 = r.trans(0).unwrap();
        let t1 = r.trans(1).unwrap();
        let t2 = r.trans(2).unwrap();

        assert_eq!(t0.src_dim(), 4);
        assert_eq!(t1.src_dim(), 6);
        assert_eq!(t2.src_dim(), 4);

        assert_eq!(t0.tgt_dim(), 1);
        assert_eq!(t1.tgt_dim(), 0);
        assert_eq!(t2.tgt_dim(), 1);

        assert_eq!(t0.forward_mat() * t0.backward_mat(), SpMat::id(1));
        assert_eq!(t1.forward_mat() * t1.backward_mat(), SpMat::id(0));
        assert_eq!(t2.forward_mat() * t2.backward_mat(), SpMat::id(1));

        let v = SpVec::unit(1, 0);
        let w = t2.backward(&v);
        let z = c[2].devectorize(&w);
        
        assert!(c.d(2, &z).is_zero());
    }

    #[test]
    fn t2_trans() {
        let c = GenericChainComplex::<i32>::t2();
        let r = ChainReducer::reduce(&c, true);

        let t0 = r.trans(0).unwrap();
        let t1 = r.trans(1).unwrap();
        let t2 = r.trans(2).unwrap();

        assert_eq!(t0.src_dim(), 9);
        assert_eq!(t1.src_dim(), 27);
        assert_eq!(t2.src_dim(), 18);

        assert_eq!(t0.tgt_dim(), 1);
        assert_eq!(t1.tgt_dim(), 2);
        assert_eq!(t2.tgt_dim(), 1);

        assert_eq!(t0.forward_mat() * t0.backward_mat(), SpMat::id(1));
        assert_eq!(t1.forward_mat() * t1.backward_mat(), SpMat::id(2));
        assert_eq!(t2.forward_mat() * t2.backward_mat(), SpMat::id(1));

        let v = SpVec::unit(1, 0);
        let w = t2.backward(&v);
        let z = c[2].devectorize(&w);

        assert!(c.d(2, &z).is_zero());

        let a = SpVec::unit(2, 0);
        let a = t1.backward(&a);
        let a = c[1].devectorize(&a);
        
        let b = SpVec::unit(2, 1);
        let b = t1.backward(&b);
        let b = c[1].devectorize(&b);
        
        assert!(c.d(1, &a).is_zero());
        assert!(c.d(1, &b).is_zero());
    }

    #[test]
    fn rp2_trans() {
        let c = GenericChainComplex::<i32>::rp2();
        let r = ChainReducer::reduce(&c, true);

        let t0 = r.trans(0).unwrap();
        let t1 = r.trans(1).unwrap();
        let t2 = r.trans(2).unwrap();

        assert_eq!(t0.src_dim(), 6);
        assert_eq!(t1.src_dim(), 15);
        assert_eq!(t2.src_dim(), 10);

        assert_eq!(t0.tgt_dim(), 1);
        assert_eq!(t1.tgt_dim(), 1);
        assert_eq!(t2.tgt_dim(), 1);

        assert_eq!(t0.forward_mat() * t0.backward_mat(), SpMat::id(1));
        assert_eq!(t1.forward_mat() * t1.backward_mat(), SpMat::id(1));
        assert_eq!(t2.forward_mat() * t2.backward_mat(), SpMat::id(1));

        let v = SpVec::unit(1, 0);
        let w = t2.backward(&v);
        let x = c[2].devectorize(&w);
        let dx = c.d(2, &x);
        let dw = c[1].vectorize(&dx);
        let dv = t1.forward(&dw);

        assert_eq!(dv.to_dense()[0].abs(), 2);

        let v = SpVec::unit(1, 0);
        let w = t1.backward(&v);
        let z = c[1].devectorize(&w);
        
        assert!(c.d(1, &z).is_zero());
    }
}