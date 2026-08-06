//! [`HomologyCalc`]: computes the homology at one grading index from two
//! consecutive differential matrices via Smith normal form.

use std::marker::PhantomData;
use log::*;

use yui_core::abst::{EucRing, EucRingOps};
use yui_matrix::MatTrait;
use yui_matrix::sparse::*;
use yui_matrix::sparse::snf::{sp_snf, SpSnf};

/// `(rank, torsion_coefficients, optional_basis_change)` returned by
/// [`HomologyCalc::calculate`].
pub type HomologyCalcResult<R> = (usize, Vec<R>, Option<Trans<R>>);

/// Stateless namespace for SNF-based homology computation.
pub struct HomologyCalc<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    _r: PhantomData<R>
}

impl<R> HomologyCalc<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    /// Computes the homology of `d1 → d2` (i.e. `H = ker(d2) / im(d1)`) via two
    /// successive SNF reductions.
    ///
    /// ```text
    ///            d1             d2
    ///    C1 ──────────→ C2 ───────────→ C3
    ///     │              │               │
    ///     │           p1 │               │
    ///     ↓      d1'     ↓               │
    ///    C11 ─────────→ C21              │
    ///     ⊕              ⊕      d2'      │
    ///    C11'           C21'──────────→ C3
    ///                    │               │
    ///                q2⁻¹│               │
    ///                    ↓      d2''     ↓
    ///                   C22 ──────────→ C31
    ///                    ⊕               ⊕
    ///                   C22'            C31'
    ///
    ///  H2 = Ker(d2) / Im(d1)
    ///     ≅ C22' (free) ⊕ (C21 / Im(d1')) (tor)
    /// ```
    pub fn calculate(d1: SpMat<R>, d2: SpMat<R>, with_trans: bool) -> HomologyCalcResult<R> {
        assert_eq!(d1.n_rows(), d2.n_cols());

        if d1.is_zero() && d2.is_zero() {
            return Self::trivial_result(d1.n_rows(), with_trans);
        }

        debug!("calculate homology: {} -> {} -> {}", d1.n_cols(), d1.n_rows(), d2.n_rows());
        Self::log_sparsity("d1", &d1);
        Self::log_sparsity("d2", &d2);

        let (s1, s2) = Self::process_snf(d1, d2, with_trans);
        let (rank, tors) = Self::result(&s1, &s2);

        let trans = if with_trans {
            Some( Self::trans(&s1, &s2) )
        } else {
            None
        };

        (rank, tors, trans)
    }

    fn trivial_result(rank: usize, with_trans: bool) -> HomologyCalcResult<R> {
        let t = with_trans.then(|| Trans::id(rank));
        (rank, vec![], t)
    }

    // sparsity structure of the reducer residual (core = non-zero rows × cols).
    fn log_sparsity(name: &str, d: &SpMat<R>) {
        use std::collections::HashSet;
        let (m, n) = d.shape();
        let nnz = d.iter_nz().count();
        let nz_rows = d.iter_nz().map(|(i, _, _)| i).collect::<HashSet<_>>().len();
        let nz_cols = d.iter_nz().map(|(_, j, _)| j).collect::<HashSet<_>>().len();
        let core_density = if nz_rows * nz_cols > 0 {
            nnz as f64 / (nz_rows * nz_cols) as f64
        } else {
            0.0
        };
        debug!("  {name}: shape ({m}, {n}), nnz {nnz}, nz-rows {nz_rows}, nz-cols {nz_cols}, core {nz_rows}x{nz_cols} (density {core_density:.4})");
    }

    fn process_snf(d1: SpMat<R>, d2: SpMat<R>, with_trans: bool) -> (SpSnf<R>, SpSnf<R>) {
        let n = d1.n_rows();

        let s1 = sp_snf(&d1, [with_trans, true, false, false]);
        let r1 = s1.rank();

        let d2 = if r1 > 0 {
            let p1_inv = s1.pinv().unwrap();
            let t2 = p1_inv.submat_cols(r1..n);
            d2 * &t2 // d2': C21' -> C3
        } else {
            d2
        };

        let s2 = sp_snf(&d2, [false, false, with_trans, with_trans]);

        (s1, s2)
    }

    fn result(s1: &SpSnf<R>, s2: &SpSnf<R>) -> (usize, Vec<R>) {
        let n = s1.result().n_rows();
        let (r1, r2) = (s1.rank(), s2.rank());

        assert!(n >= r1 + r2);

        let rank = n - r1 - r2;

        let tors = s1.factors().into_iter().filter_map(|a| {
            if !a.is_unit() {
                Some(a.clone())
            } else {
                None
            }
        }).collect();

        (rank, tors)
    }

    fn trans(s1: &SpSnf<R>, s2: &SpSnf<R>) -> Trans<R> {
        let n = s1.result().n_rows();
        let (r1, r2) = (s1.rank(), s2.rank());
        let r = n - r1 - r2;
        let t = s1.factors().iter().filter(|a| !a.is_unit()).count();

        let p1 = s1.p().unwrap();                 // size = (n, n)
        let p11 = p1.submat_rows(r1..n);          // size = (n - r1, n)

        let p2 = s2.qinv().unwrap();              // size = (n - r1, n - r1)
        let p22 = p2.submat_rows(r2..n-r1);       // size = (n - (r1 + r2), n - r1)

        let p_free = p22 * p11;                   // size = (n - (r1 + r2), n)
        let p_tor = p1.submat_rows(r1-t..r1);     // size = (t, n)

        let p = SpMat::v_stack(p_free, p_tor);      // size = (r + t, n)

        assert_eq!(p.shape(), (r + t, n));

        let q1 = s1.pinv().unwrap();              // size = (n, n)
        let q12 = q1.submat_cols(r1..n);          // size = (n, n - r1)

        let q2 = s2.q().unwrap();                 // size = (n - r1, n - r1)
        let q22 = q2.submat_cols(r2..n-r1);       // size = (n - r1, n - (r1 + r2))


        let q_free = q12 * q22;                   // size = (n, n - (r1 + r2))
        let q_tor = q1.submat_cols(r1-t..r1);     // size = (n, t)

        let q = SpMat::h_stack(q_free, q_tor);     // size = (n, r + t)

        assert_eq!(q.shape(), (n, r + t));

        Trans::new(p, q)
    }
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use crate::GenericChainComplex1;
    use super::*;

    #[test]
    fn s2_0th() {
        let c = GenericChainComplex1::<i32>::s2();
        let d1 = c.d_matrix(1);
        let d0 = c.d_matrix(0); // zero

        let (rank, tors, t) = HomologyCalc::calculate(d1, d0, true);

        assert_eq!(rank, 1);
        assert_eq!(tors.len(), 0);

        let t = t.unwrap();
        let v = t.backward_mat().col_vec(0);
        let z = c[0].devectorize(&v);

        assert!(!v.is_zero());
        assert!(c.d(0, &z).is_zero());
        assert_eq!(t.forward(&v), SpVec::unit(1, 0));
    }

    #[test]
    fn s2_1st() {
        let c = GenericChainComplex1::<i32>::s2();
        let d1 = c.d_matrix(2);
        let d0 = c.d_matrix(1);

        let (rank, tors, _) = HomologyCalc::calculate(d1, d0, true);

        assert_eq!(rank, 0);
        assert_eq!(tors, vec![]);
    }

    #[test]
    fn s2_2nd() {
        let c = GenericChainComplex1::<i32>::s2();
        let d3 = c.d_matrix(3); // zero
        let d2 = c.d_matrix(2);

        let (rank, tors, t) = HomologyCalc::calculate(d3, d2, true);

        assert_eq!(rank, 1);
        assert_eq!(tors.len(), 0);

        let t = t.unwrap();
        let v = t.backward_mat().col_vec(0);
        let z = c[2].devectorize(&v);

        assert!(!v.is_zero());
        assert!(c.d(2, &z).is_zero());
        assert_eq!(t.forward(&v), SpVec::unit(1, 0));
    }

    #[test]
    fn t2_0th() {
        let c = GenericChainComplex1::<i32>::t2();
        let d1 = c.d_matrix(1);
        let d0 = c.d_matrix(0); // zero

        let (rank, tors, t) = HomologyCalc::calculate(d1, d0, true);

        assert_eq!(rank, 1);
        assert_eq!(tors.len(), 0);

        let t = t.unwrap();
        let v = t.backward_mat().col_vec(0);
        let z = c[0].devectorize(&v);

        assert!(!v.is_zero());
        assert!(c.d(0, &z).is_zero());
        assert_eq!(t.forward(&v), SpVec::unit(1, 0));
    }

    #[test]
    fn t2_1st() {
        let c = GenericChainComplex1::<i32>::t2();
        let d2 = c.d_matrix(2);
        let d1 = c.d_matrix(1);

        let (rank, tors, t) = HomologyCalc::calculate(d2, d1, true);

        assert_eq!(rank, 2);
        assert_eq!(tors.len(), 0);

        let t = t.unwrap();

        for i in 0..2 {
            let v = t.backward_mat().col_vec(i);
            let z = c[1].devectorize(&v);

            assert!(!v.is_zero());
            assert!(c.d(1, &z).is_zero());
            assert_eq!(t.forward(&v), SpVec::unit(2, i));
        }
    }

    #[test]
    fn t2_2nd() {
        let c = GenericChainComplex1::<i32>::t2();
        let d3 = c.d_matrix(3); // zero
        let d2 = c.d_matrix(2);

        let (rank, tors, t) = HomologyCalc::calculate(d3, d2, true);

        assert_eq!(rank, 1);
        assert_eq!(tors.len(), 0);

        let t = t.unwrap();
        let v = t.backward_mat().col_vec(0);
        let z = c[2].devectorize(&v);

        assert!(!v.is_zero());
        assert!(c.d(2, &z).is_zero());
        assert_eq!(t.forward(&v), SpVec::unit(1, 0));
    }

    #[test]
    fn rp2_0th() {
        let c = GenericChainComplex1::<i32>::rp2();
        let d1 = c.d_matrix(1);
        let d0 = c.d_matrix(0); // zero

        let (rank, tors, t) = HomologyCalc::calculate(d1, d0, true);

        assert_eq!(rank, 1);
        assert_eq!(tors.len(), 0);

        let t = t.unwrap();
        let v = t.backward_mat().col_vec(0);
        let z = c[0].devectorize(&v);

        assert!(!v.is_zero());
        assert!(c.d(0, &z).is_zero());
        assert_eq!(t.forward(&v), SpVec::unit(1, 0));
    }

    #[test]
    fn rp2_1st() {
        let c = GenericChainComplex1::<i32>::rp2();
        let d2 = c.d_matrix(2);
        let d1 = c.d_matrix(1);

        let (rank, tors, t) = HomologyCalc::calculate(d2, d1, true);

        assert_eq!(rank, 0);
        assert_eq!(tors, vec![2]);

        let t = t.unwrap();
        let v = t.backward_mat().col_vec(0);
        let z = c[1].devectorize(&v);

        assert!(!v.is_zero());
        assert!(c.d(1, &z).is_zero());
        assert_eq!(t.forward(&v), SpVec::unit(1, 0));
    }
}