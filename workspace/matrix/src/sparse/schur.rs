use yui_core::{Ring, RingOps};
use super::*;
use super::triang::{TriangularType, solve_triangular};

//                [a  b]
//                [c  d]
//            X ----------> Y
//  [1 -a⁻¹b] ^             | [1      ]
//  [     1 ] |             | [-ca⁻¹ 1]
//            |             V
//            X ----------> Y
//                [a   ] 
//                [   s]
//
// s = d - c a⁻¹ b

pub struct Schur<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    orig_shape: (usize, usize),
    compl: SpMat<R>,
    t_in: Option<SpMat<R>>,  // -a⁻¹b
    t_out: Option<SpMat<R>>, // -ca⁻¹
}

impl<'a, R> Schur<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_partial_lower(a: &SpMat<R>, r: usize, with_trans: bool) -> Self {
        assert!(r <= a.rows());
        assert!(r <= a.cols());

        let orig_shape = a.shape();
        let (compl, t_in, t_out) = Self::compute_from_partial_lower(&a, r, with_trans);

        Self { 
            orig_shape, 
            compl, 
            t_in,
            t_out
        }
    }

    // Solving
    //
    //   [a    ][x1] = [b]
    //   [c  id][x2]   [d]
    //
    // gives 
    //
    //   x1 = a⁻¹ b, 
    //   x2 = d - c a⁻¹ b
    fn compute_from_partial_lower(a: &SpMat<R>, r: usize, with_trans: bool) -> (SpMat<R>, Option<SpMat<R>>, Option<SpMat<R>>) {
        let (m, n) = a.shape();
        let p = SpMat::generate((m, m), |set| { 
            // [a; c]
            for (i, j, x) in a.view().submat_cols(0..r).iter() {
                set(i, j, x.clone());
            }
            // [0; 1]
            for i in r..m { 
                set(i, i, R::one());
            }
        });
        let bd = a.submat_cols(r..n);
        
        let x = solve_triangular(TriangularType::Lower, &p, &bd);
        let s = x.submat_rows(r..m);

        let (t_in, t_out) = if with_trans { 
            let t_in = -x.submat_rows(0..r);
    
            let i = SpMat::generate((m, r), |set| { 
                // [0; 1]
                for i in 0..m-r { 
                    set(i, i, R::one());
                }
            });
    
            // -c a⁻¹
            let t_out = solve_triangular(TriangularType::Lower, &p, &i).submat_rows(r..m);

            (Some(t_in), Some(t_out))
        } else { 
            (None, None)
        };

        (s, t_in, t_out)
    }

    pub fn orig_shape(&self) -> (usize, usize) { 
        self.orig_shape
    }

    pub fn complement(&self) -> &SpMat<R> {
        &self.compl
    }

    pub fn complement_into(self) -> SpMat<R> {
        self.compl
    }

    // v -> [-a⁻¹b v; v]
    pub fn trans_in(&self, v: SpVec<R>) -> SpVec<R> { 
        let r = self.compl.shape().1;

        assert_eq!(v.dim(), r);

        let Some(t_in) = &self.t_in else { 
            panic!()
        };
        let v1 = t_in * &v;

        v1.stack(&v)
    }

    // [v1; v2] -> v2 - c a⁻¹ v1
    pub fn trans_out(&self, v: SpVec<R>) -> SpVec<R> { 
        let m = self.orig_shape.0;

        assert_eq!(v.dim(), m);

        let Some(t_out) = &self.t_out else { 
            panic!()
        };
        let r = self.compl.shape().0;
        let v1 = v.subvec(0..r);
        let v2 = v.subvec(r..m);

        v2 + t_out * v1
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn schur_compl() {
        let a = SpMat::from_vec((6, 5), vec![
            1, 0, 0, 1, 3,
            2,-1, 0, 2, 2,
            3, 2, 1, 0, 3,
            4, 2, 4,-3, 0,
            5, 3, 5, 2, 2,
            6, 2,-3, 1, 8
        ]);
        let s = Schur::from_partial_lower(&a, 3, false);

        assert_eq!(s.complement(), &SpMat::from_vec((3,2), vec![5,36,12,45,-14,-60]));
        assert_eq!(s.t_in,  None);
        assert_eq!(s.t_out, None);
    }

    #[test]
    fn trans_vec() {
        let a = SpMat::from_vec((6, 5), vec![
            1, 0, 0, 1, 3,
            2,-1, 0, 2, 2,
            3, 2, 1, 0, 3,
            4, 2, 4,-3, 0,
            5, 3, 5, 2, 2,
            6, 2,-3, 1, 8
        ]);
        let s = Schur::from_partial_lower(&a, 3, true);

        assert_eq!(s.complement(), &SpMat::from_vec((3,2), vec![5,36,12,45,-14,-60]));
        assert_eq!(s.t_in,  Some(SpMat::from_vec((3,2), vec![-1,-3,0,-4,3,14])));
        assert_eq!(s.t_out, Some(SpMat::from_vec((3,3), vec![20,-6,-4,24,-7,-5,-31,8,3])));

        let v = SpVec::from(vec![1,2]);
        let w = s.trans_in(v);
        assert_eq!(w, SpVec::from(vec![-7,-8,31,1,2]));

        let v = SpVec::from(vec![1,2,0,-3,2,1]);
        let w = s.trans_out(v);
        assert_eq!(w, SpVec::from(vec![5,12,-14]));
    }
}