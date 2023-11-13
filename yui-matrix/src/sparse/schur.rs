use yui::{Ring, RingOps};
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
    t_in: Option<SpMat<R>>,  // [-a⁻¹b; 1]
    t_out: Option<SpMat<R>>, // [-ca⁻¹  1]
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

    fn compute_from_partial_lower(a: &SpMat<R>, r: usize, with_trans: bool) -> (SpMat<R>, Option<SpMat<R>>, Option<SpMat<R>>) {
        let (m, n) = a.shape();

        //  [a   ][x] = [b]  ==>  [x] = [a⁻¹b     ]
        //  [c  1][y]   [d]       [y]   [d - ca⁻¹b]
        //  ~~~~~~
        //    p
        
        let p = SpMat::from_entries((m, m), Iterator::chain(
            a.view().submat_cols(0..r).iter().map(|(i, j, x)| (i, j, x.clone())),
            (r..m).map(|i| (i, i, R::one()))
        ));
        let bd = a.submat_cols(r..n);
        
        let x = solve_triangular(TriangularType::Lower, &p, &bd);
        let s = x.submat_rows(r..m);

        let (t_in, t_out) = if with_trans { 
            let id = |n| SpMat::id(n);

            // t_in = [-a⁻¹b]
            //        [  1  ]
            let n_ainvb = -x.submat_rows(0..r);
            let t_in = n_ainvb.stack(&id(n - r));
    
            // [a   ][x] = [1]  ==>  [x] = [a⁻¹  ] 
            // [c  1][y]   [0]       [y]   [-ca⁻¹]
            let i = SpMat::id(r).stack(&SpMat::zero((m - r, r)));
            let y = solve_triangular(TriangularType::Lower, &p, &i).submat_rows(r..m);

            // t_out = [-ca⁻¹  1]
            let t_out = y.concat(&id(m - r));

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

    pub fn trans_in(&self) -> Option<&SpMat<R>> { 
        self.t_in.as_ref()
    }

    pub fn trans_out(&self) -> Option<&SpMat<R>> { 
        self.t_out.as_ref()
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
        assert_eq!(s.trans_in(),  Some(&SpMat::from_vec((5,2), vec![-1,-3,0,-4,3,14,1,0,0,1])));
        assert_eq!(s.trans_out(), Some(&SpMat::from_vec((3,6), vec![20,-6,-4,1,0,0,24,-7,-5,0,1,0,-31,8,3,0,0,1])));

        let v = SpVec::from(vec![1,2]);
        let w = s.trans_in().unwrap() * v;
        assert_eq!(w, SpVec::from(vec![-7,-8,31,1,2]));

        let v = SpVec::from(vec![1,2,0,-3,2,1]);
        let w = s.trans_out().unwrap() * v;
        assert_eq!(w, SpVec::from(vec![5,12,-14]));
    }
}