use log::info;
use yui::{Ring, RingOps};
use super::*;
use super::triang::{TriangularType, solve_triangular, solve_triangular_left};

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

const LOG_THRESHOLD: usize = 10_000;

pub struct Schur<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    orig_shape: (usize, usize),
    compl: SpMat<R>,
    t_in: Option<SpMat<R>>,  // [-a⁻¹b; 1]
    t_out: Option<SpMat<R>>, // [-ca⁻¹  1]
}

impl<R> Schur<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_partial_lower(a: &SpMat<R>, r: usize, with_trans: bool) -> Self {
        assert!(r <= a.nrows());
        assert!(r <= a.ncols());

        let orig_shape = a.shape();
        let (compl, t_in, t_out) = Self::compute_from_partial_lower(a, r, with_trans);

        Self { 
            orig_shape, 
            compl, 
            t_in,
            t_out
        }
    }

    fn compute_from_partial_lower(a: &SpMat<R>, r: usize, with_trans: bool) -> (SpMat<R>, Option<SpMat<R>>, Option<SpMat<R>>) {
        let report = should_report(a);
        if report { 
            info!("schur, a: {:?}, r: {} ..", a.shape(), r);
        }

        let (m, n) = a.shape();

        //  [a   ][a⁻¹b   ] = [b]
        //  [c  1][d-ca⁻¹b]   [d]
        //  ~~~~~~
        //    = p
        
        let p = SpMat::from_entries((m, m), Iterator::chain(
            a.view().submat_cols(0..r).iter().map(|(i, j, x)| (i, j, x.clone())),
            (r..m).map(|i| (i, i, R::one()))
        ));
        
        let bd = a.submat_cols(r..n);
        let pinvbd = solve_triangular(TriangularType::Lower, &p, &bd);
        let s = pinvbd.submat_rows(r..m);

        if report { 
            info!("schur: {:?}", s.shape());
        }

        let (t_in, t_out) = if with_trans { 
            let id = |n| SpMat::id(n);

            // t_in = [-a⁻¹b]
            //        [  1  ]
            
            let ainvb = pinvbd.submat_rows(0..r);
            let t_in = (-ainvb).stack(&id(n - r));

            // [-ca⁻¹  1][a   ] = [0  1]
            //           [c  1] 
            // ~~~~~~~~~~
            //   = t_out

            let i = SpMat::zero((m - r, r)).concat(&SpMat::id(m - r));
            let t_out = solve_triangular_left(TriangularType::Lower, &p, &i);

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

fn should_report<R>(a: &SpMat<R>) -> bool { 
    usize::min(a.nrows(), a.ncols()) > LOG_THRESHOLD && log::max_level() >= log::LevelFilter::Info
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn schur() {
        let a = SpMat::from_dense_data((6, 5), [
            1, 0, 0, 1, 3,
            2,-1, 0, 2, 2,
            3, 2, 1, 0, 3,
            4, 2, 4,-3, 0,
            5, 3, 5, 2, 2,
            6, 2,-3, 1, 8
        ]);
        let sch = Schur::from_partial_lower(&a, 3, false);
        let s = sch.complement();

        assert_eq!(s, &SpMat::from_dense_data((3,2), [5,36,12,45,-14,-60]));
        assert!(sch.trans_in().is_none());
        assert!(sch.trans_out().is_none());
    }

    #[test]
    fn schur_with_trans() {
        let a = SpMat::from_dense_data((6, 5), [
            1, 0, 0, 1, 3,
            2,-1, 0, 2, 2,
            3, 2, 1, 0, 3,
            4, 2, 4,-3, 0,
            5, 3, 5, 2, 2,
            6, 2,-3, 1, 8
        ]);
        let sch = Schur::from_partial_lower(&a, 3, true);
        let s = sch.complement();

        assert_eq!(s, &SpMat::from_dense_data((3,2), [5,36,12,45,-14,-60]));
        assert!(sch.trans_in().is_some());
        assert!(sch.trans_out().is_some());

        let t_in = sch.trans_in().unwrap();
        let t_out = sch.trans_out().unwrap();

        assert_eq!(t_in,  &SpMat::from_dense_data((5,2), [-1,-3,0,-4,3,14,1,0,0,1]));
        assert_eq!(t_out, &SpMat::from_dense_data((3,6), [20,-6,-4,1,0,0,24,-7,-5,0,1,0,-31,8,3,0,0,1]));
        assert_eq!(&(t_out * &a * t_in), s);
    }
}