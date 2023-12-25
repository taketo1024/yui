use log::trace;
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

pub struct Schur<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    orig_shape: (usize, usize),
    compl: SpMat<R>,
    t_in: Option<SpMat<R>>,  // [-a⁻¹b; 1]
    t_out: Option<SpMat<R>>, // [-ca⁻¹  1]
}

impl<R> Schur<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn from_partial_triangular(t: TriangularType, a: &SpMat<R>, r: usize, with_trans: bool) -> Self {
        assert!(r <= a.nrows());
        assert!(r <= a.ncols());

        let orig_shape = a.shape();
        let (compl, t_in, t_out) = match t {
            TriangularType::Upper => todo!(),
            TriangularType::Lower => Self::compute_from_partial_lower(a, r, with_trans)
        };

        Self { 
            orig_shape, 
            compl, 
            t_in,
            t_out
        }
    }

    fn compute_from_partial_lower(abcd: &SpMat<R>, r: usize, with_trans: bool) -> (SpMat<R>, Option<SpMat<R>>, Option<SpMat<R>>) {
        use TriangularType::Lower as L;
        let id = |n| SpMat::<R>::id(n);
        let zero = |m, n| SpMat::<R>::zero((m, n));

        trace!("schur, a: {:?}, r: {} ..", abcd.shape(), r);

        //  [a   ] [a⁻¹b   ] = [b]
        //  [c  1] [d-ca⁻¹b]   [d]
        //  ~~~~~~ ~~~~~~~~~   ~~~
        //   = p    = pinvbd   = bd

        let (m, n) = abcd.shape();

        let pinvbd = { 
            let p = abcd.submat_cols(0..r).concat(
                &zero(r, m - r).stack(&id(m - r))
            );
            let bd = abcd.submat_cols(r..n);
            solve_triangular(L, &p, &bd) // px = bd
        };
        
        let s = pinvbd.submat_rows(r..m);

        trace!("schur: {:?}", s.shape());

        // t_in = [-a⁻¹b]
        //        [  1  ]
        let t_in = if with_trans { 
            let ainvb = pinvbd.submat_rows(0..r);
            let t_in = (-ainvb).stack(&id(n - r));
            Some(t_in)
        } else { 
            None
        };

        std::mem::drop(pinvbd);

        // t_out = [-ca⁻¹  1]
        let t_out = if with_trans { 
            let a = abcd.submat(0..r, 0..r);
            let c = abcd.submat(r..m, 0..r);
            let cainv = solve_triangular_left(L, &a, &c); // xa = c
            let t_out = (-cainv).concat(&id(m - r));

            Some(t_out)
        } else { 
            None
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
    fn schur() {
        let a = SpMat::from_dense_data((6, 5), [
            1, 0, 0, 1, 3,
            2,-1, 0, 2, 2,
            3, 2, 1, 0, 3,
            4, 2, 4,-3, 0,
            5, 3, 5, 2, 2,
            6, 2,-3, 1, 8
        ]);
        let sch = Schur::from_partial_triangular(TriangularType::Lower, &a, 3, false);
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
        let sch = Schur::from_partial_triangular(TriangularType::Lower, &a, 3, true);
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