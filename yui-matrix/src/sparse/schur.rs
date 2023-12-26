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
    pub fn from_partial_triangular(t: TriangularType, abcd: &SpMat<R>, r: usize, with_trans: bool) -> Self {
        assert!(r <= abcd.nrows());
        assert!(r <= abcd.ncols());

        trace!("schur, a: {:?}, r: {} ..", abcd.shape(), r);

        let id = |n| SpMat::<R>::id(n);
        
        let (m, n) = abcd.shape();
        let [a, b, c, d] = abcd.divide4((r, r));

        let ainvb = solve_triangular(t, &a, &b); // ax = b
        let s = Self::compute_schur(&ainvb, &c, &d);

        trace!("schur: {:?}", s.shape());

        // t_in = [-a⁻¹b]
        //        [  1  ]
        
        let t_in = if with_trans { 
            let t_in = (-ainvb).stack(&id(n - r));
            Some(t_in)
        } else { 
            None
        };

        // t_out = [-ca⁻¹  1]
        let t_out = if with_trans { 
            let mut t_out = -solve_triangular_left(t, &a, &c); // (-x)a = c
            t_out.extend_cols(id(m - r));
            Some(t_out)
        } else { 
            None
        };

        Self { 
            orig_shape: (m, n), 
            compl: s, 
            t_in,
            t_out
        }
    }

    fn compute_schur(ainvb: &SpMat<R>, c: &SpMat<R>, d: &SpMat<R>) -> SpMat<R> { 
        d - c * ainvb
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

    pub fn transpose(&self) -> Self { 
        Self {
            orig_shape: (self.orig_shape.1, self.orig_shape.0),
            compl: self.compl.transpose(),
            t_in: self.t_out.as_ref().map(|t| t.transpose()),
            t_out: self.t_in.as_ref().map(|t| t.transpose()),
        }
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn schur_lower() {
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

        assert_eq!(s, &SpMat::from_dense_data((3,2), [
             5,  36, 
             12, 45,
            -14,-60
        ]));
        assert!(sch.trans_in().is_none());
        assert!(sch.trans_out().is_none());
    }

    #[test]
    fn schur_lower_with_trans() {
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

        assert_eq!(s, &SpMat::from_dense_data((3,2), [
             5,  36, 
             12, 45,
            -14,-60
        ]));
        assert!(sch.trans_in().is_some());
        assert!(sch.trans_out().is_some());

        let t_in = sch.trans_in().unwrap();
        let t_out = sch.trans_out().unwrap();

        assert_eq!(t_in,  &SpMat::from_dense_data((5,2), [
            -1, -3,
             0, -4,
             3, 14,
             1,  0,
             0,  1
        ]));
        assert_eq!(t_out, &SpMat::from_dense_data((3,6), [
             20, -6, -4, 1, 0, 0,
             24, -7, -5, 0, 1, 0,
            -31,  8,  3, 0, 0, 1
        ]));
        assert_eq!(&(t_out * &a * t_in), s);
    }

    #[test]
    fn schur_upper() {
        let a = SpMat::from_dense_data((5, 6), [
            1, 2, 3, 4, 5, 6,
            0, -1, 2, 2, 3, 2,
            0, 0, 1, 4, 5, -3,
            1, 2, 0, -3, 2, 1,
            3, 2, 3, 0, 2, 8,
        ]);
        let sch = Schur::from_partial_triangular(TriangularType::Upper, &a, 3, false);
        let s = sch.complement();

        assert_eq!(s, &SpMat::from_dense_data((2, 3), [
            5, 12,-14,
            36,45,-60
        ]));
        assert!(sch.trans_in().is_none());
        assert!(sch.trans_out().is_none());
    }

    #[test]
    fn schur_upper_with_trans() {
        let a = SpMat::from_dense_data((5, 6), [
            1, 2, 3, 4, 5, 6,
            0, -1, 2, 2, 3, 2,
            0, 0, 1, 4, 5, -3,
            1, 2, 0, -3, 2, 1,
            3, 2, 3, 0, 2, 8,
        ]);
        let sch = Schur::from_partial_triangular(TriangularType::Upper, &a, 3, true);
        let s = sch.complement();

        assert_eq!(s, &SpMat::from_dense_data((2, 3), [
            5, 12,-14,
            36,45,-60
        ]));
        assert!(sch.trans_in().is_some());
        assert!(sch.trans_out().is_some());

        let t_in = sch.trans_in().unwrap();
        let t_out = sch.trans_out().unwrap();

        assert_eq!(t_in,  &SpMat::from_dense_data((6,3), [
            20, 24, -31,
            -6, -7,   8,
            -4, -5,   3,
             1,  0,   0,
             0,  1,   0,
             0,  0,   1
        ]));
        assert_eq!(t_out, &SpMat::from_dense_data((2, 5), [
            -1,  0,  3, 1, 0,
            -3, -4, 14, 0, 1
        ]));
        assert_eq!(&(t_out * &a * t_in), s);
    }
}