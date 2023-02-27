use core::panic;
use std::cmp::min;
use log::{info, trace};
use crate::math::traits::{EucRing, EucRingOps};
use super::dense::*;
use super::lll::{LLLRing, LLLRingOps, lll_hnf_in_place};

pub type SnfFlags = [bool; 4];

pub fn snf<R>(target: &DnsMat<R>, flags: SnfFlags) -> SnfResult<R>
where R: EucRing, for<'a> &'a R: EucRingOps<R> {
    let copy = target.clone();
    snf_in_place(copy, flags)
}

pub fn snf_in_place<R>(target: DnsMat<R>, flags: SnfFlags) -> SnfResult<R>
where R: EucRing, for<'a> &'a R: EucRingOps<R> {
    info!("start snf: {:?}, flags: {:?}.\n{}", target.shape(), flags, target);

    let mut calc = SnfCalc::new(target, flags);

    calc.process();

    info!("snf done.\n{}", calc.target);

    calc.result()
}

#[derive(Debug)]
pub struct SnfResult<R>
where R: EucRing, for<'a> &'a R: EucRingOps<R> { 
    result: DnsMat<R>,
    p:    Option<DnsMat<R>>,
    pinv: Option<DnsMat<R>>,
    q:    Option<DnsMat<R>>,
    qinv: Option<DnsMat<R>>
}

impl<R> SnfResult<R>
where R: EucRing, for<'a> &'a R: EucRingOps<R> { 
    pub fn result(&self) -> &DnsMat<R> { 
        &self.result
    }

    pub fn p(&self) -> Option<&DnsMat<R>> {
        self.p.as_ref()
    }

    pub fn pinv(&self) -> Option<&DnsMat<R>> {
        self.pinv.as_ref()
    }

    pub fn q(&self) -> Option<&DnsMat<R>> {
        self.q.as_ref()
    }

    pub fn qinv(&self) -> Option<&DnsMat<R>> {
        self.qinv.as_ref()
    }

    pub fn trans(&self) -> [Option<&DnsMat<R>>; 4] {
        [self.p.as_ref(),
         self.pinv.as_ref(),
         self.q.as_ref(),
         self.qinv.as_ref()]
    }

    pub fn destruct(self) -> (DnsMat<R>, [Option<DnsMat<R>>; 4]) {
        (self.result, [self.p, self.pinv, self.q, self.qinv])
    }

    pub fn rank(&self) -> usize {
        let n = min(self.result.rows(), self.result.cols());
        for i in 0..n { 
            if self.result[[i, i]].is_zero() { 
                return i
            }
        }
        return n
    }

    pub fn factors(&self) -> Vec<&R> { 
        let n = min(self.result.rows(), self.result.cols());
        (0..n).filter_map(|i| { 
            let a = &self.result[[i, i]];
            if !a.is_zero() { 
                Some(a)
            } else {
                None
            }
         }).collect()
    }
}

#[derive(Debug)]
pub struct SnfCalc<R>
where R: EucRing, for<'a> &'a R: EucRingOps<R> {
    target: DnsMat<R>,
    p:    Option<DnsMat<R>>,
    pinv: Option<DnsMat<R>>,
    q:    Option<DnsMat<R>>,
    qinv: Option<DnsMat<R>>
}

impl<R> SnfCalc<R>
where R: EucRing, for<'a> &'a R: EucRingOps<R> {
    pub fn new(target: DnsMat<R>, flags: SnfFlags) -> Self { 
        let id_opt = |size, flag| {
            if flag{ Some(DnsMat::id(size)) } else { None }
        };

        let (m, n) = target.shape();
        let p    = id_opt(m, flags[0]);
        let pinv = id_opt(m, flags[1]);
        let q    = id_opt(n, flags[2]);
        let qinv = id_opt(n, flags[3]);

        SnfCalc{ target, p, pinv, q, qinv }
    }

    pub fn result(self) -> SnfResult<R> {
        SnfResult { 
            result: self.target, 
            p: self.p,
            pinv: self.pinv,
            q: self.q,
            qinv: self.qinv
        }
    }

    pub fn process(&mut self) { 
        if self.target.is_zero() { 
            return
        }
        
        self.preprocess();
        self.eliminate_all();
        self.diag_normalize();
    }

    fn preprocess(&mut self) {
        use num_bigint::BigInt;
        use crate::math::types::quad_int::{GaussInt, EisenInt};
        preprocess_lll_for!(self, 
            i64, i128, BigInt, 
            GaussInt<i64>, GaussInt<i128>, GaussInt<BigInt>, 
            EisenInt<i64>, EisenInt<i128>, EisenInt<BigInt>
        );
    }

    fn eliminate_all(&mut self) {
        let (m, n) = self.target.shape();
        let mut i = 0;

        for j in 0..n { 
            if i >= m { break }
            if self.eliminate_step(i, j) { 
                i += 1;
            }
        }
    }

    fn eliminate_step(&mut self, i: usize, j: usize) -> bool {
        // select pivot
        let Some(i_p) = self.select_pivot(i, j) else { 
            return false 
        };

        trace!("select-pivot: ({i_p}, {j})");

        // swap rows
        if i_p > i { 
            self.swap_rows(i, i_p);
        }

        // swap cols
        if j > i { 
            self.swap_cols(i, j);
        }

        // normalize pivot
        let u = self.target[[i, i]].normalizing_unit();
        if !u.is_one() { 
            self.mul_col(i, &u);
        }

        // eliminate row and col
        self.eliminate_at(i, i);

        return true
    }

    fn row_nz(&self, i: usize) -> usize { 
        self.target.array().row(i).iter().filter(|a| !a.is_zero()).count()
    }

    fn col_nz(&self, j: usize) -> usize { 
        self.target.array().column(j).iter().filter(|a| !a.is_zero()).count()
    }

    fn swap_rows(&mut self, i: usize, j: usize) {
        self.target.swap_rows(i, j);
        self.p.as_mut().map( |p| p.swap_rows(i, j) );
        self.pinv.as_mut().map( |pinv| pinv.swap_cols(i, j) );

        trace!("swap-rows: ({i}, {j})\n{}", self.target);
    }

    fn swap_cols(&mut self, i: usize, j: usize) {
        self.target.swap_cols(i, j);
        self.q.as_mut().map( |q| q.swap_cols(i, j) );
        self.qinv.as_mut().map( |qinv| qinv.swap_rows(i, j) );

        trace!("swap-cols: ({i}, {j})\n{}", self.target);
    }

    fn mul_row(&mut self, i: usize, u: &R) {
        self.target.mul_row(i, u);
        self.p.as_mut().map( |p| p.mul_row(i, u));
        self.pinv.as_mut().map( |pinv| {
            let Some(uinv) = &u.inv() else { panic!("`u` is not invertible.") };
            pinv.mul_col(i, uinv) 
        });

        trace!("mul-row: {i} by {u})\n{}", self.target);
    }
    
    fn mul_col(&mut self, i: usize, u: &R) {
        self.target.mul_col(i, u);
        self.q.as_mut().map( |q| q.mul_col(i, u) );
        self.qinv.as_mut().map( |qinv| {
            let Some(uinv) = &u.inv() else { panic!("`u` is not invertible.") };
            qinv.mul_row(i, uinv) 
        });

        trace!("mul-col: {i} by {u})\n{}", self.target);
    }

    // Multiply [a, b; c, d] from left, assuming det = 1.
    pub fn left_elementary(&mut self, comps: [&R; 4], i: usize, j: usize) { 
        let [a, b, c, d] = comps;
        debug_assert!((a * d - b * c).is_one());

        self.target.left_elementary(comps, i, j);
        self.p.as_mut().map( |p| p.left_elementary(comps, i, j) ); 
        self.pinv.as_mut().map(|pinv| { 
            let inv_t = [d, &-c, &-b, a];
            pinv.right_elementary(inv_t, i, j) 
        }); 

        trace!("left-elem: [{a}, {b}; {c}, {d}] for rows ({i}, {j})).\n{}", self.target);
    }

    // Multiply [a, c; b, d] from right, assuming det = 1. 
    pub fn right_elementary(&mut self, comps: [&R; 4], i: usize, j: usize) { 
        let [a, b, c, d] = comps;
        debug_assert!((a * d - b * c).is_one());
        
        self.target.right_elementary(comps, i, j);
        self.q.as_mut().map( |q| q.right_elementary(comps, i, j) ); 
        self.qinv.as_mut().map(|qinv| { 
            let inv_t = [d, &-c, &-b, a];
            qinv.left_elementary(inv_t, i, j) 
        }); 

        trace!("right-elem: [{a}, {b}; {c}, {d}] for cols ({i}, {j})).\n{}", self.target);
    }

    fn select_pivot(&self, below_i: usize, j: usize) -> Option<usize> { 
        // find row `i` below `below_i` with minimum nnz. 
        (below_i..self.target.rows())
            .filter( |i| !self.target[[*i, j]].is_zero() )
            .map( |i| (i, self.row_nz(i)) )
            .min_by( |e1, e2| e1.1.cmp(&e2.1) )
            .map( |(i, _)| i )
    }

    fn eliminate_at(&mut self, i: usize, j: usize) {
        assert!(!self.target[[i, j]].is_zero());

        while self.row_nz(i) > 1 || self.col_nz(j) > 1 { 
            let modified = self.eliminate_col(i, j)
                         | self.eliminate_row(i, j);
            if !modified {
                panic!("Detect endless loop");
            }
        }
    }

    fn eliminate_row(&mut self, i: usize, j: usize) -> bool { 
        let mut modified = false;

        for j1 in 0..self.target.cols() {
            if j == j1 || self.target[[i, j1]].is_zero() { continue }

            // d = sx + ty,
            // a = x/d,
            // b = y/d.
        
            // [x y][s -b] = [d 0]
            //      [t  a]   

            let x = &self.target[[i, j ]];
            let y = &self.target[[i, j1]];

            let (d, s, t) = Self::gcdx(x, y);
            let (a, b) = (x / &d, y / &d);

            self.right_elementary(
                [&s, &t, &-b, &a], 
                j, j1
            );
            modified = true
        }

        modified
    }
    
    fn eliminate_col(&mut self, i: usize, j: usize) -> bool { 
        let mut modified = false;

        for i1 in 0..self.target.rows() {
            if i == i1 || self.target[[i1, j]].is_zero() { continue }

            // d = sx + ty,
            // a = x/d,
            // b = y/d.
        
            // [ s t][x] < i  = [d]
            // [-b a][y] < i1   [0]

            let x = &self.target[[i , j]];
            let y = &self.target[[i1, j]];

            let (d, s, t) = Self::gcdx(x, y);
            let (a, b) = (x / &d, y / &d);

            self.left_elementary(
                [&s, &t, &-b, &a], 
                i, i1
            );
            modified = true
        }
        
        modified
    }
    
    fn diag_normalize(&mut self) {
        debug_assert!(self.target.is_diag());

        let r = min(self.target.rows(), self.target.cols());
        if r == 0 { 
            return
        }

        loop { 
            let mut done = true;
            for i in 0..r-1 { 
                done &= self.diag_normalize_step(i);
            }
            if done { break }
        }

        for i in 0..r { 
            let a = &self.target[[i, i]];
            let u = a.normalizing_unit();
            if !u.is_one() {
                self.mul_row(i, &u);
            }
        }
    }

    fn diag_normalize_step(&mut self, i: usize) -> bool {
        let x = &self.target[[i, i]];
        let y = &self.target[[i + 1, i + 1]];

        if 
             x.is_zero() && y.is_zero() || 
            !x.is_zero() && x.divides(y)
        { 
            return true
        }

        // sx + ty = d, a = x/d, b = y/d.
        //
        // [1   1 ][x   ][s  -b] = [d      ]
        // [-tb sa][   y][t   a]   [   xy/d]

        let (d, s, t) = Self::gcdx(x, y);
        let (a, b) = (x / &d, y / &d);
        let (tb, sa) = (&t * &b, &s * &a);

        self.left_elementary(
            [&R::one(), &R::one(), &-tb, &sa], 
            i, i + 1
        );
        self.right_elementary(
            [&s, &t, &-b, &a], 
            i, i + 1
        );

        false
    }

    fn gcdx(x: &R, y: &R) -> (R, R, R) { 
        let (mut d, mut s, mut t) = EucRing::gcdx(x, y);

        let u = d.normalizing_unit();
        if !u.is_one() {
            (d, s, t) = (&d * &u, &s * &u, &t * &u);
        }

        let a = x / &d;
        if a.is_unit() { 
            (d, a, R::zero())
        } else {
            (d, s, t)
        }
    }
}

impl<R> SnfCalc<R>
where R: LLLRing, for<'a> &'a R: LLLRingOps<R> {
    fn preprocess_lll(&mut self) {
        info!("start lll-preprocess, type = {}", std::any::type_name::<R>());

        let flag = [self.p.is_some(), self.pinv.is_some()];
        
        let b = std::mem::take(&mut self.target);
        let (res, p, pinv) = lll_hnf_in_place(b, flag);

        self.target = res;
        self.p = p;
        self.pinv = pinv;

        info!("preprocess done.\n{}", self.target);
    }
}

macro_rules! preprocess_lll_expand {
    ($any:ident) => {};
    ($any:ident, $t:ty $(,$next:ty)*) => {{
        if let Some(_self) = $any.downcast_mut::<SnfCalc<$t>>() {
            _self.preprocess_lll()
        } else {
            preprocess_lll_expand!($any $(,$next)*);
        }
    }};
}

macro_rules! preprocess_lll_for {
    ($self:ident, $t:ty $(,$next:ty)*) => {{
        let any: &mut dyn std::any::Any = $self;
        preprocess_lll_expand!(any, $t, $($next),*);
    }};
}

use {preprocess_lll_for, preprocess_lll_expand};

#[cfg(test)]
mod tests {
    use ndarray::array;
    use super::*;

    #[test]
    fn init() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6]]);
        let calc = SnfCalc::new(a, [true; 4]);
        let (res, [p, pinv, q, qinv]) = calc.result().destruct();

        assert_eq!(res, DnsMat::from(array![[1,2,3], [4,5,6]]));
        assert_eq!(p,    Some(DnsMat::id(2)));
        assert_eq!(pinv, Some(DnsMat::id(2)));
        assert_eq!(q,    Some(DnsMat::id(3)));
        assert_eq!(qinv, Some(DnsMat::id(3)));
    }

    #[test]
    fn init_no_pq() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6]]);
        let calc = SnfCalc::new(a, [false; 4]);
        
        let (res, [p, pinv, q, qinv]) = calc.result().destruct();

        assert_eq!(res, DnsMat::from(array![[1,2,3], [4,5,6]]));
        assert_eq!(p,    None);
        assert_eq!(pinv, None);
        assert_eq!(q,    None);
        assert_eq!(qinv, None);
    }

    #[test]
    fn row_nz() { 
        let a = DnsMat::from(array![[1,0,0], [0,5,6], [7,8,9]]);
        let calc = SnfCalc::new(a, [false; 4]);
        assert_eq!(calc.row_nz(0), 1);
        assert_eq!(calc.row_nz(1), 2);
        assert_eq!(calc.row_nz(2), 3);
    }

    #[test]
    fn col_nz() { 
        let a = DnsMat::from(array![[1,0,0], [0,5,6], [7,8,9]]);
        let calc = SnfCalc::new(a, [false; 4]);
        assert_eq!(calc.col_nz(0), 2);
        assert_eq!(calc.col_nz(1), 2);
        assert_eq!(calc.col_nz(2), 2);
    }

    #[test]
    fn swap_rows() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.swap_rows(0, 1);

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res,  DnsMat::from(array![[4,5,6], [1,2,3], [7,8,9]]));
        assert_eq!(p * a.clone(), res);
        assert_eq!(pinv * res, a);
        assert!(q.is_id());
        assert!(qinv.is_id());
    }

    #[test]
    fn swap_cols() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.swap_cols(0, 1);

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res,  DnsMat::from(array![[2,1,3], [5,4,6], [8,7,9]]));
        assert!(p   .is_id());
        assert!(pinv.is_id());
        assert_eq!(a.clone() * q, res);
        assert_eq!(res * qinv, a);
    }

    #[test]
    fn mul_row() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.mul_row(0, &-1);
        
        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res,  DnsMat::from(array![[-1,-2,-3], [4,5,6], [7,8,9]]));
        assert_eq!(p * a.clone(), res);
        assert_eq!(pinv * res, a);
        assert!(q.is_id());
        assert!(qinv.is_id());
    }

    #[test]
    fn mul_col() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.mul_col(0, &-1);
        
        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res,  DnsMat::from(array![[-1,2,3], [-4,5,6], [-7,8,9]]));
        assert!(p.is_id());
        assert!(pinv.is_id());
        assert_eq!(a.clone() * q, res);
        assert_eq!(res * qinv, a);
    }

    #[test]
    fn left_elementary() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let e = [&3,&2,&4,&3]; // det = 1
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.left_elementary(e, 0, 1);

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res,  DnsMat::from(array![[11,16,21], [16,23,30], [7,8,9]]));
        assert_eq!(p * a.clone(), res);
        assert_eq!(pinv * res, a);
        assert!(q.is_id());
        assert!(qinv.is_id());
    }

    #[test]
    fn right_elementary() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let e = [&3,&2,&4,&3]; // det = 1
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.right_elementary(e, 0, 1);

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res,  DnsMat::from(array![[7,10,3], [22,31,6], [37,52,9]]));
        assert!(p.is_id());
        assert!(pinv.is_id());
        assert_eq!(a.clone() * q, res);
        assert_eq!(res * qinv, a);
    }

    #[test]
    fn gcdx() {
        let (x, y) = (14, -52);
        let (d, s, t) = SnfCalc::gcdx(&x, &y);
        assert_eq!(d, 2);
        assert_eq!(s * x + t * y, d);

        let (x, y) = (2, 52);
        let (d, s, t) = SnfCalc::gcdx(&x, &y);
        assert_eq!(d, 2);
        assert_eq!(s, 1);
        assert_eq!(t, 0);

        let (x, y) = (-2, 52);
        let (d, s, t) = SnfCalc::gcdx(&x, &y);
        assert_eq!(d, 2);
        assert_eq!(s, -1);
        assert_eq!(t, 0);
    }

    #[test]
    fn eliminate_row1() {
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.eliminate_row(0,0);

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res[[0, 0]], 1);
        assert_eq!(res[[0, 1]], 0);
        assert_eq!(res[[0, 2]], 0);
        assert!(p.is_id());
        assert!(pinv.is_id());
        assert_eq!(a.clone() * q, res);
        assert_eq!(res * qinv, a);
    }

    #[test]
    fn eliminate_row2() {
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.eliminate_row(1,1);

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res[[1, 0]], 0);
        assert_eq!(res[[1, 1]], 1);
        assert_eq!(res[[1, 2]], 0);
        assert!(p.is_id());
        assert!(pinv.is_id());
        assert_eq!(a.clone() * q, res);
        assert_eq!(res * qinv, a);
    }

    #[test]
    fn eliminate_col1() {
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.eliminate_col(0,0);

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res[[0, 0]], 1);
        assert_eq!(res[[1, 0]], 0);
        assert_eq!(res[[2, 0]], 0);
        assert_eq!(p * a.clone(), res);
        assert_eq!(pinv * res, a);
        assert!(q.is_id());
        assert!(qinv.is_id());
    }

    #[test]
    fn eliminate_col2() {
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.eliminate_col(1,1);

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res[[0, 1]], 0);
        assert_eq!(res[[1, 1]], 1);
        assert_eq!(res[[2, 1]], 0);
        assert_eq!(p * a.clone(), res);
        assert_eq!(pinv * res, a);
        assert!(q.is_id());
        assert!(qinv.is_id());
    }

    #[test]
    fn eliminate_at1() {
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.eliminate_at(0,0);

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res[[0, 0]], 1);
        assert_eq!(res[[0, 1]], 0);
        assert_eq!(res[[0, 2]], 0);
        assert_eq!(res[[1, 0]], 0);
        assert_eq!(res[[2, 0]], 0);
        assert_eq!(p * a.clone() * q, res);
        assert_eq!(pinv * res * qinv, a.clone());
    }

    #[test]
    fn eliminate_at2() {
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.eliminate_at(1,1);

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res[[1, 1]], 1);
        assert_eq!(res[[1, 0]], 0);
        assert_eq!(res[[1, 2]], 0);
        assert_eq!(res[[0, 1]], 0);
        assert_eq!(res[[2, 1]], 0);
        assert_eq!(p * a.clone() * q, res);
        assert_eq!(pinv * res * qinv, a.clone());
    }

    #[test]
    fn select_pivot() {
        let a = DnsMat::from(array![[1,0,1], [0,1,0], [0,1,1]]);
        let calc = SnfCalc::new(a.clone(), [true; 4]);

        assert_eq!(calc.select_pivot(0, 0), Some(0));
        assert_eq!(calc.select_pivot(1, 0), None);
        assert_eq!(calc.select_pivot(2, 0), None);
        assert_eq!(calc.select_pivot(0, 1), Some(1));
        assert_eq!(calc.select_pivot(1, 1), Some(1));
        assert_eq!(calc.select_pivot(2, 1), Some(2));
        assert_eq!(calc.select_pivot(0, 2), Some(0));
        assert_eq!(calc.select_pivot(1, 2), Some(2));
        assert_eq!(calc.select_pivot(2, 2), Some(2));
    }

    #[test]
    fn eliminate_all1() {
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.eliminate_all();

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res, DnsMat::from(array![[1,0,0],[0,3,0],[0,0,0]]));
        assert_eq!(p * a.clone() * q, res);
        assert_eq!(pinv * res * qinv, a.clone());
    }

    #[test]
    fn eliminate_all2() {
        let a = DnsMat::from(array![
            [1, 0, 1, 0, 0, 1, 1, 0, 1],
            [0, 1, 3, 1, 0, 1, 0, 2, 0],
            [0, 0, 1, 1, 0, 0, 0, 5, 1],
            [0, 1, 1, 0, 3, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 0, 1, 0, 1],
            [1, 0, 2, 0, 1, 1, 0, 1, 1]
        ]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.eliminate_all();

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res, DnsMat::diag((6, 9), vec![1,1,1,1,1,1]));
        assert_eq!(p * a.clone() * q, res);
        assert_eq!(pinv * res * qinv, a.clone());
    }

    #[test]
    fn eliminate_all3() {
        let a: DnsMat<i64> = DnsMat::from(array![
            [-20, -7, -27, 2, 29], 
            [17, 8, 14, -4, -10], 
            [13, 8, 10, -4, -6], 
            [-9, -2, -14, 0, 16], 
            [5, 0, 5, -1, -4]
        ]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.eliminate_all();

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );

        assert_eq!(res, DnsMat::diag((5, 5), vec![1,1,1,2,60]));
        assert_eq!(p * a.clone() * q, res);
        assert_eq!(pinv * res * qinv, a.clone());
    }

    #[test]
    fn diag_normalize1() {
        let a = DnsMat::diag((5, 5), vec![4, 24, -2, 1, 72]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.diag_normalize();

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );
        
        assert_eq!(p * a.clone() * q, res);
        assert_eq!(pinv * res * qinv, a.clone());
    }

    #[test]
    fn diag_normalize2() {
        let a = DnsMat::diag((5, 5), vec![0, -3, 54, 92, -4]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.diag_normalize();

        let (res, trans) = calc.result().destruct();
        let [p, pinv, q, qinv] = trans.map( |p| p.unwrap() );
        
        assert_eq!(p * a.clone() * q, res);
        assert_eq!(pinv * res * qinv, a.clone());
    }

    #[test]
    fn lll_preprocess_i64() { 
        use super::super::lll::tests::helper::assert_is_hnf;

        let a: DnsMat<i64> = DnsMat::from(array![
            [1, 0, 1, 0, 0, 1, 1, 0, 1],
            [0, 1, 3, 1, 0, 1, 0, 2, 0],
            [0, 0, 1, 1, 0, 0, 0, 5, 1],
            [0, 1, 1, 0, 3, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 0, 1, 0, 1],
            [1, 0, 2, 0, 1, 1, 0, 1, 1]
        ]);
        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.preprocess();

        let (res, trans) = calc.result().destruct();
        let [p, pinv, _, _] = trans.map( |p| p.unwrap() );

        assert_is_hnf(&res);

        assert_eq!(p * a.clone(), res);
        assert_eq!(pinv * res, a.clone());
    }

    #[test]
    fn lll_preprocess_bigint() { 
        use super::super::lll::tests::helper::assert_is_hnf;
        use num_bigint::BigInt;

        let a: DnsMat<BigInt> = DnsMat::from(array![
            [1, 0, 1, 0, 0, 1, 1, 0, 1],
            [0, 1, 3, 1, 0, 1, 0, 2, 0],
            [0, 0, 1, 1, 0, 0, 0, 5, 1],
            [0, 1, 1, 0, 3, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 0, 1, 0, 1],
            [1, 0, 2, 0, 1, 1, 0, 1, 1]
        ].map(|&a| a.into()));

        let mut calc = SnfCalc::new(a.clone(), [true; 4]);
        calc.preprocess();

        let (res, trans) = calc.result().destruct();
        let [p, pinv, _, _] = trans.map( |p| p.unwrap() );

        assert_is_hnf(&res);

        assert_eq!(p * a.clone(), res);
        assert_eq!(pinv * res, a.clone());
    }
}