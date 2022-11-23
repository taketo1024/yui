use core::panic;

use crate::math::traits::EucRing;
use super::DnsMat;

pub type SnfFlags = [bool; 4];
#[derive(Clone)]
pub struct SnfResult<R: EucRing> { 
    result: DnsMat<R>,
    p:    Option<DnsMat<R>>,
    pinv: Option<DnsMat<R>>,
    q:    Option<DnsMat<R>>,
    qinv: Option<DnsMat<R>>
}

pub fn snf<R: EucRing>(target: &DnsMat<R>, flags: SnfFlags) -> SnfResult<R> {
    let copy = target.clone();
    snf_in_place(copy, flags)
}

pub fn snf_in_place<R: EucRing>(target: DnsMat<R>, flags: SnfFlags) -> SnfResult<R> {
    let mut calc = SnfCalc::new(target, flags);
    
    calc.eliminate_all();

    calc.result()
}

// -- private -- //

struct SnfCalc<R: EucRing> { 
    target: DnsMat<R>,
    p:    Option<DnsMat<R>>,
    pinv: Option<DnsMat<R>>,
    q:    Option<DnsMat<R>>,
    qinv: Option<DnsMat<R>>
}

impl<R: EucRing> SnfCalc<R> {
    fn new(target: DnsMat<R>, flags: SnfFlags) -> Self { 
        let eye_opt = |size, flag| {
            if flag{ Some(DnsMat::eye(size)) } else { None }
        };

        let (m, n) = target.shape();
        let p    = eye_opt(m, flags[0]);
        let pinv = eye_opt(m, flags[1]);
        let q    = eye_opt(n, flags[2]);
        let qinv = eye_opt(n, flags[3]);

        SnfCalc{ target, p, pinv, q, qinv }
    }

    fn result(self) -> SnfResult<R> {
        SnfResult { result: self.target, p: self.p, pinv: self.pinv, q: self.q, qinv: self.qinv }
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

        // swap rows
        if i_p > i { 
            self.swap_rows(i, i_p);
        }

        // swap cols
        if j > i { 
            self.swap_cols(i, j);
        }

        // normalize pivot
        let (u, uinv) = self.target[[i, i]].normalizing_unit();
        if !u.is_one() { 
            self.mul_col(i, u, uinv);
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
    }

    fn swap_cols(&mut self, i: usize, j: usize) {
        self.target.swap_cols(i, j);
        self.q.as_mut().map( |q| q.swap_cols(i, j) );
        self.qinv.as_mut().map( |qinv| qinv.swap_rows(i, j) );
    }

    fn mul_row(&mut self, i: usize, u: R, uinv: R) {
        assert_eq!(u.clone() * uinv.clone(), R::one());
        self.target.mul_row(i, u.clone());
        self.p.as_mut().map( |p| p.mul_row(i, u));
        self.pinv.as_mut().map( |pinv| pinv.mul_col(i, uinv) );
    }
    
    fn mul_col(&mut self, i: usize, u: R, uinv: R) {
        assert_eq!(u.clone() * uinv.clone(), R::one());
        self.target.mul_col(i, u.clone());
        self.q.as_mut().map( |q| q.mul_col(i, u) );
        self.qinv.as_mut().map( |qinv| qinv.mul_row(i, uinv) );
    }

    // Multiply [a, b; c, d] from left, assuming det = 1.
    pub fn left_elementary(&mut self, comps: [R; 4], i: usize, j: usize) { 
        self.target.left_elementary(comps.clone(), i, j);
        self.p.as_mut().map( |p| p.left_elementary(comps.clone(), i, j) ); 
        self.pinv.as_mut().map(|pinv| { 
            let [a, b, c, d] = comps;
            let inv_t = [d, -c, -b, a];
            pinv.right_elementary(inv_t, i, j) 
        }); 
    }

    // Multiply [a, c; b, d] from right, assuming det = 1. 
    pub fn right_elementary(&mut self, comps: [R; 4], i: usize, j: usize) { 
        self.target.right_elementary(comps.clone(), i, j);
        self.q.as_mut().map( |q| q.right_elementary(comps.clone(), i, j) ); 
        self.qinv.as_mut().map(|qinv| { 
            let [a, b, c, d] = comps;
            let inv_t = [d, -c, -b, a];
            qinv.left_elementary(inv_t, i, j) 
        }); 
    }

    fn select_pivot(&self, below_i: usize, j: usize) -> Option<usize> { 
        // find row `i` below `below_i` with minimum nnz. 
        (below_i..self.target.nrows())
            .map( |i| (i, self.row_nz(i)) )
            .min_by( |e1, e2| e1.1.cmp(&e2.1) )
            .map( |e| e.0 )
    }

    fn eliminate_at(&mut self, i: usize, j: usize) {
        assert!(!self.target[[i, j]].is_zero());

        while self.row_nz(i) > 1 || self.col_nz(j) > 1 { 
            let modified = self.eliminate_row(i, j) || self.eliminate_col(i, j);
            if !modified {
                panic!("Detect endless loop");
            }
        }
    }

    fn eliminate_row(&mut self, i: usize, j: usize) -> bool { 
        let mut modified = false;

        for i1 in 0..self.target.nrows() {
            if i == i1 || self.target[[i1, j]].is_zero() { continue }

            // d = sx + ty,
            // a = x/d,
            // b = y/d.
        
            // [ s t][x] < i  = [d]
            // [-b a][y] < i1   [0]

            let x = self.target[[i , j]].clone();
            let y = self.target[[i1, j]].clone();

            let (d, s, t) = EucRing::gcdx(x.clone(), y.clone());
            let (a, b) = (x / d.clone(), y / d.clone());

            self.left_elementary([s, t, -b, a], i, i1);
            modified = true
        }
        
        modified
    }
    
    fn eliminate_col(&mut self, i: usize, j: usize) -> bool { 
        let mut modified = false;

        for j1 in 0..self.target.ncols() {
            if j == j1 || self.target[[i, j1]].is_zero() { continue }

            // d = sx + ty,
            // a = x/d,
            // b = y/d.
        
            // [x y][s -b] = [d 0]
            //      [t  a]   

            let x = self.target[[i, j ]].clone();
            let y = self.target[[i, j1]].clone();

            let (d, s, t) = EucRing::gcdx(x.clone(), y.clone());
            let (a, b) = (x / d.clone(), y / d.clone());

            self.right_elementary([s, t, -b, a], j, j1);
            modified = true
        }
        
        modified
    }
    
}

#[cfg(test)]
mod tests {
    use ndarray::array;

    use super::*;

    #[test]
    fn init() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6]]);
        let calc = SnfCalc::new(a, [true; 4]);

        assert_eq!(calc.target, DnsMat::from(array![[1,2,3], [4,5,6]]));
        assert_eq!(calc.p,    Some(DnsMat::eye(2)));
        assert_eq!(calc.pinv, Some(DnsMat::eye(2)));
        assert_eq!(calc.q,    Some(DnsMat::eye(3)));
        assert_eq!(calc.qinv, Some(DnsMat::eye(3)));
    }

    #[test]
    fn init_no_pq() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6]]);
        let calc = SnfCalc::new(a, [false; 4]);
        assert_eq!(calc.p,    None);
        assert_eq!(calc.pinv, None);
        assert_eq!(calc.q,    None);
        assert_eq!(calc.qinv, None);
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
        let mut calc = SnfCalc::new(a, [true; 4]);
        calc.swap_rows(0, 1);

        assert_eq!(calc.target, DnsMat::from(array![[4,5,6], [1,2,3], [7,8,9]]));
        assert_eq!(calc.p,    Some(DnsMat::from(array![[0,1,0], [1,0,0], [0,0,1]])));
        assert_eq!(calc.pinv, Some(DnsMat::from(array![[0,1,0], [1,0,0], [0,0,1]])));
        assert!(calc.q   .map(|x| x.is_eye()).unwrap_or(false));
        assert!(calc.qinv.map(|x| x.is_eye()).unwrap_or(false));
    }

    #[test]
    fn swap_cols() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let mut calc = SnfCalc::new(a, [true; 4]);
        calc.swap_cols(0, 1);

        assert_eq!(calc.target, DnsMat::from(array![[2,1,3], [5,4,6], [8,7,9]]));
        assert_eq!(calc.q,    Some(DnsMat::from(array![[0,1,0], [1,0,0], [0,0,1]])));
        assert_eq!(calc.qinv, Some(DnsMat::from(array![[0,1,0], [1,0,0], [0,0,1]])));
        assert!(calc.p   .map(|x| x.is_eye()).unwrap_or(false));
        assert!(calc.pinv.map(|x| x.is_eye()).unwrap_or(false));
    }

    #[test]
    fn mul_row() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let mut calc = SnfCalc::new(a, [true; 4]);
        calc.mul_row(0, -1, -1);

        assert_eq!(calc.target, DnsMat::from(array![[-1,-2,-3], [4,5,6], [7,8,9]]));
        assert_eq!(calc.p,    Some(DnsMat::from(array![[-1,0,0], [0,1,0], [0,0,1]])));
        assert_eq!(calc.pinv, Some(DnsMat::from(array![[-1,0,0], [0,1,0], [0,0,1]])));
        assert!(calc.q   .map(|x| x.is_eye()).unwrap_or(false));
        assert!(calc.qinv.map(|x| x.is_eye()).unwrap_or(false));
    }

    #[test]
    fn mul_col() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let mut calc = SnfCalc::new(a, [true; 4]);
        calc.mul_col(0, -1, -1);

        assert_eq!(calc.target, DnsMat::from(array![[-1,2,3], [-4,5,6], [-7,8,9]]));
        assert!(calc.p   .map(|x| x.is_eye()).unwrap_or(false));
        assert!(calc.pinv.map(|x| x.is_eye()).unwrap_or(false));
        assert_eq!(calc.q,    Some(DnsMat::from(array![[-1,0,0], [0,1,0], [0,0,1]])));
        assert_eq!(calc.qinv, Some(DnsMat::from(array![[-1,0,0], [0,1,0], [0,0,1]])));
    }

    #[test]
    fn left_elementary() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let e = [3,2,4,3]; // det = 1
        let mut calc = SnfCalc::new(a, [true; 4]);
        calc.left_elementary(e, 0, 1);

        assert_eq!(calc.target, DnsMat::from(array![[11,16,21], [16,23,30], [7,8,9]]));
        assert_eq!(calc.p,    Some(DnsMat::from(array![[3,2,0], [4,3,0], [0,0,1]])));
        assert_eq!(calc.pinv, Some(DnsMat::from(array![[3,-2,0], [-4,3,0], [0,0,1]])));
        assert!(calc.q   .map(|x| x.is_eye()).unwrap_or(false));
        assert!(calc.qinv.map(|x| x.is_eye()).unwrap_or(false));
    }
    #[test]
    fn right_elementary() { 
        let a = DnsMat::from(array![[1,2,3], [4,5,6], [7,8,9]]);
        let e = [3,2,4,3]; // det = 1
        let mut calc = SnfCalc::new(a, [true; 4]);
        calc.right_elementary(e, 0, 1);

        assert_eq!(calc.target, DnsMat::from(array![[7,10,3], [22,31,6], [37,52,9]]));
        assert!(calc.p   .map(|x| x.is_eye()).unwrap_or(false));
        assert!(calc.pinv.map(|x| x.is_eye()).unwrap_or(false));
        assert_eq!(calc.q,    Some(DnsMat::from(array![[3,4,0], [2,3,0], [0,0,1]])));
        assert_eq!(calc.qinv, Some(DnsMat::from(array![[3,-4,0], [-2,3,0], [0,0,1]])));
    }
}