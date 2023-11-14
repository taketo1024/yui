// "Extended GCD and Hermite Normal Form Algorithms via Lattice Basis Reduction",
// George Havas, Bohdan S. Majewski, and Keith R. Matthews
// https://projecteuclid.org/journals/experimental-mathematics/volume-7/issue-2/Extended-GCD-and-Hermite-normal-form-algorithms-via-lattice-basis/em/1048515660.full
//
// see also: "Keith Matthews' LLL page", 
// http://www.numbertheory.org/lll.html
//
// "A generalization of the LLL-algorithm over euclidean rings or orders",
// Huguette Napias
// https://www.jstor.org/stable/43974220

use std::fmt::Debug;
use ndarray::{Array1, ArrayView1, Array2, ArrayView2};
use log::{trace, info};
use num_bigint::BigInt;

use yui::{Ring, RingOps, EucRing, EucRingOps, DivRound, Integer, IntOps};
use yui::{QuadInt, GaussInt, EisenInt};
use crate::dense::*;

pub fn lll<R>(b: &Mat<R>, with_trans: bool) -> (Mat<R>, Option<Mat<R>>)
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    lll_in_place(b.clone(), with_trans)
}

pub fn lll_in_place<R>(b: Mat<R>, with_trans: bool) -> (Mat<R>, Option<Mat<R>>)
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    info!("lll: {:?}", b.shape());
    
    let mut calc = LLLCalc::new(b, with_trans);
    calc.process();
    calc.result()
}

pub fn lll_hnf<R>(b: &Mat<R>, with_trans: [bool; 2]) -> (Mat<R>, Option<Mat<R>>, Option<Mat<R>>)
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    lll_hnf_in_place(b.clone(), with_trans)
}

pub fn lll_hnf_in_place<R>(b: Mat<R>, with_trans: [bool; 2]) -> (Mat<R>, Option<Mat<R>>, Option<Mat<R>>)
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    info!("lll-hnf: {:?}", b.shape());
    
    let mut calc = LLLHNFCalc::new(b, with_trans);
    calc.process();
    calc.result()
}

// -- LLLRing trait -- //

pub trait LLLRingOps<T>: EucRingOps<T> {}

pub trait LLLRing: EucRing + LLLRingOps<Self> + DivRound
where for<'x> &'x Self: LLLRingOps<Self> {
    type Int: PartialOrd + Ord + Debug;
    fn alpha() -> (Self, Self);
    fn as_int(&self) -> Option<Self::Int>;
    fn conj(&self) -> Self;
    fn norm(&self) -> Self { 
        self * self.conj()
    }
}

// -- implementations of LLLRing -- //

macro_rules! impl_for_int {
    ($type:ty) => {
        impl LLLRingOps<Self> for $type {}
        impl<'a> LLLRingOps<$type> for &'a $type {}

        impl LLLRing for $type {
            type Int = Self;
            fn alpha() -> (Self, Self) {
                (Self::from(3), Self::from(4))
            }

            fn as_int(&self) -> Option<Self::Int> { 
                Some(self.clone())
            }
            
            fn conj(&self) -> Self { 
                self.clone()
            }
        }
    };
}

impl_for_int!(i32);
impl_for_int!(i64);
impl_for_int!(i128);
impl_for_int!(BigInt);

macro_rules! impl_for_quad_int {
    ($type:ident, $p:literal, $q:literal) => {
        impl<I> LLLRingOps<Self> for $type<I>
        where I: Integer, for<'x> &'x I: IntOps<I> {}

        impl<'a, I> LLLRingOps<$type<I>> for &'a $type<I>
        where I: Integer, for<'x> &'x I: IntOps<I> {}

        impl<I> LLLRing for $type<I>
        where I: Integer, for<'x> &'x I: IntOps<I> {
            type Int = I;

            fn alpha() -> (Self, Self) {
                (Self::from_real(I::from($p)), Self::from_real(I::from($q)))
            }

            fn as_int(&self) -> Option<Self::Int> {
                if self.right().is_zero() { 
                    Some(self.left().clone())
                } else { 
                    None
                }
            }

            fn conj(&self) -> Self {
                QuadInt::conj(self)
            }
        }
    }
}

impl_for_quad_int!(GaussInt, 3, 4);
impl_for_quad_int!(EisenInt, 2, 3);

// -- private implementation -- //

type Row = usize;
type Col = usize;

#[derive(Debug)]
struct LLLCalc<R>
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    data: LLLData<R>
}

impl<R> LLLCalc<R>
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    fn new(target: Mat<R>, with_trans: bool) -> Self {
        let mut data = LLLData::new(target, [with_trans, false]);
        data.setup();

        LLLCalc { data }
    }

    pub fn process(&mut self) { 
        assert!(self.data.step == 1);
        let m = self.data.rows();

        while self.data.step < m { 
            self.iterate()
        }
    }

    fn iterate(&mut self) { 
        trace!("step: {}.\n{}", self.data.step, self.data.dump());

        let k = self.data.step;

        self.data.reduce(k - 1, k);

        if self.data.lovasz_ok(k) { 
            for i in (0..k-1).rev() { 
                self.data.reduce(i, k);
            }
            self.data.next();
        } else { 
            self.data.swap(k);
            self.data.back();
        }
    }

    pub fn result(self) -> (Mat<R>, Option<Mat<R>>) { 
        let (target, p, _) = self.data.result();
        (target, p)
    }
}

#[derive(Debug)]
struct LLLHNFCalc<R>
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    data: LLLData<R>
}

impl<R> LLLHNFCalc<R>
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    fn new(target: Mat<R>, with_trans: [bool; 2]) -> Self { 
        let data = LLLData::new(target, with_trans);
        LLLHNFCalc { data }
    }

    fn process(&mut self) { 
        assert!(self.data.step > 0);
        let m = self.data.rows();

        while self.data.step < m { 
            self.iterate();
        }
    }

    fn iterate(&mut self) { 
        trace!("step: {}.\n{}", self.data.step, self.data.dump());

        let k = self.data.step;

        self.reduce(k - 1, k);

        if self.is_ok(k) { 
            for i in (0..k-1).rev() { 
                self.reduce(i, k);
            }
            self.data.next();
        } else { 
            self.data.swap(k);
            self.data.back();
        }
    }

    fn result(self) -> (Mat<R>, Option<Mat<R>>, Option<Mat<R>>) { 
        let m = self.data.rows();
        let (mut target, mut p, mut pinv) = self.data.result();

        for i in 0..m/2 {
            let j = m - i - 1;
            if i == j { break } 

            target.swap_rows(i, j);
            if let Some(p) = p.as_mut() { 
                p.swap_rows(i, j) 
            }
            if let Some(pinv) = pinv.as_mut() { 
                pinv.swap_cols(i, j) 
            }
        }

        (target, p, pinv)
    }

    fn reduce(&mut self, i: Row, k: Row) {
        assert!(i < k);

        let j = self.data.nz_col_in(i);

        if let Some(j) = j { 
            let u = self.data.target[[i, j]].normalizing_unit();
            if !u.is_one() { 
                self.data.mul_row(i, &u);
            }

            let a = &self.data.target;
            let a0 = &a[[i, j]];
            let a1 = &a[[k, j]];
            let q = a1.div_round(a0);

            if !q.is_zero() { 
                self.data.add_row_to(i, k, &-q) // a[k] -= q & a[i]
            }
        } else { 
            self.data.reduce(i, k)
        }
    }

    fn is_ok(&self, k: usize) -> bool {
        assert!(k > 0);

        let j = self.data.nz_col_in(k - 1);
        let l = self.data.nz_col_in(k);

        match (j, l) { 
            (Some(j), Some(l)) => j > l,
            (Some(_), None)    => false,
            (None,    Some(_)) => true,
            (None,    None)    => self.data.lovasz_ok(k)
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
struct LLLData<R>
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    target: Mat<R>,
    p:    Option<Mat<R>>,
    pinv: Option<Mat<R>>,
    det: Vec<R>,        // D[i] = det(b_1, ..., b_i)^2 = Π^i |b^*_j|^2. 
    lambda: Mat<R>,  // l[i,j] = D[j] * p_ij (0 <= j < i)
    step: usize,
}

impl<R> LLLData<R>
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    fn new(target: Mat<R>, flags: [bool; 2]) -> Self { 
        let m = target.rows();
        let p =    if flags[0] { Some(Mat::id(m)) } else { None };
        let pinv = if flags[1] { Some(Mat::id(m)) } else { None };
        let det = vec![R::one(); m];
        let lambda = Mat::zero((m, m)); // lower-triangular

        LLLData { target, p, pinv, det, lambda, step: 1 }
    }

    fn result(self) -> (Mat<R>, Option<Mat<R>>, Option<Mat<R>>) {
        (self.target, self.p, self.pinv)
    }

    fn setup(&mut self) { 
        let b = self.target.array();
        let (_, l, d) = orthogonalize(&b.view());

        self.lambda = Mat::from(l);
        self.det = d;
    }

    // Lovasz condition:
    //   |b^*_k|^2 >= (α - m[k,k-1]^2)|b^*_{k-1}|^2
    // ⇔ |b^*{k-1}|^{-2} * |b^*_k|^2 + m[k,k-1]^2 >= α
    // ⇔ (d[k-2]/d[k-1]) * (d[k]/d[k-1]) + (l[k,k-1]/d[k-1])^2 >= p/q
    // ⇔ q * (d[k-2] * d[k] + (l[k,k-1])^2) >= p d[k-1]^2
    
    fn lovasz_ok(&self, k: usize) -> bool { 
        assert!(k > 0);

        let (d, l) = (&self.det, &self.lambda);
        let (p, q) = &R::alpha();

        let one = R::one();
        let d0 = if k >= 2 { &d[k - 2] } else { &one };
        let d1 = &d[k - 1];
        let d2 = &d[k];
        let l0 = &l[[k, k - 1]];

        let lhs = q * (d0 * d2 + l0.norm());
        let rhs = p * (d1 * d1);

        let Some(lhs) = lhs.as_int() else { panic!() };
        let Some(rhs) = rhs.as_int() else { panic!() };

        lhs >= rhs
    }

    fn reduce(&mut self, i: Row, k: Row) { 
        assert!(i < k);

        let (d, l) = (&self.det, &self.lambda);
        let l_ki = &l[[k, i]];
        let d_i = &d[i];
        let q = l_ki.div_round(d_i);

        if !q.is_zero() { 
            self.add_row_to(i, k, &-q) // a[k] -= q & a[i]
        }
    }

    fn swap(&mut self, k: Row) { 
        assert!(k > 0);

        // b[k-1, ..] <--> b[k, ..]
        self.target.swap_rows(k - 1, k);
        if let Some(p) = self.p.as_mut() { 
            p.swap_rows(k-1, k) 
        }
        if let Some(pinv) = self.pinv.as_mut() { 
            pinv.swap_cols(k-1, k) 
        }

        //                   k 
        //      |                     |
        //  k-1 |  . . . . 0          |
        //  k   |  . . . . * 0        |
        //  k+1 |          . . 0      |
        //      |          . .        |
        //      |          . .        |
        //

        // λ[k-1, ..] <--> λ[k, ..] 
        for j in 0..k-1 { 
            let slice = ndarray::s![.., j];
            let l_j = &mut self.lambda.array_mut().slice_mut(slice);
            l_j.swap(k - 1, k);
        }

        let one = R::one();
        let d = &self.det;
        let d0 = if k >= 2 { &d[k - 2] } else { &one };
        let d1 = &d[k - 1];
        let d2 = &d[k];

        let m = self.lambda.cols();

        // λ[.., k-1] <--> λ[.., k]
        for i in k+1..m { 
            let l = &self.lambda;
            let l0 = &l[[k, k-1]];
            let l1 = &l[[i, k-1]];
            let l2 = &l[[i, k]];

            let s = l0.conj() * l1 + l2 * d0;
            let t = l1 * d2 - l2 * l0;

            self.lambda[[i, k-1]] = &s / d1;
            self.lambda[[i, k]]   = &t / d1;
        }

        // λ[k, k-1] remains unchanged.

        let l0 = &self.lambda[[k,k-1]];

        self.det[k-1] = (d0 * d2 + l0.norm()) / d1;
        self.lambda[[k, k-1]] = l0.conj();

        trace!("swap-rows ({},{}).\n{}", k-1, k, self.dump());
    }

    fn mul_row(&mut self, i: Row, r: &R) { 
        assert!(r.is_unit());

        self.target.mul_row(i, r);
        if let Some(p) = self.p.as_mut() { 
            p.mul_row(i, r) 
        }

        if let Some(pinv) = self.pinv.as_mut() { 
            let rinv = r.inv().unwrap();
            pinv.mul_col(i, &rinv) 
        }

        self.lambda.mul_row(i, r);
        self.lambda.mul_col(i, &r.conj());

        trace!("mul {} to row {}.\n{}", r, i, self.dump());
    }

    fn add_row_to(&mut self, i: Row, k: Row, r: &R) {
        assert!(i < k);
        
        self.target.add_row_to(i, k, r);
        if let Some(p) = self.p.as_mut() { 
            p.add_row_to(i, k, r) 
        }

        if let Some(pinv) = self.pinv.as_mut() { 
            let nr = -r;
            pinv.add_col_to(k, i, &nr) 
        }

        self.lambda[[k, i]] += r * &self.det[i];

        for j in 0..i { 
            let a = r * &self.lambda[[i, j]];
            self.lambda[[k, j]] += a;
        }

        trace!("add-row {} to {}, mul {}.\n{}", i, k, r, self.dump());
    }

    fn nz_col_in(&self, i: Row) -> Option<Col> {
        self.target.array().row(i).iter().enumerate().filter_map(|(j, a)| { 
            if !a.is_zero() { Some(j) } else { None }
        }).next()
    }

    fn next(&mut self) { 
        self.step += 1;
    }

    fn back(&mut self) { 
        if self.step > 1 { 
            self.step -= 1;
        }
    }

    fn rows(&self) -> usize { 
        self.target.rows()
    }

    fn dump(&self) -> String { 
        format!("{}", self.target)
        // format!("target:\n{},\nlambda:\n{},\ndet: {:?}.", self.target, self.lambda, self.det)
    }
}

// -- helper funcs -- //

fn orthogonalize<R>(b: &ArrayView2<R>) -> (Array2<R>, Array2<R>, Vec<R>)
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    let m = b.nrows();
    let one = R::one();

    let mut c = b.to_owned();
    let mut d = vec![R::one(); m];
    let mut l = Array2::zeros((m, m));

    d[0] = h_dot(&c.row(0), &c.row(0));

    for i in 1..m { 
        for j in 0..i { 
            // c_i = (d[j] * c_i - l[i,j] * c_j) / d[j - 1];
            
            let l0 = h_dot(&b.row(i), &c.row(j));
            let d0 = if j > 0 { &d[j - 1] } else { &one };
            let d1 = &d[j];
            
            let c_i = smul(d1, &c.row(i)) - smul(&l0, &c.row(j));
            let c_i = sdiv(&c_i.view(), d0);
            
            l[[i, j]] = l0;
            c.row_mut(i).assign(&c_i);
        }

        let c_i = &c.row(i);
        let d0 = &d[i - 1];
        
        d[i] = &h_dot(c_i, c_i) / d0;
    }

    (c, l, d)
}

fn smul<'a, R>(r: &'a R, rhs: &ArrayView1<'a, R>) -> Array1<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    rhs.map(|x| r * x)
}

fn sdiv<'a, R>(lhs: &ArrayView1<'a, R>, r: &'a R) -> Array1<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    lhs.map(|x| x / r)
}

fn h_dot<'a, R>(lhs: &ArrayView1<'a, R>, rhs: &ArrayView1<'a, R>) -> R
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    assert_eq!(lhs.dim(), rhs.dim());
    ndarray::Zip::from(lhs).and(rhs).fold(R::zero(), |mut acc, a, b| { 
        acc += a * b.conj();
        acc
    })
}

#[cfg(test)]
pub(super) mod tests {
    use super::*;
    use ndarray::array;
 
    #[test]
    fn test_large_orth_basis() {
        let a = array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ];
        let (c, l, d) = orthogonalize(&a.view());

        assert_eq!(c.row(0), array![1, -1, 3].view());
        assert_eq!(c.row(1), array![-5, 16, 7].view());
        assert_eq!(c.row(2), array![15, 6, -3].view());

        assert_eq!(d[0], 11);
        assert_eq!(d[1], 30);
        assert_eq!(d[2], 9);

        assert_eq!(l, array![
            [0,  0, 0],
            [16, 0, 0],
            [17,69, 0]
        ]);
    }

    #[test]
    fn data_init() { 
        let a = Mat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);

        let data = LLLData::new(a.clone(), [false, false]);
        assert_eq!(data.target, a);
        assert_eq!(data.p, None);
        assert_eq!(data.pinv, None);

        let data = LLLData::new(a.clone(), [true, false]);
        assert_eq!(data.target, a);
        assert_eq!(data.p, Some(Mat::id(3)));
        assert_eq!(data.pinv, None);

        let data = LLLData::new(a.clone(), [false, true]);
        assert_eq!(data.target, a);
        assert_eq!(data.p, None);
        assert_eq!(data.pinv, Some(Mat::id(3)));
    }
    
    #[test]
    fn setup() { 
        let a = Mat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut data = LLLData::new(a, [false, false]);

        data.setup();

        assert_eq!(data.det.len(), 3);
        assert_eq!(data.det[0], 11);
        assert_eq!(data.det[1], 30);
        assert_eq!(data.det[2], 9);

        assert_eq!(data.lambda, Mat::from(array![
            [0,  0, 0],
            [16, 0, 0],
            [17,69, 0]
        ]));
    }

    #[test]
    fn swap() { 
        let a0 = Mat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut data0 = LLLData::new(a0, [false, false]);
        data0.setup();
        data0.swap(1);
        data0.swap(2);
        
        // compare data
        let a1 = Mat::from(array![
            [1, 0, 5],
            [1, 2, 6],
            [1,-1, 3]
        ]);
        let mut data1 = LLLData::new(a1, [false, false]);
        data1.setup();

        assert_eq!(data0, data1);
    }

    #[test]
    fn swap_trans() { 
        let a0 = Mat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut data0 = LLLData::new(a0.clone(), [true, true]);
        data0.setup();
        data0.swap(1);
        data0.swap(2);

        // compare data
        let a1 = Mat::from(array![
            [1, 0, 5],
            [1, 2, 6],
            [1,-1, 3]
        ]);
        
        let p = data0.p.unwrap().clone();
        let pinv = data0.pinv.unwrap().clone();

        assert_eq!(p * a0.clone(), a1.clone());
        assert_eq!(pinv * a1.clone(), a0.clone());
    }

    #[test]
    fn add_row_to() { 
        let a0 = Mat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut data0 = LLLData::new(a0, [false, false]);
        data0.setup();
        data0.add_row_to(0, 1, &2);
        data0.add_row_to(1, 2, &-3);
        
        // compare data
        let a1 = Mat::from(array![
            [1,-1, 3],
            [3,-2,11],
            [-8,8,-27]
        ]);
        let mut data1 = LLLData::new(a1, [false, false]);
        data1.setup();

        assert_eq!(data0, data1);
    }

    #[test]
    fn add_row_to_trans() { 
        let a0 = Mat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut data0 = LLLData::new(a0.clone(), [true, true]);
        data0.setup();
        data0.add_row_to(0, 1, &2);
        data0.add_row_to(1, 2, &-3);
        
        // compare data
        let a1 = Mat::from(array![
            [1,-1, 3],
            [3,-2,11],
            [-8,8,-27]
        ]);

        let p = data0.p.unwrap().clone();
        let pinv = data0.pinv.unwrap().clone();
 
        assert_eq!(p * a0.clone(), a1.clone());
        assert_eq!(pinv * a1.clone(), a0.clone());
     }

    #[test]
    fn mul_row() { 
        let a0 = Mat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut data0 = LLLData::new(a0, [false, false]);
        data0.setup();
        data0.mul_row(1, &-1);
        data0.mul_row(2, &-1);
        
        // compare data
        let a1 = Mat::from(array![
            [1,-1, 3],
            [-1, 0, -5],
            [-1, -2, -6]
        ]);
        let mut data1 = LLLData::new(a1, [false, false]);
        data1.setup();

        assert_eq!(data0, data1);
    }

    #[test]
    fn mul_row_trans() { 
       let a0 = Mat::from(array![
           [1,-1, 3],
           [1, 0, 5],
           [1, 2, 6]
       ]);
       let mut data0 = LLLData::new(a0.clone(), [true, true]);
       data0.setup();
       data0.mul_row(1, &-1);
       data0.mul_row(2, &-1);
       
       // compare data
       let a1 = Mat::from(array![
           [1,-1, 3],
           [-1, 0, -5],
           [-1, -2, -6]
       ]);

       let p = data0.p.unwrap().clone();
       let pinv = data0.pinv.unwrap().clone();

       assert_eq!(p * a0.clone(), a1.clone());
       assert_eq!(pinv * a1.clone(), a0.clone());
    }

    #[test]
    fn lll() { 
        let a = Mat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut calc = LLLCalc::new(a.clone(), true);
        calc.process();

        let (res, Some(p)) = calc.result() else { panic!() };

        assert_eq!(res, Mat::from(array![
            [0, 1, -1],
            [1, 0, -1],
            [1, 1, 1]
        ]));

        helper::assert_is_reduced( &res.array().view() );

        assert_eq!(p * a, res);
    }

     #[test]
     fn lll_gcdx() { 
        // MEMO: γ = 10
        let a = Mat::from(array![
            [1, 0, 0, 40],
            [0, 1, 0, 60],
            [0, 0, 1, 90]
        ]);
        let mut calc = LLLCalc::new(a.clone(), true);
        calc.process();

        let (res, Some(p)) = calc.result() else { panic!() };

        assert_eq!(res, Mat::from(array![
            [3, -2, 0, 0],
            [0, 3, -2, 0],
            [-2, 0, 1, 10]
        ]));

        helper::assert_is_reduced( &res.array().view() );

        assert_eq!(p * a, res);
      }

     #[test]
     fn hnf() { 
        let a: Mat<i64> = Mat::from(array![
            [8,    44,   43],
            [4,    10,   43],
            [56, -550, -328],
            [76,   10,   42]
        ]);
        let mut calc = LLLHNFCalc::new(a.clone(), [true, true]);
        calc.process();

        let (res, Some(p), Some(pinv)) = calc.result() else { panic!() };

        helper::assert_is_hnf(&res);
        
        assert_eq!(p.clone() * a, res);
        assert_eq!(p * pinv, Mat::id(4));
    }

    #[test]
    fn setup_gauss() { 
        type A = GaussInt<i64>;

        let i = A::new;
        let a = Mat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(3, 3), i(-2, 4), i(6, 2)],
            [i(2, 2), i(-8, 0), i(-9, 1)],
        ]);
        let mut data = LLLData::new(a, [false, false]);

        data.setup();

        assert_eq!(data.det.len(), 3);
        assert_eq!(data.det[0], i(129, 0));
        assert_eq!(data.det[1], i(7436, 0));
        assert_eq!(data.det[2], i(161408, 0));

        assert_eq!(data.lambda, Mat::from(array![
            [i(0, 0),     i(0, 0),      i(0, 0)],
            [i(49, 15),   i(0, 0),      i(0, 0)],
            [i(-114, 48), i(1770,3162), i(0, 0)]
        ]));
    }

    #[test]
    fn swap_gauss() { 
        type A = GaussInt<i64>;
        let i = A::new;
        let a0 = Mat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(3, 3), i(-2, 4), i(6, 2)],
            [i(2, 2), i(-8, 0), i(-9, 1)],
        ]);
        let mut data0 = LLLData::new(a0, [false, false]);
        data0.setup();
        data0.swap(1);
        data0.swap(2);
        
        // compare data
        let a1 = Mat::from(array![
            [i(3, 3), i(-2, 4), i(6, 2)],
            [i(2, 2), i(-8, 0), i(-9, 1)],
            [i(-2, 3), i(7, 3), i(7, 3)],
        ]);
        let mut data1 = LLLData::new(a1, [false, false]);
        data1.setup();

        assert_eq!(data0, data1);
    }

    #[test]
    fn add_row_to_gauss() { 
        type A = GaussInt<i64>;
        let i = A::new;
        let a0 = Mat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(3, 3), i(-2, 4), i(6, 2)],
            [i(2, 2), i(-8, 0), i(-9, 1)],
        ]);
        let mut data0 = LLLData::new(a0, [false, false]);
        data0.setup();
        data0.add_row_to(0, 1, &i(1, 1));
        data0.add_row_to(1, 2, &i(-3, 2));
        
        // compare data
        let a1 = Mat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(-2, 4), i(2, 14), i(10, 12)],
            [i(0, -14), i(-42, -38), i(-63, -15)],
        ]);
        let mut data1 = LLLData::new(a1, [false, false]);
        data1.setup();

        assert_eq!(data0, data1);
    }
 
    #[test]
    fn mul_row_gauss() { 
        type A = GaussInt<i64>;
        let i = A::new;
        let a0 = Mat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(3, 3), i(-2, 4), i(6, 2)],
            [i(2, 2), i(-8, 0), i(-9, 1)],
        ]);

        let mut data0 = LLLData::new(a0.clone(), [false, false]);
        data0.setup();
        data0.mul_row(1, &i(0, 1));
        data0.mul_row(2, &i(0, -1));

        // compare data
        let a1 = Mat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(-3, 3), i(-4, -2), i(-2, 6)],
            [i(2, -2), i(0, 8), i(1, 9)],
        ]);
        let mut data1 = LLLData::new(a1, [false, false]);
        data1.setup();

        assert_eq!(data0, data1);
    }

    #[test]
    fn hnf_gauss() { 
        type A = GaussInt<i64>;
        let i = A::new;

        let a: Mat<A> = Mat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(3, 3), i(-2, 4), i(6, 2)],
            [i(2, 2), i(-8, 0), i(-9, 1)],
        ]);

        let mut calc = LLLHNFCalc::new(a.clone(), [true, true]);
        calc.process();

        let (res, Some(p), Some(pinv)) = calc.result() else { panic!() };

        helper::assert_is_hnf(&res);
        
        assert_eq!(p.clone() * a, res);
        assert_eq!(p * pinv, Mat::id(3));
    }

    #[test]
    fn hnf_eisen() { 
        type A = EisenInt<i64>;
        let i = A::new;

        let a: Mat<A> = Mat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(3, 3), i(-2, 4), i(6, 2)],
            [i(2, 2), i(-8, 0), i(-9, 1)],
        ]);

        let mut calc = LLLHNFCalc::new(a.clone(), [true, true]);
        calc.process();

        let (res, Some(p), Some(pinv)) = calc.result() else { panic!() };

        helper::assert_is_hnf(&res);

        assert_eq!(p.clone() * a, res);
        assert_eq!(p * pinv, Mat::id(3));
    }

    #[test]
    fn hnf_rand() {
        let d = 0.5;
        let shape = (8, 8);
        let a = crate::sparse::SpMat::<i64>::rand(shape, d).to_dense();

        let mut calc = LLLHNFCalc::new(a.clone(), [true, true]);
        calc.process();

        let (res, Some(p), Some(pinv)) = calc.result() else { panic!() };

        helper::assert_is_hnf(&res);
        
        assert_eq!(p.clone() * a, res);
        assert_eq!(p * pinv, Mat::id(shape.0));
    }

    pub(in super::super) mod helper { 
        use super::*;
        use std::ops::Div;
        use yui::Ratio;

        pub fn assert_is_hnf<R>(b: &Mat<R>)
        where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
            use ndarray::s;

            let (m, n) = b.shape();
            let b = b.array();

            let mut j0 = 0;
            for i0 in 0..m {
                let j1 = (0..n).find(|&j| !b[[i0, j]].is_zero()).unwrap_or(n);

                // assert: b[i0.., j0..j1] = 0
                let s = s![i0..m, j0..j1];
                assert!( b.slice(s).iter().all(|x| x.is_zero()) );

                if j1 < n { 
                    j0 = j1;

                    // assert: b[i0+1.., j0] = 0
                    let s = s![i0+1..m, j0];
                    assert!( b.slice(s).iter().all(|x| x.is_zero()) );

                    let a = &b[[i0, j0]];
                    assert!( !a.is_zero() );
                    assert!( a.normalizing_unit().is_one() );

                    let a0 = (a * &a.conj()).as_int().unwrap();
                    let s = s![0..i0, j0];
                    assert!( b.slice(s).iter().all(|x| {
                        let x0 = (x * &x.conj()).as_int().unwrap();
                        x0 < a0
                    }));
                } else {
                    break
                }
            }
        }
    
        pub fn assert_is_reduced<R>(b: &ArrayView2<R>)
        where R: Integer + LLLRing, for<'x> &'x R: IntOps<R> + LLLRingOps<R> {
            let m = b.nrows();
    
            let (c, l) = gram_schmidt(b);
            let alpha = Ratio::from(R::alpha());
            let thr = Ratio::new(R::one(), R::from(2));
    
            let size_reduced = l.iter().all(|r| r.abs() <= thr);
            let lovasz_ok = (1..m).all(|i| {
                let c0 = &c.row(i - 1);
                let c1 = &c.row(i);
                let m = &l[[i, i - 1]];
                is_lovasz_ok(c0, c1, m, &alpha)
            });
    
            assert!(size_reduced);
            assert!(lovasz_ok);
        }
    
        fn is_lovasz_ok<R>(c0: &ArrayView1<Ratio<R>>, c1: &ArrayView1<Ratio<R>>, m: &Ratio<R>, alpha: &Ratio<R>) -> bool
        where R: Integer, for<'x> &'x R: IntOps<R> {
            let r0 = dot::<Ratio<R>>(c0, c0);
            let r1 = dot::<Ratio<R>>(c1, c1);
            r1 >= (alpha - (m * m)) * r0
        }
    
        fn gram_schmidt<R>(b: &ArrayView2<R>) -> (Array2<Ratio<R>>, Array2<Ratio<R>>)
        where R: Integer, for<'x> &'x R: IntOps<R> {
            let m = b.nrows();
    
            let mut c = b.map(|x| Ratio::from_numer(x.clone()));
            let mut l = Array2::zeros((m, m));
    
            for i in 1..m { 
                for j in 0..i {
                    let (c_i, c_j) = (c.row(i), c.row(j));
                    let p_ij = proj_coeff::<Ratio<R>>(&c_j, &c_i);
                    let v_ij = smul::<Ratio<R>>(&p_ij, &c_j);
    
                    let mut c_i = c.row_mut(i);
                    c_i -= &v_ij;
    
                    l[[i, j]] = p_ij;
                }
            }
    
            (c, l)
        }
    
        fn proj_coeff<'a, R>(base: &ArrayView1<'a, R>, other: &ArrayView1<'a, R>) -> R
        where R: Ring + Div<Output = R>, for<'x> &'x R: RingOps<R> {
            let p = dot(base, other);
            let q = dot(base, base);
            p / q
        }
    
        fn dot<'a, R>(lhs: &ArrayView1<'a, R>, rhs: &ArrayView1<'a, R>) -> R
        where R: Ring, for<'x> &'x R: RingOps<R> {
            assert_eq!(lhs.dim(), rhs.dim());
            ndarray::Zip::from(lhs).and(rhs).fold(R::zero(), |mut acc, a, b| { 
                acc += a * b;
                acc
            })
        }
    }
}