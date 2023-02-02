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

use ndarray::{Array1, ArrayView1, Array2, ArrayView2};
use log::trace;
use num_bigint::BigInt;

use crate::math::traits::{Ring, RingOps, EucRing, EucRingOps, DivRound};
use crate::math::ext::int_ext::{Integer, IntOps};
use crate::math::types::quad_int::{GaussInt, QuadInt, EisenInt};
use super::DnsMat;

pub trait LLLRingOps<T>: EucRingOps<T> {}

pub trait LLLRing: EucRing + LLLRingOps<Self> + DivRound
where for<'x> &'x Self: LLLRingOps<Self> {
    type Int: PartialOrd + Ord;
    fn alpha() -> (Self, Self);
    fn as_int(&self) -> Option<Self::Int>;
    fn conj(&self) -> Self;
}

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
                (Self::from(I::from($p)), Self::from(I::from($q)))
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

type Row = usize;
type Col = usize;

#[derive(Debug)]
pub struct LLLCalc<R>
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    data: LLLData<R>
}

impl<R> LLLCalc<R>
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    pub fn new(target: DnsMat<R>) -> Self {
        let mut data = LLLData::new(target);
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

    pub fn result(self) -> DnsMat<R> { 
        self.data.target
    }
}

#[derive(Debug)]
pub struct LLLHNFCalc<R>
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    data: LLLData<R>
}

impl<R> LLLHNFCalc<R>
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    pub fn new(target: DnsMat<R>) -> Self { 
        let data = LLLData::new(target);
        LLLHNFCalc { data }
    }

    pub fn process(&mut self) { 
        assert!(self.data.step > 0);
        let m = self.data.rows();

        while self.data.step < m { 
            self.iterate();
        }

        self.finalize();
    }

    fn iterate(&mut self) { 
        trace!("step: {}", self.data.step);

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

    pub fn result(self) -> DnsMat<R> { 
        self.data.target
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

    fn finalize(&mut self) { 

    }
}

#[derive(Debug, PartialEq, Eq)]
struct LLLData<R>
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    target: DnsMat<R>,
    det: Vec<R>,        // D[i] = det(b_1, ..., b_i)^2 = Π^i |b^*_j|^2. 
    lambda: DnsMat<R>,  // l[i,j] = D[j] * p_ij (0 <= j < i)
    step: usize
}

impl<R> LLLData<R>
where R: LLLRing, for<'x> &'x R: LLLRingOps<R> {
    fn new(target: DnsMat<R>) -> Self { 
        let m = target.nrows();
        let det = vec![R::one(); m];
        let lambda = DnsMat::zero((m, m)); // lower-triangular

        LLLData { target, det, lambda, step: 1 }
    }

    fn setup(&mut self) { 
        let b = self.target.array();
        let (_, l, d) = orthogonalize(&b.view());

        self.lambda = DnsMat::from(l);
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

        let lhs = q * &(d0 * d2 + l0 * &l0.conj());
        let rhs = p * &(d1 * d1);

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

        let m = self.lambda.ncols();

        // λ[.., k-1] <--> λ[.., k]
        for i in k+1..m { 
            let l = &self.lambda;
            let l0 = &l[[k, k-1]];
            let l1 = &l[[i, k-1]];
            let l2 = &l[[i, k]];

            let s = &l0.conj() * l1 + l2 * d0;
            let t = l1 * d2 - l2 * l0;

            self.lambda[[i, k-1]] = &s / d1;
            self.lambda[[i, k]]   = &t / d1;
        }

        // λ[k, k-1] remains unchanged.

        let l0 = &self.lambda[[k,k-1]];

        self.det[k-1] = &(d0 * d2 + l0 * &l0.conj()) / d1;
        self.lambda[[k, k-1]] = l0.conj();

        trace!("swap {},{}.\n{}", k-1, k, self.target);
    }

    fn mul_row(&mut self, i: Row, r: &R) { 
        assert!(r.is_unit());

        self.target.mul_row(i, r);
        self.lambda.mul_row(i, r);
        self.lambda.mul_col(i, &r.conj());

        trace!("mul {} to row {}.\n{}", r, i, self.target);
    }

    fn add_row_to(&mut self, i: Row, k: Row, r: &R) {
        assert!(i < k);
        self.target.add_row_to(i, k, r);

        self.lambda[[k, i]] += r * &self.det[i];

        for j in 0..i { 
            let a = r * &self.lambda[[i, j]];
            self.lambda[[k, j]] += a;
        }

        trace!("add-row {} to {}, mul {}.\n{}", i, k, r, self.target);
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
        self.target.nrows()
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
        
        d[i] = &h_dot(&c_i, &c_i) / d0;
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
        acc.mul_acc(a, &b.conj());
        acc
    })
}

#[cfg(test)]
mod tests {
    use ndarray::array;
    use crate::math::matrix::DnsMat;
    use super::*;
 
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
    fn setup() { 
        let a = DnsMat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut data = LLLData::new(a);

        data.setup();

        assert_eq!(data.det.len(), 3);
        assert_eq!(data.det[0], 11);
        assert_eq!(data.det[1], 30);
        assert_eq!(data.det[2], 9);

        assert_eq!(data.lambda, DnsMat::from(array![
            [0,  0, 0],
            [16, 0, 0],
            [17,69, 0]
        ]));
    }

    #[test]
    fn swap() { 
        let a0 = DnsMat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut data0 = LLLData::new(a0);
        data0.setup();
        data0.swap(1);
        data0.swap(2);
        
        // compare data
        let a1 = DnsMat::from(array![
            [1, 0, 5],
            [1, 2, 6],
            [1,-1, 3]
        ]);
        let mut data1 = LLLData::new(a1);
        data1.setup();

        assert_eq!(data0, data1);
     }

    #[test]
    fn add_row_to() { 
        let a0 = DnsMat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut data0 = LLLData::new(a0);
        data0.setup();
        data0.add_row_to(0, 1, &2);
        data0.add_row_to(1, 2, &-3);
        
        // compare data
        let a1 = DnsMat::from(array![
            [1,-1, 3],
            [3,-2,11],
            [-8,8,-27]
        ]);
        let mut data1 = LLLData::new(a1);
        data1.setup();

        assert_eq!(data0, data1);
     }

     #[test]
     fn mul_row() { 
        let a0 = DnsMat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut data0 = LLLData::new(a0);
        data0.setup();
        data0.mul_row(1, &-1);
        data0.mul_row(2, &-1);
        
        // compare data
        let a1 = DnsMat::from(array![
            [1,-1, 3],
            [-1, 0, -5],
            [-1, -2, -6]
        ]);
        let mut data1 = LLLData::new(a1);
        data1.setup();

        assert_eq!(data0, data1);
    }

     #[test]
     fn lll() { 
        let a = DnsMat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut calc = LLLCalc::new(a);
        calc.process();
        let res = calc.result();

        assert_eq!(res, DnsMat::from(array![
            [0, 1, -1],
            [1, 0, -1],
            [1, 1, 1]
        ]));
        assert!( helper::is_reduced( &res.array().view() ) );
    }

     #[test]
     fn lll_gcdx() { 
        // MEMO: γ = 10
        let a = DnsMat::from(array![
            [1, 0, 0, 40],
            [0, 1, 0, 60],
            [0, 0, 1, 90]
        ]);
        let mut calc = LLLCalc::new(a.clone());
        calc.process();
        let res = calc.result();

        assert_eq!(res, DnsMat::from(array![
            [3, -2, 0, 0],
            [0, 3, -2, 0],
            [-2, 0, 1, 10]
        ]));
        assert!( helper::is_reduced( &res.array().view() ));
      }

     #[test]
     fn hnf() { 
        let a: DnsMat<i64> = DnsMat::from(array![
            [8,    44,   43],
            [4,    10,   43],
            [56, -550, -328],
            [76,   10,   42]
        ]);
        let mut calc = LLLHNFCalc::new(a.clone());
        calc.process();
        let res = calc.result();

        assert_eq!(res, DnsMat::from(array![
            [0, 0, 0],
            [0, 0, 5],
            [0, 6,-2],
            [4,-2, 2]
        ]));
    }

    #[test]
    fn setup_gauss() { 
        type A = GaussInt<i64>;

        let i = |a, b| A::new(a, b);
        let a = DnsMat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(3, 3), i(-2, 4), i(6, 2)],
            [i(2, 2), i(-8, 0), i(-9, 1)],
        ]);
        let mut data = LLLData::new(a);

        data.setup();

        assert_eq!(data.det.len(), 3);
        assert_eq!(data.det[0], i(129, 0));
        assert_eq!(data.det[1], i(7436, 0));
        assert_eq!(data.det[2], i(161408, 0));

        assert_eq!(data.lambda, DnsMat::from(array![
            [i(0, 0),     i(0, 0),      i(0, 0)],
            [i(49, 15),   i(0, 0),      i(0, 0)],
            [i(-114, 48), i(1770,3162), i(0, 0)]
        ]));
    }

    #[test]
    fn swap_gauss() { 
        type A = GaussInt<i64>;
        let i = |a, b| A::new(a, b);
        let a0 = DnsMat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(3, 3), i(-2, 4), i(6, 2)],
            [i(2, 2), i(-8, 0), i(-9, 1)],
        ]);
        let mut data0 = LLLData::new(a0);
        data0.setup();
        data0.swap(1);
        data0.swap(2);
        
        // compare data
        let a1 = DnsMat::from(array![
            [i(3, 3), i(-2, 4), i(6, 2)],
            [i(2, 2), i(-8, 0), i(-9, 1)],
            [i(-2, 3), i(7, 3), i(7, 3)],
        ]);
        let mut data1 = LLLData::new(a1);
        data1.setup();

        assert_eq!(data0, data1);
    }

    #[test]
    fn add_row_to_gauss() { 
        type A = GaussInt<i64>;
        let i = |a, b| A::new(a, b);
        let a0 = DnsMat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(3, 3), i(-2, 4), i(6, 2)],
            [i(2, 2), i(-8, 0), i(-9, 1)],
        ]);
        let mut data0 = LLLData::new(a0);
        data0.setup();
        data0.add_row_to(0, 1, &i(1, 1));
        data0.add_row_to(1, 2, &i(-3, 2));
        
        // compare data
        let a1 = DnsMat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(-2, 4), i(2, 14), i(10, 12)],
            [i(0, -14), i(-42, -38), i(-63, -15)],
        ]);
        let mut data1 = LLLData::new(a1);
        data1.setup();

        assert_eq!(data0, data1);
    }
 
    #[test]
    fn mul_row_gauss() { 
        type A = GaussInt<i64>;
        let i = |a, b| A::new(a, b);
        let a0 = DnsMat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(3, 3), i(-2, 4), i(6, 2)],
            [i(2, 2), i(-8, 0), i(-9, 1)],
        ]);

        let mut data0 = LLLData::new(a0.clone());
        data0.setup();
        data0.mul_row(1, &i(0, 1));
        data0.mul_row(2, &i(0, -1));

        // compare data
        let a1 = DnsMat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(-3, 3), i(-4, -2), i(-2, 6)],
            [i(2, -2), i(0, 8), i(1, 9)],
        ]);
        let mut data1 = LLLData::new(a1);
        data1.setup();

        assert_eq!(data0, data1);
    }

    #[test]
    fn hnf_gauss() { 
        type A = GaussInt<i64>;
        let i = |a, b| A::new(a, b);

        let a: DnsMat<A> = DnsMat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(3, 3), i(-2, 4), i(6, 2)],
            [i(2, 2), i(-8, 0), i(-9, 1)],
        ]);

        let mut calc = LLLHNFCalc::new(a);
        calc.process();
        let _res = calc.result();

        // TODO
    }

    #[test]
    fn hnf_eisen() { 
        type A = EisenInt<i64>;
        let i = |a, b| A::new(a, b);

        let a: DnsMat<A> = DnsMat::from(array![
            [i(-2, 3), i(7, 3), i(7, 3)],
            [i(3, 3), i(-2, 4), i(6, 2)],
            [i(2, 2), i(-8, 0), i(-9, 1)],
        ]);

        let mut calc = LLLHNFCalc::new(a);
        calc.process();
        let _res = calc.result();

        // TODO
    }

    mod helper { 
        use std::ops::Div;
        use num_rational::Ratio;
        use num_traits::Signed;
        use super::*;
        use crate::math::ext::int_ext::{Integer, IntOps};
    
        pub fn is_reduced<R>(b: &ArrayView2<R>) -> bool
        where R: Integer + LLLRing, for<'x> &'x R: IntOps<R> + LLLRingOps<R> {
            let m = b.nrows();
    
            let (c, l) = gram_schmidt(b);
            let alpha = Ratio::from(R::alpha());
            let thr = Ratio::new(R::one(), R::from(2));
    
            let size_reduced = l.iter().all(|r| &r.abs() <= &thr);
            let lovasz_ok = (1..m).all(|i| {
                let c0 = &c.row(i - 1);
                let c1 = &c.row(i);
                let m = &l[[i, i - 1]];
                is_lovasz_ok(c0, c1, m, &alpha)
            });
    
            size_reduced && lovasz_ok
        }
    
        pub fn is_lovasz_ok<R>(c0: &ArrayView1<Ratio<R>>, c1: &ArrayView1<Ratio<R>>, m: &Ratio<R>, alpha: &Ratio<R>) -> bool
        where R: Integer, for<'x> &'x R: IntOps<R> {
            let r0 = dot::<Ratio<R>>(c0, c0);
            let r1 = dot::<Ratio<R>>(c1, c1);
            r1 >= (alpha - m * m) * r0
        }
    
        fn gram_schmidt<R>(b: &ArrayView2<R>) -> (Array2<Ratio<R>>, Array2<Ratio<R>>)
        where R: Integer, for<'x> &'x R: IntOps<R> {
            let m = b.nrows();
    
            let mut c = b.map(|x| Ratio::from(x.clone()));
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
            let p = dot(&base, &other);
            let q = dot(&base, &base);
            p / q
        }
    
        fn dot<'a, R>(lhs: &ArrayView1<'a, R>, rhs: &ArrayView1<'a, R>) -> R
        where R: Ring, for<'x> &'x R: RingOps<R> {
            assert_eq!(lhs.dim(), rhs.dim());
            ndarray::Zip::from(lhs).and(rhs).fold(R::zero(), |mut acc, a, b| { 
                acc.mul_acc(a, b);
                acc
            })
        }
    }
}