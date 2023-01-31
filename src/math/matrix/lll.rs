// "Extended GCD and Hermite Normal Form Algorithms via Lattice Basis Reduction",
// George Havas, Bohdan S. Majewski, and Keith R. Matthews
// https://projecteuclid.org/journals/experimental-mathematics/volume-7/issue-2/Extended-GCD-and-Hermite-normal-form-algorithms-via-lattice-basis/em/1048515660.full

use std::ops::{Mul, Add, Div};

use ndarray::{ArrayBase, Array1, ArrayView1, array, Array2, ArrayView2};
use num_rational::Ratio;
use num_traits::Signed;
use crate::math::ext::ratio_ext::*;
use crate::math::{ext::int_ext::{Integer, IntOps}, traits::{Ring, RingOps, EucRing, EucRingOps}};
use super::DnsMat;

type Row = usize;
type Col = usize;

#[derive(Debug)]
pub struct LLLCalc<R>
where R: Integer, for<'x> &'x R: IntOps<R> {
    data: LLLData<R>
}

impl<R> LLLCalc<R>
where R: Integer, for<'x> &'x R: IntOps<R> {
    pub fn new(target: DnsMat<R>, alpha: (R, R)) -> Self {
        let mut data = LLLData::new(target, alpha);
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
where R: Integer, for<'x> &'x R: IntOps<R> {
    data: LLLData<R>
}

impl<R> LLLHNFCalc<R>
where R: Integer, for<'x> &'x R: IntOps<R> {
    pub fn new(target: DnsMat<R>) -> Self { 
        let alpha = (R::from(3), R::from(4)); // α = 3/4
        let data = LLLData::new(target, alpha);
        
        LLLHNFCalc { data }
    }

    pub fn process(&mut self) { 
        assert!(self.data.step > 0);
        let m = self.data.rows();

        while self.data.step < m { 
            println!("step: {}", self.data.step);

            let k = self.data.step;

            self.data.reduce(k - 1, k);

            if self.should_swap(k) { 
                self.data.swap(k);
                self.data.back();
            } else { 
                if k >= 2 { 
                    for i in (0..(k - 2)).rev() { 
                        self.data.reduce(i, k);
                    }
                }
                self.data.next();
            }
        }

        self.finalize();
    }

    fn should_swap(&self, k: usize) -> bool {
        assert!(k > 0);

        let i = k - 1;
        let (j1, j2) = (self.data.nz_col_in(i), self.data.nz_col_in(k));

        match (j1, j2) { 
            (Some(j1), Some(j2)) => j1 <= j2,
            (Some(_),  None)     => true,
            (None,     Some(_))  => false,
            (None,     None)     => !self.data.lovasz_ok(k)
        }
    }

    fn finalize(&mut self) { 

    }
}

#[derive(Debug, PartialEq, Eq)]
struct LLLData<R>
where R: Integer, for<'x> &'x R: IntOps<R> {
    target: DnsMat<R>,
    det: Vec<R>,        // D[i] = det(b_1, ..., b_i)^2 = Π^i |b^*_j|^2. 
    lambda: DnsMat<R>,  // l[i,j] = D[j] * p_ij (0 <= j < i)
    alpha: (R, R),
    step: usize
}

impl<R> LLLData<R>
where R: Integer, for<'x> &'x R: IntOps<R> {
    fn new(target: DnsMat<R>, alpha: (R, R)) -> Self { 
        // 1/4 < α < 1 
        let (p, q) = alpha.clone();
        assert!(&q < &(&R::from(4) * &p));
        assert!(&p < &q);

        let m = target.nrows();
        let det = vec![R::one(); m];
        let lambda = DnsMat::zero((m, m)); // lower-triangular

        LLLData { target, det, lambda, alpha, step: 1 }
    }

    #[allow(unused)]
    fn setup(&mut self) { 
        let m = self.rows();
        self.setup_upto(m - 1);
    }

    fn setup_upto(&mut self, i0: Row) { 
        assert!(i0 < self.rows());

        let s = ndarray::s![0..=i0, ..];
        let b = self.target.array().slice(s);
        let (_, l, d) = large_orth_basis(&b);

        self.lambda = DnsMat::from(l);
        self.det = d;
    }

    fn lovasz_ok(&self, k: usize) -> bool { 
        assert!(k > 0);

        let (d, l) = (&self.det, &self.lambda);
        let (p, q) = &self.alpha;

        let one = R::one();
        let d0 = if k >= 2 { &d[k - 2] } else { &one };
        let d1 = &d[k - 1];
        let d2 = &d[k];
        let l0 = &l[[k, k - 1]];

        q * &(d0 * d2 + l0 * l0) >= p * &(d1 * d1) 

        // q * (d[k-2] * d[k] + (l[k,k-1])^2) >= p d[k-1]^2
        // ⇔ (d[k-2]/d[k-1]) * (d[k]/d[k-1]) + (l[k,k-1]/d[k-1])^2 >= p/q
        // ⇔ |b^*{k-1}|^{-2} * |b^*_k|^2 + m[k,k-1]^2 >= α
        // ⇔ |b^*_k|^2 >= (α - m[k,k-1]^2)|b^*_{k-1}|^2 : Lovasz condition
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

            let s = l1 * l0 + l2 * d0;
            let t = l1 * d2 - l2 * l0;

            self.lambda[[i, k-1]] = &s / d1;
            self.lambda[[i, k]]   = &t / d1;
        }

        // λ[k, k-1] remains unchanged.

        let l0 = &self.lambda[[k,k-1]];
        self.det[k-1] = &(d0 * d2 + l0 * l0) / d1;
    }

    fn neg_row(&mut self, i: Row) { 
        let n_one = -R::one();
        self.target.mul_row(i, &n_one);
        self.lambda.mul_row(i, &n_one);
        self.lambda.mul_col(i, &n_one);
    }

    fn add_row_to(&mut self, i: Row, k: Row, r: &R) {
        assert!(i < k);
        self.target.add_row_to(i, k, r);

        self.lambda[[k, i]] += r * &self.det[i];

        for j in 0..i { 
            let a = r * &self.lambda[[i, j]];
            self.lambda[[k, j]] += a;
        }
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
        self.print_current();
    }

    fn rows(&self) -> usize { 
        self.target.nrows()
    }

    fn cols(&self) -> usize { 
        self.target.ncols()
    }

    fn print_current(&self) {
        println!("step: {}", self.step);
        println!("target:\n{}", self.target);
    }
}

// -- helper funcs -- //

#[allow(unused)]
fn is_reduced<R>(b: &ArrayView2<R>, alpha: &(R, R)) -> bool
where R: Integer, for<'x> &'x R: IntOps<R> {
    let m = b.nrows();

    let (c, l) = orth_basis(b);
    let alpha = Ratio::from(alpha.clone());
    let thr = Ratio::new(R::one(), R::from(2));

    let size_reduced = l.iter().all(|r| &r.abs() <= &thr);
    let lovasz_ok = (1..m).all(|i| {
        let c0 = &c.row(i - 1);
        let c1 = &c.row(i);
        let m = &l[[i, i - 1]];
        is_lovasz_ok(c0, c1, m, &alpha)
    });

    dbg!(size_reduced, lovasz_ok);

    size_reduced && lovasz_ok
}

#[allow(unused)]
fn is_lovasz_ok<R>(c0: &ArrayView1<Ratio<R>>, c1: &ArrayView1<Ratio<R>>, m: &Ratio<R>, alpha: &Ratio<R>) -> bool
where R: Integer, for<'x> &'x R: IntOps<R> {
    let r0 = dot::<Ratio<R>>(c0, c0);
    let r1 = dot::<Ratio<R>>(c1, c1);
    r1 >= (alpha - m * m) * r0
}

#[allow(unused)]
fn orth_basis<R>(b: &ArrayView2<R>) -> (Array2<Ratio<R>>, Array2<Ratio<R>>)
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

fn large_orth_basis<R>(b: &ArrayView2<R>) -> (Array2<R>, Array2<R>, Vec<R>)
where R: Integer, for<'x> &'x R: IntOps<R> {
    let (m, n) = (b.nrows(), b.ncols());
    let one = R::one();

    let mut c = b.to_owned();
    let mut d = vec![R::one(); m];
    let mut l = Array2::zeros((m, m));

    d[0] = dot(&c.row(0), &c.row(0));

    for i in 1..m { 
        for j in 0..i { 
            let l0 = dot(&b.row(i), &c.row(j));
            let d0 = if j > 0 { &d[j - 1] } else { &one };
            let d1 = &d[j];
            
            let c_i = Array1::from_iter( (0..n).map(|k| { 
                &(d1 * &c[[i, k]] - &l0 * &c[[j, k]]) / d0
            }));

            l[[i, j]] = l0;
            c.row_mut(i).assign(&c_i);
        }

        let c_i = &c.row(i);
        let d0 = &d[i - 1];
        
        d[i] = &dot(&c_i, &c_i) / d0;
    }

    (c, l, d)
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

fn smul<'a, R>(r: &'a R, rhs: &ArrayView1<'a, R>) -> Array1<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    rhs.map(|x| r * x)
}

#[cfg(test)]
mod tests {
    use ndarray::array;
    use num_traits::Zero;
    use crate::math::matrix::DnsMat;
    use super::*;
 
    #[test]
    fn test_orth_basis() {
        let a = array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ];
        let (c, _) = orth_basis(&a.view());

        assert!(dot(&c.row(0), &c.row(1)).is_zero());
        assert!(dot(&c.row(0), &c.row(2)).is_zero());
        assert!(dot(&c.row(1), &c.row(2)).is_zero());
        
        let det: Ratio<i32> = (0..3).map(|i| dot(&c.row(i), &c.row(i))).product();
        assert!(det.is_integer());
        assert_eq!(det.to_integer(), 9);

        assert_eq!(dot(&c.row(0), &c.row(0)), Ratio::new(11, 1));
        assert_eq!(dot(&c.row(1), &c.row(1)), Ratio::new(30, 11));
        assert_eq!(dot(&c.row(2), &c.row(2)), Ratio::new(3, 10));
    }

    #[test]
    fn test_large_orth_basis() {
        let a = array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ];
        let (c, l, d) = large_orth_basis(&a.view());

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
        let alpha = (3, 4);
        let mut data = LLLData::new(a, alpha);

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
    fn setup_upto() { 
        let a = DnsMat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let alpha = (3, 4);
        let mut data = LLLData::new(a, alpha);

        data.setup_upto(1);

        assert_eq!(data.det.len(), 2);
        assert_eq!(data.det[0], 11);
        assert_eq!(data.det[1], 30);

        assert_eq!(data.lambda, DnsMat::from(array![
            [0,  0],
            [16, 0],
        ]));
    }

    #[test]
    fn swap() { 
        let mut a = DnsMat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let alpha = (3, 4);

        let mut data = LLLData::new(a.clone(), alpha);
        data.setup();
        
        // first swap
        data.swap(1);
        
        // compare data
        a.swap_rows(0, 1);

        let mut data2 = LLLData::new(a.clone(), alpha);
        data2.setup();

        assert_eq!(&data, &data2);

        // second swap
        data.swap(2);

        // compare data
        a.swap_rows(1, 2);

        let mut data3 = LLLData::new(a.clone(), alpha);
        data3.setup();

        assert_eq!(&data, &data3);
     }

    #[test]
    fn add_row_to() { 
        let mut a = DnsMat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let alpha = (3, 4);

        let mut data = LLLData::new(a.clone(), alpha);
        data.setup();
        
        // first addition
        data.add_row_to(0, 1, &2);
        
        // compare data
        a.add_row_to(0, 1, &2);

        let mut data2 = LLLData::new(a.clone(), alpha);
        data2.setup();

        assert_eq!(&data, &data2);

        // second addition
        data.add_row_to(1, 2, &-3);

        // compare data
        a.add_row_to(1, 2, &-3);
        
        let mut data3 = LLLData::new(a.clone(), alpha);
        data3.setup();

        assert_eq!(&data, &data3);
     }

     #[test]
     fn lll() { 
        let a = DnsMat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let alpha = (3, 4);

        let mut calc = LLLCalc::new(a, alpha);
        calc.process();
        let res = calc.result();

        assert_eq!(res, DnsMat::from(array![
            [0, 1, -1],
            [1, 0, -1],
            [1, 1, 1]
        ]));
        assert!(is_reduced(&res.array().view(), &alpha));
    }

     #[test]
     fn lll_gcdx() { 
        // MEMO: γ = 10
        let a = DnsMat::from(array![
            [1, 0, 0, 40],
            [0, 1, 0, 60],
            [0, 0, 1, 90]
        ]);
        let alpha = (3, 4);

        let mut calc = LLLCalc::new(a.clone(), alpha);
        calc.process();
        let res = calc.result();

        assert_eq!(res, DnsMat::from(array![
            [3, -2, 0, 0],
            [0, 3, -2, 0],
            [-2, 0, 1, 10]
        ]));
        assert!(is_reduced(&res.array().view(), &alpha));
      }
}