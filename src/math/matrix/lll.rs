// "Extended GCD and Hermite Normal Form Algorithms via Lattice Basis Reduction",
// George Havas, Bohdan S. Majewski, and Keith R. Matthews
// https://projecteuclid.org/journals/experimental-mathematics/volume-7/issue-2/Extended-GCD-and-Hermite-normal-form-algorithms-via-lattice-basis/em/1048515660.full

use std::ops::{Mul, Add, Div};

use ndarray::{ArrayBase, Array1, ArrayView1, array, Array2, ArrayView2};
use num_rational::Ratio;
use crate::math::ext::ratio_ext::*;
use crate::math::{ext::int_ext::{Integer, IntOps}, traits::{Ring, RingOps, EucRing, EucRingOps}};
use super::DnsMat;

#[derive(Debug)]
pub struct LLLCalc<R>
where R: Integer, for<'x> &'x R: IntOps<R> {
    data: LLLData<R>
}

type Row = usize;
type Col = usize;

impl<R> LLLCalc<R>
where R: Integer, for<'x> &'x R: IntOps<R> {
    pub fn new(target: DnsMat<R>) -> Self { 
        let alpha = (R::from(3), R::from(4)); // α = 3/4
        let data = LLLData::new(target, alpha);
        
        LLLCalc { data }
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
            (None,     None)     => !self.data.lovasz_cond(k)
        }
    }

    fn finalize(&mut self) { 

    }
}

#[derive(Debug)]
struct LLLData<R>
where R: Integer, for<'x> &'x R: IntOps<R> {
    target: DnsMat<R>,
    d: Vec<R>,     // d[0] = 1, d[i] = |b^*_i|^2 * d[i-1]. 
    l: DnsMat<R>,  // l[i,j] = d[j] * m[i, j], m[i, j] = (b_i.b^*_j)/(b^*_j.b^*_j) (0 <= j < i)
    alpha: (R, R),
    step: usize
}

impl<R> LLLData<R>
where R: Integer, for<'x> &'x R: IntOps<R> {
    fn new(target: DnsMat<R>, alpha: (R, R)) -> Self { 
        let m = target.nrows();
        let d = vec![R::one(); m];
        let l = DnsMat::zero((m, m)); // lower-triangular

        LLLData { target, d, l, alpha, step: 1 }
    }

    fn setup(&mut self) { 
        let m = self.rows();
        self.setup_upto(m - 1);
    }

    fn setup_upto(&mut self, i0: Row) { 
        assert!(i0 < self.rows());

        let s = ndarray::s![0..=i0, ..];
        let b = self.target.array().slice(s);
        let (_, l, d) = large_orth_basis(&b);

        self.l = DnsMat::from(l);
        self.d = d;
    }

    fn lovasz_cond(&self, k: usize) -> bool { 
        assert!(k > 0);

        let (d, l) = (&self.d, &self.l);
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
        let j = self.nz_col_in(i);

        if let Some(j) = j { 
            if self.target[[i, j]].is_negative() { 
                self.neg_row(i);
            }
        }

        let q = if let Some(j) = j { 
            let a = &self.target;
            let a_ij = &a[[i, j]];
            let a_kj = &a[[k, j]];
            a_kj.div_round(a_ij)
        } else {
            let (d, l) = (&self.d, &self.l);
            let l_ki = &l[[k, i]];
            let d_i = &d[i];
            l_ki.div_round(d_i)
        };

        if !q.is_zero() { 
            self.add_row_to(i, k, &-q) // a[k] -= q & a[i]
        }
    }

    fn swap(&mut self, k: Row) { 
        assert!(k > 0);
        let m = self.rows();

        self.target.swap_rows(k - 1, k);

        for j in 0..k-1 { 
            let slice = ndarray::s![.., j];
            let l_j = &mut self.l.array_mut().slice_mut(slice);
            l_j.swap(k - 1, k);
        }

        let one = R::one();
        let d = &self.d;
        let d0 = if k >= 2 { &d[k - 2] } else { &one };
        let d1 = &d[k - 1];
        let d2 = &d[k];

        for i in k+1..m { 
            let l = &self.l;
            let t = &l[[i, k-1]] * d2 - &l[[i, k]] * &l[[k, k-1]];
            let s = &(&l[[i, k-1]] * &l[[k, k-1]] + &l[[i, k]] * d0) / d1;

            self.l[[i, k-1]] = s;
            self.l[[i, k]]   = &t / d1;
        }

        let l0 = &self.l[[k,k-1]];
        self.d[k-1] = &(d0 * d2 + l0 * l0) / d1;
    }

    fn neg_row(&mut self, i: Row) { 
        let n_one = -R::one();
        self.target.mul_row(i, &n_one);
        self.l.mul_row(i, &n_one);
        self.l.mul_col(i, &n_one);
    }

    fn add_row_to(&mut self, i: Row, k: Row, r: &R) {
        self.target.add_row_to(i, k, r);

        self.l[[k, i]] += r * &self.d[i];

        if i > 0 { 
            for j in 0..i-1 { 
                let a = r * &self.l[[i, j]];
                self.l[[k, j]] += a;
            }
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
        dbg!(self.step);
        dbg!(&self.target);
        dbg!(&self.d);
        dbg!(&self.l);
        dbg!();
    }
}

// -- helper funcs -- //

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

fn sdiv<'a, R>(lhs: &ArrayView1<'a, R>, r: &'a R) -> Array1<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    lhs.map(|x| x / r)
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
        let mut data = LLLData::new(a, (3, 4));
        data.setup();

        assert_eq!(data.d.len(), 3);
        assert_eq!(data.d[0], 11);
        assert_eq!(data.d[1], 30);
        assert_eq!(data.d[2], 9);

        assert_eq!(data.l, DnsMat::from(array![
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
        let mut data = LLLData::new(a, (3, 4));
        data.setup_upto(1);

        assert_eq!(data.d.len(), 2);
        assert_eq!(data.d[0], 11);
        assert_eq!(data.d[1], 30);

        assert_eq!(data.l, DnsMat::from(array![
            [0,  0],
            [16, 0],
        ]));
    }

    #[test]
    fn swap() { 
        let a = DnsMat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut data = LLLData::new(a, (3, 4));
        data.swap(2);

        assert_eq!(&data.target, &DnsMat::from(array![
            [1,-1, 3],
            [1, 2, 6],
            [1, 0, 5]
        ]));
    }
}