// "Extended GCD and Hermite Normal Form Algorithms via Lattice Basis Reduction",
// George Havas, Bohdan S. Majewski, and Keith R. Matthews
// https://projecteuclid.org/journals/experimental-mathematics/volume-7/issue-2/Extended-GCD-and-Hermite-normal-form-algorithms-via-lattice-basis/em/1048515660.full

use std::ops::Mul;

use ndarray::{ArrayBase, Array1, ArrayView1, array};
use crate::math::{ext::int_ext::{Integer, IntOps}, traits::{Ring, RingOps}};
use super::DnsMat;

#[derive(Debug)]
pub struct LLLCalc<R>
where R: Integer, for<'x> &'x R: IntOps<R> {
    target: DnsMat<R>,
    d: Vec<R>,     // d[0] = 1, d[i] = |b^*_i|^2 * d[i-1]. 
    l: DnsMat<R>,  // l[i,j] = d[j] * m[i, j], m[i, j] = (b_i.b^*_j)/(b^*_j.b^*_j) (0 <= j < i)
    alpha: (R, R),
    step: usize
}

type Row = usize;
type Col = usize;

impl<R> LLLCalc<R>
where R: Integer, for<'x> &'x R: IntOps<R> {
    pub fn new(target: DnsMat<R>) -> Self { 
        let m = target.nrows();
        let d = vec![R::one(); m];
        let l = DnsMat::zero((m, m)); // lower-triangular
        let alpha = (R::from(3), R::from(4)); // α = 3/4
        LLLCalc { target, d, l, alpha, step: 1 }
    }

    pub fn perform_lll(&mut self) { 
        let m = self.rows();

        while self.step < m { 

        }
    }

    fn lll_setup(&mut self) { 
        let (m, n) = self.target.shape();

        let b = self.target.array();
        let mut c = b.clone();

        self.d[0] = dot(c.row(0), c.row(0));

        for i in 1..m { 
            for j in 0..i { 
                let l0 = dot(b.row(i), c.row(j));
                
                let one = R::one();
                let d0 = if j > 0 { &self.d[j - 1] } else { &one };
                let d1 = &self.d[j];

                let c_i = Array1::from_iter( (0..n).map(|k| { 
                    &(d1 * &c[[i, k]] - &l0 * &c[[j, k]]) / d0
                }));
                c.row_mut(i).assign(&c_i);

                self.l[[i, j]] = l0;
            }

            let d0 = &self.d[i];
            self.d[i] = &dot(c.row(i), c.row(i)) / d0;
        }
    }

    pub fn process(&mut self) { 
        assert!(self.step > 0);
        let m = self.rows();

        while self.step < m { 
            println!("step: {}", self.step);

            let k = self.step;

            self.reduce(k - 1, k);

            if self.should_swap(k- 1, k) { 
                self.swap(k - 1, k);
                self.back();
            } else { 
                if k >= 2 { 
                    for i in (0..(k - 2)).rev() { 
                        self.reduce(i, k);
                    }
                }
                self.next();
            }
        }

        self.finalize();
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

    fn should_swap(&self, i: usize, k: usize) -> bool {
        assert_eq!(i, k - 1);

        let (j1, j2) = (self.nz_col_in(i), self.nz_col_in(k));
        match (j1, j2) { 
            (Some(j1), Some(j2)) => j1 <= j2,
            (Some(_),  None)     => true,
            (None,     Some(_))  => false,
            (None,     None)     => !self.lovasz_cond(i, k)
        }
    }

    fn lovasz_cond(&self, i: usize, k: usize) -> bool { 
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

    fn finalize(&mut self) { 

    }

    fn swap(&mut self, i: Row, k: Row) { 
        assert_eq!(i, k - 1);

        let m = self.rows();

        self.target.swap_rows(i, k);

        for j in 0..k-1 { 
            let slice = ndarray::s![.., j];
            let l_j = &mut self.l.array_mut().slice_mut(slice);
            l_j.swap(i, k);
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

    fn neg_row(&mut self, k: Row) { 
        let n_one = -R::one();
        self.target.mul_row(k, &n_one);
        self.l.mul_row(k, &n_one);
        self.l.mul_col(k, &n_one);
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

    fn print_current(&self) {
        dbg!(self.step);
        dbg!(&self.target);
        dbg!(&self.d);
        dbg!(&self.l);
        dbg!();
    }
}

fn dot<'a, R>(lhs: ArrayView1<'a, R>, rhs: ArrayView1<'a, R>) -> R
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert_eq!(lhs.dim(), rhs.dim());
    ndarray::Zip::from(lhs).and(rhs).fold(R::zero(), |mut acc, a, b| { 
        acc.mul_acc(a, b);
        acc
    })
}

#[cfg(test)]
mod tests {
    use ndarray::array;
    use crate::math::matrix::DnsMat;
    use super::*;
 
    #[test]
    fn swap() { 
        let a = DnsMat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut c = LLLCalc::new(a);
        c.swap(1, 2);

        assert_eq!(&c.target, &DnsMat::from(array![
            [1,-1, 3],
            [1, 2, 6],
            [1, 0, 5]
        ]));

        dbg!(&c.l);
    }

    #[test]
    fn lll_setup() { 
        let a = DnsMat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut c = LLLCalc::new(a);
        c.lll_setup();

        dbg!(c);
    }


    #[test]
    fn gcd() { 
        let a = DnsMat::from(array![
            [4],
            [6],
            [9]
        ]);
        let mut c = LLLCalc::new(a);
        c.process();

        dbg!(c.target);
    }

    #[test]
    fn test1() { 
        let a = DnsMat::from(array![
            [1,-1, 3],
            [1, 0, 5],
            [1, 2, 6]
        ]);
        let mut c = LLLCalc::new(a);
        c.process();

        dbg!(c.target);
    }
}