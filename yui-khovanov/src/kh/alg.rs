use num_traits::Zero;
use yui::lc::Lc;
use yui::{CloneAnd, Ring, RingOps};

use crate::kh::r#gen::KhGen;
use crate::kh::KhTensor;

#[derive(Clone)]
pub struct KhAlg<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    h: R,
    t: R
}

impl<R> KhAlg<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(h: &R, t: &R) -> Self { 
        Self { 
            h: h.clone(), 
            t: t.clone() 
        }
    }

    pub fn h(&self) -> &R { 
        &self.h
    }

    pub fn t(&self) -> &R { 
        &self.t
    }

    pub fn ht(&self) -> (&R, &R) { 
        (&self.h, &self.t)
    }

    pub fn mul(&self, x: KhGen, y: KhGen) -> Lc<KhGen, R> {
        use KhGen::{I, X};
        let (h, t) = (self.h(), self.t());

        match (x, y) { 
            (I, I) => Lc::from(I),
            (X, I) | (I, X) => Lc::from(X),
            (X, X) => match (h.is_zero(), t.is_zero()) { 
                (true,  true ) => Lc::zero(),
                (false, true ) => Lc::from((X, h.clone())),
                (true,  false) => Lc::from((I, t.clone())),
                (false, false) => Lc::from_iter([(X, h.clone()), (I, t.clone())])
            } 
        }
    }

    pub fn mul_tensor(&self, x: &KhTensor, in_index: (usize, usize), out_index: usize) -> Lc<KhTensor, R> { 
        assert_ne!(in_index.0, in_index.1);

        let (i, j) = if in_index.0 < in_index.1 {
            in_index
        } else { 
            (in_index.1, in_index.0)
        };
        let k = out_index;

        self.mul(x[i], x[j]).into_map_gens(|a| { 
            x.clone_and(|y| {
                y.remove(j);
                y.remove(i);
                y.insert(k, a);
            })
        })
    }

    pub fn comul(&self, x: KhGen) -> Lc<KhTensor, R> {
        use KhGen::{I, X};
        let (h, t) = (self.h(), self.t());
        let tsr = |x, y| KhTensor::from_iter([x, y]);

        match x { 
            I => if h.is_zero() { 
                Lc::from_iter([
                    (tsr(X, I), R::one()), 
                    (tsr(I, X), R::one())
                ])
            } else {
                Lc::from_iter([
                    (tsr(X, I), R::one()), 
                    (tsr(I, X), R::one()),
                    (tsr(I, I), -h)
                ])
            },
            X => if t.is_zero() { 
                Lc::from(
                    (tsr(X, X), R::one()) 
                )
            } else { 
                Lc::from_iter([
                    (tsr(X, X), R::one()), 
                    (tsr(I, I), t.clone())
                ])
            }
        }
    }

    pub fn comul_tensor(&self, x: &KhTensor, in_index: usize, out_index: (usize, usize)) -> Lc<KhTensor, R> { 
        assert_ne!(out_index.0, out_index.1);

        let i = in_index;
        let (j, k) = if out_index.0 < out_index.1 {
            out_index
        } else { 
            (out_index.1, out_index.0)
        };

        self.comul(x[i]).into_map_gens(|a| 
            x.clone_and(|y| { 
                y.remove(i);
                y.insert(j, a[0]);
                y.insert(k, a[1]);
            })
        )
    }

    pub fn sigma(&self, x: &KhGen) -> Lc<KhGen, R> {
        match x { 
            KhGen::I => Lc::from((KhGen::I, R::one())),
            KhGen::X => Lc::from_iter([
                (KhGen::X, -R::one()), 
                (KhGen::I, self.h().clone())
            ])
        }
    }
}

#[cfg(test)]
pub mod tests {
    use num_traits::Zero;
    use yui::lc::Lc;

    use crate::kh::KhTensor;

    use super::{KhGen, KhAlg};
 
    #[test]
    fn alg_gen() { 
        use KhGen::{I, X};
        assert_eq!(I.deg(), 0);
        assert_eq!(X.deg(), -2);
    }

    #[test]
    fn str_prod_kh() { 
        use KhGen::{I, X};
        let a = KhAlg::new(&0, &0);
        assert_eq!(a.mul(I, I), Lc::from((I, 1)));
        assert_eq!(a.mul(X, I), Lc::from((X, 1)));
        assert_eq!(a.mul(I, X), Lc::from((X, 1)));
        assert_eq!(a.mul(X, X), Lc::zero());
    }

    #[test]
    fn str_coprod_kh() { 
        use KhGen::{I, X};
        let a = KhAlg::new(&0, &0);
        assert_eq!(a.comul(I), Lc::from_iter([
            (KhTensor::from([X, I]), 1),
            (KhTensor::from([I, X]), 1),
        ]));
        assert_eq!(a.comul(X), Lc::from(
            (KhTensor::from_iter([X, X]), 1)
        ));
    }
    #[test]
    fn str_prod_bn() { 
        use KhGen::{I, X};
        let a = KhAlg::new(&1, &0);
        assert_eq!(a.mul(I, I), Lc::from((I, 1)));
        assert_eq!(a.mul(X, I), Lc::from((X, 1)));
        assert_eq!(a.mul(I, X), Lc::from((X, 1)));
        assert_eq!(a.mul(X, X), Lc::from((X, 1)));
    }

    #[test]
    fn str_coprod_bn() { 
        use KhGen::{I, X};
        let a = KhAlg::new(&1, &0);
        assert_eq!(a.comul(I), Lc::from_iter([
            (KhTensor::from([X, I]), 1),
            (KhTensor::from([I, X]), 1),
            (KhTensor::from([I, I]), -1),
        ]));
        assert_eq!(a.comul(X), Lc::from(
            (KhTensor::from_iter([X, X]), 1)
        ));
    }
    #[test]
    fn str_prod_lee() { 
        use KhGen::{I, X};
        let a = KhAlg::new(&0, &1);
        assert_eq!(a.mul(I, I), Lc::from((I, 1)));
        assert_eq!(a.mul(X, I), Lc::from((X, 1)));
        assert_eq!(a.mul(I, X), Lc::from((X, 1)));
        assert_eq!(a.mul(X, X), Lc::from((I, 1)));
    }

    #[test]
    fn str_coprod_lee() { 
        use KhGen::{I, X};
        let a = KhAlg::new(&0, &1);
        assert_eq!(a.comul(I), Lc::from_iter([
            (KhTensor::from([X, I]), 1),
            (KhTensor::from([I, X]), 1),
        ]));
        assert_eq!(a.comul(X), Lc::from_iter([
            (KhTensor::from([X, X]), 1),
            (KhTensor::from([I, I]), 1),
        ]));
    }

    #[test]
    fn mul_x_at() { 
        use KhGen::{I, X};
        let a = KhAlg::new(&2, &1);
        let x = KhTensor::from_iter([I, X, I]);

        assert_eq!(
            x.apply_at(0, |&x| a.mul(x, X)), 
            Lc::from(KhTensor::from_iter([X, X, I]))
        );
        assert_eq!(
            x.apply_at(1, |&x| a.mul(x, X)), 
            Lc::from_iter([
                (KhTensor::from([I, X, I]), 2),
                (KhTensor::from([I, I, I]), 1),
            ]
        ));
        assert_eq!(
            x.apply_at(2, |&x| a.mul(x, X)), 
            Lc::from(KhTensor::from_iter([I, X, X]))
        );
    }
}