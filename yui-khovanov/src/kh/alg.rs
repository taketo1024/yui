use std::fmt::Display;
use num_traits::Zero;
use yui::lc::{Gen, Lc};
use yui::{Elem, Ring, RingOps};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub enum KhAlgGen { 
    #[default] 
    I, 
    X
}

impl KhAlgGen { 
    #[allow(non_snake_case)]
    pub fn is_X(&self) -> bool { 
        self == &KhAlgGen::X
    }

    pub fn is_1(&self) -> bool { 
        self == &KhAlgGen::I
    }

    pub fn deg(&self) -> isize {
        match self { 
            KhAlgGen::I => 0,
            KhAlgGen::X => -2
        }
    }
}

impl Display for KhAlgGen {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self { 
            KhAlgGen::I => f.write_str("1"),
            KhAlgGen::X => f.write_str("X")
        }
    }
}

impl Elem for KhAlgGen {
    fn math_symbol() -> String {
        format!("A")
    }
}

impl Gen for KhAlgGen {}

#[derive(Clone)]
pub struct KhAlgStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    h: R,
    t: R
}

impl<R> KhAlgStr<R>
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

    pub fn prod(&self, x: KhAlgGen, y: KhAlgGen) -> Lc<KhAlgGen, R> {
        use KhAlgGen::{I, X};
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

    pub fn coprod(&self, x: KhAlgGen) -> Vec<(KhAlgGen, KhAlgGen, R)> {
        use KhAlgGen::{I, X};
        let (h, t) = (self.h(), self.t());

        match x { 
            I => if h.is_zero() { 
                vec![(X, I, R::one()), (I, X, R::one())]
            } else {
                vec![(X, I, R::one()), (I, X, R::one()), (I, I, -h.clone())]
            },
            X => if t.is_zero() { 
                vec![(X, X, R::one())]
            } else { 
                vec![(X, X, R::one()), (I, I, t.clone())]
            }
        }
    }
}

#[cfg(test)]
pub mod tests {
    use num_traits::Zero;
    use yui::lc::Lc;

    use super::{KhAlgGen, KhAlgStr};
 
    #[test]
    fn alg_gen() { 
        use KhAlgGen::{I, X};
        assert_eq!(I.deg(), 0);
        assert_eq!(X.deg(), -2);
    }

    #[test]
    fn str_prod_kh() { 
        use KhAlgGen::{I, X};
        let a = KhAlgStr::new(&0, &0);
        assert_eq!(a.prod(I, I), Lc::from((I, 1)));
        assert_eq!(a.prod(X, I), Lc::from((X, 1)));
        assert_eq!(a.prod(I, X), Lc::from((X, 1)));
        assert_eq!(a.prod(X, X), Lc::zero());
    }

    #[test]
    fn str_coprod_kh() { 
        use KhAlgGen::{I, X};
        let a = KhAlgStr::new(&0, &0);
        assert_eq!(a.coprod(I), vec![(X, I, 1), (I, X, 1)]);
        assert_eq!(a.coprod(X), vec![(X, X, 1)]);
    }
    #[test]
    fn str_prod_bn() { 
        use KhAlgGen::{I, X};
        let a = KhAlgStr::new(&1, &0);
        assert_eq!(a.prod(I, I), Lc::from((I, 1)));
        assert_eq!(a.prod(X, I), Lc::from((X, 1)));
        assert_eq!(a.prod(I, X), Lc::from((X, 1)));
        assert_eq!(a.prod(X, X), Lc::from((X, 1)));
    }

    #[test]
    fn str_coprod_bn() { 
        use KhAlgGen::{I, X};
        let a = KhAlgStr::new(&1, &0);
        assert_eq!(a.coprod(I), vec![(X, I, 1), (I, X, 1), (I, I, -1)]);
        assert_eq!(a.coprod(X), vec![(X, X, 1)]);
    }
    #[test]
    fn str_prod_lee() { 
        use KhAlgGen::{I, X};
        let a = KhAlgStr::new(&0, &1);
        assert_eq!(a.prod(I, I), Lc::from((I, 1)));
        assert_eq!(a.prod(X, I), Lc::from((X, 1)));
        assert_eq!(a.prod(I, X), Lc::from((X, 1)));
        assert_eq!(a.prod(X, X), Lc::from((I, 1)));
    }

    #[test]
    fn str_coprod_lee() { 
        use KhAlgGen::{I, X};
        let a = KhAlgStr::new(&0, &1);
        assert_eq!(a.coprod(I), vec![(X, I, 1), (I, X, 1)]);
        assert_eq!(a.coprod(X), vec![(X, X, 1), (I, I, 1)]);
    }
}