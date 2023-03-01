use std::fmt::Display;
use std::ops::{Mul, MulAssign};
use itertools::join;
use auto_impl_ops::auto_ops;
use crate::links::links::State;
use crate::math::types::lin_comb::{FreeGen, LinComb};
use yui_core::Elem;
use crate::utils::format::subscript;
use yui_core::{Ring, RingOps};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum KhAlgGen { 
    I, X
}

impl KhAlgGen { 
    pub fn q_deg(&self) -> isize {
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

#[derive(Clone)]
pub struct KhAlgStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    h: R,
    t: R
}

impl<R> KhAlgStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(h:R, t:R) -> Self { 
        Self { h, t }
    }

    pub fn h(&self) -> &R { 
        &self.h
    }

    pub fn t(&self) -> &R { 
        &self.t
    }

    pub fn prod(&self, x: KhAlgGen, y: KhAlgGen) -> Vec<(KhAlgGen, R)> {
        use KhAlgGen::{I, X};
        let (h, t) = (self.h(), self.t());

        let res = match (x, y) { 
            (I, I) => vec![(I, R::one())],
            (X, I) | (I, X) => vec![(X, R::one())],
            (X, X) => vec![(X, h.clone()), (I, t.clone())]
        };

        res.into_iter().filter(|(_, a)| !a.is_zero()).collect()
    }

    pub fn coprod(&self, x: KhAlgGen) -> Vec<(KhAlgGen, KhAlgGen, R)> {
        use KhAlgGen::{I, X};
        let (h, t) = (self.h(), self.t());

        let res = match x { 
            I => vec![(X, I, R::one()), (I, X, R::one()), (I, I, -h.clone())],
            X => vec![(X, X, R::one()), (I, I, t.clone())]
        };

        res.into_iter().filter(|(_, _, a)| !a.is_zero()).collect()
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct KhEnhState { 
    state: State,
    label: Vec<KhAlgGen>
}

impl KhEnhState {
    pub fn new(state: State, label: Vec<KhAlgGen>) -> KhEnhState { 
        KhEnhState { state, label }
    }

    pub fn state(&self) -> &State {
        &self.state
    }

    pub fn label(&self) -> &Vec<KhAlgGen> {
        &self.label
    }

    pub fn label_at(&self, i: usize) -> KhAlgGen {
        self.label[i]
    }

    pub fn q_deg(&self) -> isize { 
        let q = self.label.iter().map(|x| x.q_deg()).sum::<isize>();
        let r = self.label.len() as isize;
        let s = self.state.weight() as isize;
        q + r + s
    }

    pub fn append(&mut self, other: KhEnhState) { 
        let KhEnhState { state, mut label } = other;
        self.state.append(state);
        self.label.append(&mut label);
    }
}

#[auto_ops]
impl MulAssign<&KhEnhState> for KhEnhState {
    fn mul_assign(&mut self, rhs: &KhEnhState) {
        self.append(rhs.clone())
    }
}

impl Elem for KhEnhState { 
    fn set_symbol() -> String { 
        String::from("Kh")
    }
}

impl Display for KhEnhState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let label = join(self.label.iter(), "");
        let state = join(self.state.values().map(|i| subscript(i.as_u8() as isize)), "");
        write!(f, "{}{}", label, state)
    }
}

impl FreeGen for KhEnhState {}

pub type KhChain<R> = LinComb<KhEnhState, R>;

#[cfg(test)]
pub mod tests {
    use super::{KhAlgGen, KhAlgStr};
 
    #[test]
    fn alg_gen() { 
        use KhAlgGen::{I, X};
        assert_eq!(I.q_deg(), 0);
        assert_eq!(X.q_deg(), -2);
    }

    #[test]
    fn str_prod_kh() { 
        use KhAlgGen::{I, X};
        let a = KhAlgStr::new(0, 0);
        assert_eq!(a.prod(I, I), vec![(I, 1)]);
        assert_eq!(a.prod(X, I), vec![(X, 1)]);
        assert_eq!(a.prod(I, X), vec![(X, 1)]);
        assert_eq!(a.prod(X, X), vec![]);
    }

    #[test]
    fn str_coprod_kh() { 
        use KhAlgGen::{I, X};
        let a = KhAlgStr::new(0, 0);
        assert_eq!(a.coprod(I), vec![(X, I, 1), (I, X, 1)]);
        assert_eq!(a.coprod(X), vec![(X, X, 1)]);
    }
    #[test]
    fn str_prod_bn() { 
        use KhAlgGen::{I, X};
        let a = KhAlgStr::new(1, 0);
        assert_eq!(a.prod(I, I), vec![(I, 1)]);
        assert_eq!(a.prod(X, I), vec![(X, 1)]);
        assert_eq!(a.prod(I, X), vec![(X, 1)]);
        assert_eq!(a.prod(X, X), vec![(X, 1)]);
    }

    #[test]
    fn str_coprod_bn() { 
        use KhAlgGen::{I, X};
        let a = KhAlgStr::new(1, 0);
        assert_eq!(a.coprod(I), vec![(X, I, 1), (I, X, 1), (I, I, -1)]);
        assert_eq!(a.coprod(X), vec![(X, X, 1)]);
    }
    #[test]
    fn str_prod_lee() { 
        use KhAlgGen::{I, X};
        let a = KhAlgStr::new(0, 1);
        assert_eq!(a.prod(I, I), vec![(I, 1)]);
        assert_eq!(a.prod(X, I), vec![(X, 1)]);
        assert_eq!(a.prod(I, X), vec![(X, 1)]);
        assert_eq!(a.prod(X, X), vec![(I, 1)]);
    }

    #[test]
    fn str_coprod_lee() { 
        use KhAlgGen::{I, X};
        let a = KhAlgStr::new(0, 1);
        assert_eq!(a.coprod(I), vec![(X, I, 1), (I, X, 1)]);
        assert_eq!(a.coprod(X), vec![(X, X, 1), (I, I, 1)]);
    }
}