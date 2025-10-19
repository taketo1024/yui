use std::fmt::Display;
use std::ops::{Add, AddAssign, Index};
use itertools::{join, Itertools};
use auto_impl_ops::auto_ops;
use yui::util::format::subscript;
use yui::{Elem, Ring, RingOps};
use yui::bitseq::BitSeq;
use yui::lc::{Gen, Lc};
use yui_link::State;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub enum KhGen { 
    #[default] 
    I, 
    X
}

impl KhGen { 
    #[allow(non_snake_case)]
    pub fn is_X(&self) -> bool { 
        self == &KhGen::X
    }

    pub fn is_1(&self) -> bool { 
        self == &KhGen::I
    }

    pub fn deg(&self) -> isize {
        match self { 
            KhGen::I => 0,
            KhGen::X => -2
        }
    }
}

impl Display for KhGen {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self { 
            KhGen::I => f.write_str("1"),
            KhGen::X => f.write_str("X")
        }
    }
}

impl Elem for KhGen {
    fn math_symbol() -> String {
        format!("A")
    }
}

impl Gen for KhGen {}

#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct KhTensor(
    pub(crate) BitSeq
);

impl KhTensor { 
    pub fn empty() -> Self { 
        Self(BitSeq::empty())
    }

    pub fn is_empty(&self) -> bool { 
        self.0.is_empty()
    }

    pub fn len(&self) -> usize { 
        self.0.len()
    }

    pub fn iter(&self) -> impl Iterator<Item = KhGen> {
        self.0.iter().map(|b| 
            if b.is_zero() { 
                KhGen::X
            } else { 
                KhGen::I
            }
        )
    }

    pub fn push(&mut self, x: KhGen) {
        if x.is_X() { 
            self.0.push_0()
        } else { 
            self.0.push_1()
        }
    }

    pub fn append(&mut self, other: KhTensor) { 
        self.0.append(other.0)
    }

    pub fn remove(&mut self, i: usize) { 
        self.0.remove(i)
    }

    pub fn insert(&mut self, i: usize, x: KhGen) { 
        if x.is_X() { 
            self.0.insert_0(i)
        } else { 
            self.0.insert_1(i)
        }
    }

    pub fn sub(&self, l: usize) -> Self { 
        Self(self.0.sub(l))
    }

    pub fn is_sub(&self, other: &Self) -> bool { 
        self.0.is_sub(&other.0)
    }

    pub fn generate(len: usize) -> impl Iterator<Item = Self> { 
        BitSeq::generate(len).map(KhTensor)
    }
}

impl From<KhGen> for KhTensor {
    fn from(x: KhGen) -> Self {
        let v = if x.is_X() { 0 } else { 1 };
        Self(BitSeq::from(v))
    }
}

impl<const N: usize> From<[KhGen; N]> for KhTensor {
    fn from(xs: [KhGen; N]) -> Self {
        Self::from_iter(xs)
    }
}

impl FromIterator<KhGen> for KhTensor {
    fn from_iter<I: IntoIterator<Item = KhGen>>(iter: I) -> Self {
        Self(BitSeq::from_iter(iter.into_iter().map(|x| 
            if x.is_X() { 0 } else { 1 }
        )))
    }
}

impl Index<usize> for KhTensor {
    type Output = KhGen;

    fn index(&self, index: usize) -> &Self::Output {
        if self.0[index].is_zero() { 
            &KhGen::X
        } else { 
            &KhGen::I
        }
    }
}

impl Display for KhTensor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        join(self.iter(), "").fmt(f)
    }
}

impl Elem for KhTensor {
    fn math_symbol() -> String {
        String::from("KhT")
    }
}

impl Gen for KhTensor {}

#[auto_ops]
impl AddAssign<KhTensor> for KhTensor {
    fn add_assign(&mut self, rhs: Self) {
        self.append(rhs);
    }
}

#[auto_ops]
impl AddAssign<KhGen> for KhTensor {
    fn add_assign(&mut self, x: KhGen) {
        self.push(x);
    }
}

#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct KhChainGen { 
    pub state: State,
    pub tensor: KhTensor,
    pub deg_shift: (isize, isize)
}

impl KhChainGen {
    pub fn new(state: State, tensor: KhTensor, deg_shift: (isize, isize)) -> KhChainGen { 
        KhChainGen { state, tensor, deg_shift }
    }

    pub fn init() -> Self { 
        KhChainGen::new(State::empty(), KhTensor::empty(), (0, 0))
    }

    pub fn h_deg(&self) -> isize { 
        let h0 = self.deg_shift.0;
        let s = self.state.weight() as isize;
        h0 + s
    }

    pub fn q_deg(&self) -> isize {
        let q0 = self.deg_shift.1;
        let d = self.tensor.iter().map(|x| x.deg()).sum::<isize>();
        let r = self.tensor.len() as isize;
        let s = self.state.weight() as isize;
        q0 + d + r + s
    }
}

impl Elem for KhChainGen { 
    fn math_symbol() -> String { 
        String::from("Kh")
    }
}

impl Display for KhChainGen {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}){}", self.tensor, self.state.iter().map(|i| subscript(i as u8)).join("") )
    }
}

impl Gen for KhChainGen {}

pub type KhChain<R> = Lc<KhChainGen, R>;
pub trait KhChainExt { 
    fn h_deg(&self) -> isize;
    fn q_deg(&self) -> isize;
}

impl<R> KhChainExt for KhChain<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn h_deg(&self) -> isize {
        self.gens().map(|x| x.h_deg()).min().unwrap_or(0)
    }
    
    fn q_deg(&self) -> isize {
        self.gens().map(|x| x.q_deg()).min().unwrap_or(0)
    }
}