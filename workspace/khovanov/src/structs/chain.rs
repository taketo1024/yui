use std::fmt::Display;
use std::ops::{Mul, MulAssign, Index};
use itertools::join;
use auto_impl_ops::auto_ops;
use yui_core::Elem;
use yui_lin_comb::{FreeGen, LinComb};
use yui_link::State;
use yui_utils::bitseq::BitSeq;

use crate::KhAlgGen;

#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct KhLabel(BitSeq);

impl KhLabel { 
    pub fn empty() -> Self { 
        Self(BitSeq::empty())
    }

    pub fn is_empty(&self) -> bool { 
        self.0.is_empty()
    }

    pub fn len(&self) -> usize { 
        self.0.len()
    }

    pub fn iter(&self) -> impl Iterator<Item = KhAlgGen> {
        self.0.iter().map(|b| 
            if b.is_zero() { 
                KhAlgGen::X
            } else { 
                KhAlgGen::I
            }
        )
    }

    pub fn push(&mut self, x: KhAlgGen) {
        if x.is_X() { 
            self.0.push_0()
        } else { 
            self.0.push_1()
        }
    }

    pub fn append(&mut self, other: KhLabel) { 
        self.0.append(other.0)
    }

    pub fn remove(&mut self, i: usize) { 
        self.0.remove(i)
    }

    pub fn insert(&mut self, i: usize, x: KhAlgGen) { 
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

    pub fn generate(len: usize) -> Vec<Self> { 
        BitSeq::generate(len).into_iter().map(|b| 
            Self(b)
        ).collect()
    }
}

impl From<KhAlgGen> for KhLabel {
    fn from(x: KhAlgGen) -> Self {
        let v = if x.is_X() { 0 } else { 1 };
        Self(BitSeq::from(v))
    }
}

impl FromIterator<KhAlgGen> for KhLabel {
    fn from_iter<I: IntoIterator<Item = KhAlgGen>>(iter: I) -> Self {
        Self(BitSeq::from_iter(iter.into_iter().map(|x| 
            if x.is_X() { 0 } else { 1 }
        )))
    }
}

impl Index<usize> for KhLabel {
    type Output = KhAlgGen;

    fn index(&self, index: usize) -> &Self::Output {
        if self.0[index].is_zero() { 
            &KhAlgGen::X
        } else { 
            &KhAlgGen::I
        }
    }
}

impl Display for KhLabel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        join(self.iter(), "").fmt(f)
    }
}

#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct KhEnhState { 
    pub state: State,
    pub label: KhLabel
}

impl KhEnhState {
    pub fn new(state: State, label: KhLabel) -> KhEnhState { 
        KhEnhState { state, label }
    }

    pub fn init() -> Self { 
        Self::new(State::empty(), KhLabel::empty())
    }

    pub fn q_deg(&self) -> isize { 
        let q = self.label.iter().map(|x| x.q_deg()).sum::<isize>();
        let r = self.label.len() as isize;
        let s = self.state.weight() as isize;
        q + r + s
    }

    pub fn append(&mut self, other: KhEnhState) { 
        let KhEnhState { state, label } = other;
        self.state.append(state);
        self.label.append(label);
    }

    pub fn is_sub(&self, other: &Self) -> bool { 
        self.state.is_sub(&other.state) && 
        self.label.is_sub(&other.label)
    }
}

#[auto_ops]
impl MulAssign<&KhEnhState> for KhEnhState {
    fn mul_assign(&mut self, rhs: &KhEnhState) {
        self.append(rhs.clone())
    }
}

impl Elem for KhEnhState { 
    fn math_symbol() -> String { 
        String::from("Kh")
    }
}

impl Display for KhEnhState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.state.is_empty() {
            write!(f, "()")
        } else if self.label.is_empty() { 
            write!(f, "({})", self.state)
        } else { 
            write!(f, "({}-{})", self.state, self.label)
        }
    }
}

impl FreeGen for KhEnhState {}

pub type KhChain<R> = LinComb<KhEnhState, R>;