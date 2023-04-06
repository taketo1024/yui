use std::fmt::Display;
use std::ops::{Mul, MulAssign, Index};
use itertools::join;
use auto_impl_ops::auto_ops;
use yui_core::Elem;
use yui_lin_comb::{FreeGen, LinComb};
use yui_link::State;
use yui_utils::bitseq::BitSeq;

use crate::KhAlgLabel;

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

    pub fn iter(&self) -> impl Iterator<Item = KhAlgLabel> {
        self.0.iter().map(|b| 
            if b.is_zero() { 
                KhAlgLabel::X
            } else { 
                KhAlgLabel::I
            }
        )
    }

    pub fn push(&mut self, x: KhAlgLabel) {
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

    pub fn insert(&mut self, i: usize, x: KhAlgLabel) { 
        if x.is_X() { 
            self.0.insert_0(i)
        } else { 
            self.0.insert_1(i)
        }
    }

    pub fn generate(len: usize) -> Vec<Self> { 
        BitSeq::generate(len).into_iter().map(|b| 
            Self(b)
        ).collect()
    }
}

impl FromIterator<KhAlgLabel> for KhLabel {
    fn from_iter<I: IntoIterator<Item = KhAlgLabel>>(iter: I) -> Self {
        Self(BitSeq::from_iter(iter.into_iter().map(|x| 
            if x.is_X() { 0 } else { 1 }
        )))
    }
}

impl Index<usize> for KhLabel {
    type Output = KhAlgLabel;

    fn index(&self, index: usize) -> &Self::Output {
        if self.0[index].is_zero() { 
            &KhAlgLabel::X
        } else { 
            &KhAlgLabel::I
        }
    }
}

impl Display for KhLabel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        join(self.iter(), "").fmt(f)
    }
}

#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct KhGen { 
    pub state: State,
    pub label: KhLabel
}

impl KhGen {
    pub fn new(state: State, label: KhLabel) -> KhGen { 
        KhGen { state, label }
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

    pub fn append(&mut self, other: KhGen) { 
        let KhGen { state, label } = other;
        self.state.append(state);
        self.label.append(label);
    }
}

#[auto_ops]
impl MulAssign<&KhGen> for KhGen {
    fn mul_assign(&mut self, rhs: &KhGen) {
        self.append(rhs.clone())
    }
}

impl Elem for KhGen { 
    fn set_symbol() -> String { 
        String::from("Kh")
    }
}

impl Display for KhGen {
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

impl FreeGen for KhGen {}

pub type KhChain<R> = LinComb<KhGen, R>;