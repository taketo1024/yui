use std::fmt::Display;
use std::ops::{Add, AddAssign, Index};
use itertools::{join, Itertools};
use auto_impl_ops::auto_ops;
use yui::util::format::subscript;
use yui::Elem;
use yui::bitseq::BitSeq;
use yui::lc::Gen;
use yui_link::State;

use crate::kh::KhAlgGen;

#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct KhLabel(
    pub(crate) BitSeq
);

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

    pub fn generate(len: usize) -> impl Iterator<Item = Self> { 
        BitSeq::generate(len).map(KhLabel)
    }
}

impl From<KhAlgGen> for KhLabel {
    fn from(x: KhAlgGen) -> Self {
        let v = if x.is_X() { 0 } else { 1 };
        Self(BitSeq::from(v))
    }
}

impl<const N: usize> From<[KhAlgGen; N]> for KhLabel {
    fn from(xs: [KhAlgGen; N]) -> Self {
        Self::from_iter(xs)
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

#[auto_ops]
impl AddAssign<KhLabel> for KhLabel {
    fn add_assign(&mut self, rhs: Self) {
        self.append(rhs);
    }
}

#[auto_ops]
impl AddAssign<KhAlgGen> for KhLabel {
    fn add_assign(&mut self, x: KhAlgGen) {
        self.push(x);
    }
}

#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct KhGen { 
    pub state: State,
    pub label: KhLabel,
    pub deg_shift: (isize, isize)
}

impl KhGen {
    pub fn new(state: State, label: KhLabel, deg_shift: (isize, isize)) -> KhGen { 
        KhGen { state, label, deg_shift }
    }

    pub fn init() -> Self { 
        KhGen::new(State::empty(), KhLabel::empty(), (0, 0))
    }

    pub fn h_deg(&self) -> isize { 
        let h0 = self.deg_shift.0;
        let s = self.state.weight() as isize;
        h0 + s
    }

    pub fn q_deg(&self) -> isize {
        let q0 = self.deg_shift.1;
        let d = self.label.iter().map(|x| x.deg()).sum::<isize>();
        let r = self.label.len() as isize;
        let s = self.state.weight() as isize;
        q0 + d + r + s
    }
}

impl Elem for KhGen { 
    fn math_symbol() -> String { 
        String::from("Kh")
    }
}

impl Display for KhGen {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}){}", self.label, self.state.iter().map(|i| subscript(i as u8)).join("") )
    }
}

impl Gen for KhGen {}