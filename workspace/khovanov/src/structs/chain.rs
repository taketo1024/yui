use std::fmt::Display;
use std::ops::{Mul, MulAssign};
use itertools::join;
use auto_impl_ops::auto_ops;
use yui_core::Elem;
use yui_lin_comb::{FreeGen, LinComb};
use yui_link::{State, Resolution};
use yui_utils::subscript;

use crate::KhAlgLabel;

#[derive(Debug, Clone, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct KhGen { 
    state: State,
    label: Vec<KhAlgLabel>
}

impl KhGen {
    pub fn new(state: State, label: Vec<KhAlgLabel>) -> KhGen { 
        KhGen { state, label }
    }

    pub fn init() -> Self { 
        Self::new(State::empty(), vec![])
    }

    pub fn state(&self) -> &State {
        &self.state
    }

    pub fn label(&self) -> &Vec<KhAlgLabel> {
        &self.label
    }

    pub fn label_at(&self, i: usize) -> KhAlgLabel {
        self.label[i]
    }

    pub fn q_deg(&self) -> isize { 
        let q = self.label.iter().map(|x| x.q_deg()).sum::<isize>();
        let r = self.label.len() as isize;
        let s = self.state.weight() as isize;
        q + r + s
    }

    pub fn append_state(&mut self, r: Resolution) { 
        self.state.push(r)
    }
    
    pub fn append_label(&mut self, x: KhAlgLabel) { 
        self.label.push(x)
    }
    
    pub fn append(&mut self, other: KhGen) { 
        let KhGen { state, mut label } = other;
        self.state.append(state);
        self.label.append(&mut label);
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
        let label = join(self.label.iter(), "");
        let state = join(self.state.iter().map(|i| subscript(i.as_u8() as isize)), "");
        write!(f, "{}{}", label, state)
    }
}

impl FreeGen for KhGen {}

pub type KhChain<R> = LinComb<KhGen, R>;