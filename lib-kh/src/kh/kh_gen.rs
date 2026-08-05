use std::fmt::Display;
use itertools::Itertools;
use yui_core::util::format::subscript;
use yui_core::abst::{MathType, Ring, RingOps};
use yui_core::lc::{LcKey, Lc};
use yui_link::State;

use super::{KhAlgGen, KhTensor};

#[derive(Clone, Copy, Default, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct KhGen {
    state: State,
    tensor: KhTensor,
}

impl KhGen {
    pub fn new(state: State, tensor: KhTensor) -> KhGen {
        KhGen { state, tensor }
    }

    pub fn init() -> Self {
        KhGen::new(State::empty(), KhTensor::empty())
    }

    pub fn state(&self) -> &State {
        &self.state
    }

    pub fn tensor(&self) -> &KhTensor {
        &self.tensor
    }

    pub fn rel_h_deg(&self) -> isize {
        self.state.weight() as isize
    }

    pub fn rel_q_deg(&self) -> isize {
        let d = self.tensor.iter().map(|x| x.deg()).sum::<isize>();
        let r = self.tensor.len() as isize;
        let s = self.state.weight() as isize;
        d + r + s
    }

    pub fn apply_at<F, R>(&self, i: usize, f: F) -> Lc<KhGen, R>
    where F: Fn(&KhAlgGen) -> Lc<KhAlgGen, R>, R: Ring, for<'x> &'x R: RingOps<R> {
        self.tensor.apply_at(i, f).map_keys(|t| {
            Self::new(self.state, t)
        })
    }

    pub fn apply_each<F, R>(&self, f: F) -> Lc<KhGen, R>
    where F: Fn(&KhAlgGen) -> Lc<KhAlgGen, R>, R: Ring, for<'x> &'x R: RingOps<R> {
        self.tensor.apply_each(f).map_keys(|t| {
            Self::new(self.state, t)
        })
    }
}

impl Display for KhGen {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}){}", self.tensor, self.state.iter().map(|i| subscript(i as u8)).join("") )
    }
}

impl MathType for KhGen {
    fn math_symbol() -> String {
        String::from("Kh")
    }
}

impl LcKey for KhGen {}
