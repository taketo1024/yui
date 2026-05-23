use std::cmp::Ordering;
use std::fmt::Display;
use std::hash::{Hash, Hasher};
use itertools::Itertools;
use yui_core::util::format::subscript;
use yui_core::{MathType, Ring, RingOps};
use yui_core::lc::{LcKey, Lc};
use yui_link::State;

use super::{KhAlgGen, KhTensor};

#[derive(Clone, Copy, Default, Debug)]
pub struct KhGen {
    state: State,
    tensor: KhTensor,
    deg_shift: (isize, isize),
}

impl KhGen {
    pub fn new(state: State, tensor: KhTensor, deg_shift: (isize, isize)) -> KhGen {
        KhGen { state, tensor, deg_shift }
    }

    pub fn init() -> Self {
        KhGen::new(State::empty(), KhTensor::empty(), (0, 0))
    }

    pub fn state(&self) -> &State {
        &self.state
    }

    pub fn tensor(&self) -> &KhTensor {
        &self.tensor
    }

    pub fn deg_shift(&self) -> (isize, isize) {
        self.deg_shift
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

    pub fn apply_at<F, R>(&self, i: usize, f: F) -> KhChain<R>
    where F: Fn(&KhAlgGen) -> Lc<KhAlgGen, R>, R: Ring, for<'x> &'x R: RingOps<R> {
        self.tensor.apply_at(i, f).map_keys(|t| {
            Self::new(self.state, t, self.deg_shift)
        })
    }

    pub fn apply_each<F, R>(&self, f: F) -> KhChain<R>
    where F: Fn(&KhAlgGen) -> Lc<KhAlgGen, R>, R: Ring, for<'x> &'x R: RingOps<R> {
        self.tensor.apply_each(f).map_keys(|t| {
            Self::new(self.state, t, self.deg_shift)
        })
    }
}

impl PartialEq for KhGen {
    fn eq(&self, other: &Self) -> bool {
        self.state == other.state && self.tensor == other.tensor
    }
}

impl Eq for KhGen {}

impl Hash for KhGen {
    fn hash<H: Hasher>(&self, hasher: &mut H) {
        self.state.hash(hasher);
        self.tensor.hash(hasher);
    }
}

impl PartialOrd for KhGen {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for KhGen {
    fn cmp(&self, other: &Self) -> Ordering {
        self.state.cmp(&other.state).then_with(|| self.tensor.cmp(&other.tensor))
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

pub type KhChain<R> = Lc<KhGen, R>;
