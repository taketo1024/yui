use std::fmt::Display;
use itertools::Itertools;
use yui_core::util::format::subscript;
use yui_core::{Elem, Ring, RingOps};
use yui_core::lc::{LcKey, Lc};
use yui_link::State;

use super::{KhAlgGen, KhTensor};

#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct KhState {
    pub state: State,
    pub tensor: KhTensor,
    pub deg_shift: (isize, isize)
}

impl KhState {
    pub fn new(state: State, tensor: KhTensor, deg_shift: (isize, isize)) -> KhState {
        KhState { state, tensor, deg_shift }
    }

    pub fn init() -> Self {
        KhState::new(State::empty(), KhTensor::empty(), (0, 0))
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

impl Elem for KhState {
    fn math_symbol() -> String {
        String::from("Kh")
    }
}

impl Display for KhState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}){}", self.tensor, self.state.iter().map(|i| subscript(i as u8)).join("") )
    }
}

impl LcKey for KhState {}

pub type KhChain<R> = Lc<KhState, R>;

pub trait KhChainExt {
    fn h_deg(&self) -> isize;
    fn q_deg(&self) -> isize;
}

impl<R> KhChainExt for KhChain<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn h_deg(&self) -> isize {
        self.keys().map(|x| x.h_deg()).min().unwrap_or(0)
    }

    fn q_deg(&self) -> isize {
        self.keys().map(|x| x.q_deg()).min().unwrap_or(0)
    }
}
