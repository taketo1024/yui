use std::fmt::Display;

use yui_core::lc::LcKey;
use yui_core::MathType;
use crate::kh::KhGen;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum KhIState { 
    B(KhGen), Q(KhGen)
}

impl KhIState { 
    pub fn h_deg(&self) -> isize { 
        match self {
            KhIState::B(x) => x.h_deg(),
            KhIState::Q(x) => x.h_deg() + 1,
        }
    }

    pub fn q_deg(&self) -> isize { 
        match self {
            KhIState::B(x) | 
            KhIState::Q(x) => x.q_deg()
        }
    }
}

impl Display for KhIState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self { 
            KhIState::B(x) => x.fmt(f),
            KhIState::Q(x) => write!(f, "Q{}", x)
        }
    }
}

impl Default for KhIState {
    fn default() -> Self {
        KhIState::B(KhGen::default())
    }
}

impl MathType for KhIState {
    fn math_symbol() -> String {
        String::new()
    }
}

impl LcKey for KhIState {}

