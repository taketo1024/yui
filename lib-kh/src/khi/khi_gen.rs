use std::fmt::Display;

use yui_core::lc::LcKey;
use yui_core::MathType;
use crate::kh::KhGen;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum KhIGen { 
    B(KhGen), Q(KhGen)
}

impl KhIGen {
    pub fn rel_h_deg(&self) -> isize {
        match self {
            KhIGen::B(x) => x.rel_h_deg(),
            KhIGen::Q(x) => x.rel_h_deg() + 1,
        }
    }

    pub fn rel_q_deg(&self) -> isize {
        match self {
            KhIGen::B(x) |
            KhIGen::Q(x) => x.rel_q_deg()
        }
    }
}

impl Display for KhIGen {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self { 
            KhIGen::B(x) => x.fmt(f),
            KhIGen::Q(x) => write!(f, "Q{}", x)
        }
    }
}

impl Default for KhIGen {
    fn default() -> Self {
        KhIGen::B(KhGen::default())
    }
}

impl MathType for KhIGen {
    fn math_symbol() -> String {
        String::new()
    }
}

impl LcKey for KhIGen {}

