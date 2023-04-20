use std::fmt::Display;
use std::iter::zip;

use itertools::Itertools;

use crate::KhEnhState;

use super::cob::Cob;
use super::mor::Mor;

// element in C as a cobordism ∅ → C.
pub struct TngElem {
    targets: Vec<KhEnhState>,
    mors: Vec<Mor>
}

impl TngElem { 
    pub fn new(targets: Vec<KhEnhState>, mors: Vec<Mor>) -> Self { 
        assert_eq!(targets.len(), mors.len());
        Self {targets, mors}
    }

    pub fn init(t: KhEnhState, c: Cob) -> Self { 
        let f = Mor::from(c);
        Self::new(vec![t], vec![f])
    }
}

impl Display for TngElem {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let maps = zip(self.targets.iter(), self.mors.iter()).map(|(t, f)| { 
            format!("{}: {}", t, f)
        }).join(", ");
        write!(f, "[{}]", maps)
    }
}