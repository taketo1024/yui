use std::collections::HashMap;
use std::fmt::Display;

use itertools::Itertools;
use num_traits::Zero;
use yui_utils::map;

use crate::{KhEnhState, KhAlgGen};

use super::cob::{Cob, Bottom, Dot};
use super::mor::{Mor, MorTrait};
use super::tng::TngComp;

// element in C as a cobordism ∅ → C.
pub struct TngElem {
    mors: HashMap<KhEnhState, Mor>
}

impl TngElem { 
    pub fn new(mors: HashMap<KhEnhState, Mor>) -> Self { 
        Self {mors}
    }

    pub fn init(t: KhEnhState, c: Cob) -> Self { 
        let f = Mor::from(c);
        Self::new(map! { t => f })
    }

    pub fn deloop(&mut self, k: &KhEnhState, c: &TngComp) {
        let targets = self.mors.keys().cloned().collect_vec();

        for t in targets { 
            if !k.is_sub(&t) { 
                continue 
            }

            let f = self.mors.remove(&t).unwrap();

            let mut t0 = t.clone();
            let mut t1 = t.clone();

            t0.label.push(KhAlgGen::X);
            t1.label.push(KhAlgGen::I);

            let f0 = f.clone().cap_off(Bottom::Tgt, c, Dot::None);
            let f1 = f.clone().cap_off(Bottom::Tgt, c, Dot::Y);

            if !f0.is_zero() {
                self.mors.insert(t0, f0);
            }
            if !f1.is_zero() { 
                self.mors.insert(t1, f1);
            }
        }
    }
}

impl Display for TngElem {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mors = self.mors.iter().map(|(t, f)| { 
            format!("{}: {}", t, f)
        }).join(", ");
        write!(f, "[{}]", mors)
    }
}