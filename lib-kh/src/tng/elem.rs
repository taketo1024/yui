use std::collections::HashMap;
use std::fmt::Display;

use itertools::Itertools;
use maplit::hashmap;
use yui_core::bitseq::Bit;
use yui_core::{Ring, RingOps};
use yui_link::{Edge, Link, Node};

use super::tng::TngComp;
use super::cob::{Dot, Cob, CobComp, LcCobTrait, LcCob};
use super::complex::TngComplexKey;
use crate::kh::KhChain;
use crate::ext::LinkExt;

#[derive(Clone)]
pub struct TngComplexElem<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    state: HashMap<Node, Bit>,
    in_cob: Cob,                        // initial cob, precomposed at the final step.
    out_cob: HashMap<TngComplexKey, LcCob<R>>, // building cob, src must always match init_cob. 
    base_pt: Option<Edge>
}

impl<R> TngComplexElem<R> 
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new(state: HashMap<Node, Bit>, in_cob: Cob, base_pt: Option<Edge>) -> Self { 
        let k0 = TngComplexKey::init();
        let f0 = LcCob::from(Cob::empty());
        let out_cob = hashmap! { k0 => f0 };
        Self{ state, in_cob, out_cob, base_pt }
    }

    pub fn state(&self) -> &HashMap<Node, Bit> {
        &self.state
    }

    pub fn base_pt(&self) -> Option<Edge> {
        self.base_pt
    }

    pub fn out_cob(&self) -> &HashMap<TngComplexKey, LcCob<R>> {
        &self.out_cob
    }

    pub fn out_cob_mut(&mut self) -> &mut HashMap<TngComplexKey, LcCob<R>> {
        &mut self.out_cob
    }

    pub fn is_evalable(&self) -> bool {
        let init = LcCob::from(self.in_cob.clone());
        self.out_cob.values().all(|c| init.is_stackable(c)) && 
        self.out_cob.values().all(|c| c.iter().all(|(c, _)| c.comps().all(|c| c.tgt().is_empty())))
    }

    pub fn eval(&self, h: &R, t: &R) -> KhChain<R> {
        assert!(self.is_evalable());

        let init = LcCob::from(self.in_cob.clone());
        let eval = self.out_cob.iter().map(|(k, retr)| {
            let x = k.as_gen();
            let f = retr * &init;
            let r = f.eval(h, t);
            (x, r)
        }).collect::<KhChain<R>>();

        eval
    }

    pub fn canon_cycles(l: &Link, base_pt: Option<Edge>) -> Vec<Self> { 
        assert!(l.is_knot());
        assert!(l.base_pt().is_some());

        let reduced = base_pt.is_some();
        let circles = l.colored_seifert_circles();

        let crossings = l.nodes().filter(|x| x.is_crossing()).cloned();
        let state = l.seifert_state();
        let state_map = Iterator::zip(crossings.into_iter(), state.iter()).collect::<HashMap<_, _>>();

        let ori = if reduced { 
            vec![true]
        } else { 
            vec![true, false]
        };

        let cycles = ori.into_iter().map(|o| {
            let cob = Cob::new(
                circles.iter().map(|(circ, col)| {
                    let marked = base_pt.map(|b| circ.contains(b)).unwrap_or(false);
                    let t = TngComp::from_path(circ.clone(), marked);
                    let mut cup = CobComp::cup(t);
                    let dot = if col.is_a() == o { 
                        Dot::X 
                    } else { 
                        Dot::Y 
                    };
                    cup.add_dot(dot);
                    cup
                })
            );
            TngComplexElem::new(state_map.clone(), cob, base_pt)
        }).collect();

        cycles
    }
}

impl<R> Display for TngComplexElem<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mors = self.out_cob.iter().sorted_by_key(|&(&k, _)| k).map(|(k, f)| { 
            format!("{}: {}", k, f)
        }).join(", ");
        write!(f, "[{}]", mors)
    }
}