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
    in_cob: Cob,                        // initial cup-cob, precomposed at the final step.
    out_cob: HashMap<TngComplexKey, LcCob<R>>, // building cob, src must always match in_cob.
    base_pt: Option<Edge>
}

impl<R> TngComplexElem<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(state: HashMap<Node, Bit>, in_cob: Cob, base_pt: Option<Edge>) -> Self {
        let k0 = TngComplexKey::init();
        let f0 = LcCob::from(Cob::empty());
        let out_cob = hashmap! { k0 => f0 };
        let elem = Self { state, in_cob, out_cob, base_pt };
        elem.verify();
        elem
    }

    pub fn state(&self) -> &HashMap<Node, Bit> {
        &self.state
    }

    // h-degree relative to `deg_shift`: the weight (number of 1-resolutions) of the state.
    pub fn rel_h_deg(&self) -> isize {
        self.state.values().filter(|b| b.is_one()).count() as isize
    }

    pub fn base_pt(&self) -> Option<Edge> {
        self.base_pt
    }

    pub fn in_cob(&self) -> &Cob {
        &self.in_cob
    }

    pub fn out_cob(&self) -> &HashMap<TngComplexKey, LcCob<R>> {
        &self.out_cob
    }

    pub(crate) fn set_out_cob<I>(&mut self, out_cob: I)
    where I: IntoIterator<Item = (TngComplexKey, LcCob<R>)> {
        self.out_cob = out_cob.into_iter().collect();
        self.verify();
    }

    // Take the `out_cob`, let `f` rebuild or edit it, store the result and re-check the invariant.
    pub(crate) fn modify_out_cob<F>(&mut self, f: F)
    where F: FnOnce(HashMap<TngComplexKey, LcCob<R>>) -> HashMap<TngComplexKey, LcCob<R>> {
        let out = std::mem::take(&mut self.out_cob);
        self.out_cob = f(out);
        self.verify();
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

    // Invariant (debug-only): no zero `out_cob` value; all `out_cob` terms share one source; every
    // closed circle of that source is one of in_cob's target circles.
    fn verify(&self) {
        if !cfg!(debug_assertions) {
            return;
        }

        // 0. no out_cob value is zero (a zero entry has no source — invisible to (1)).
        debug_assert!(self.out_cob.values().all(|f| f.any_term().is_some()), "verify: zero out_cob entry");

        // 1. all out_cob terms reconstruct the same source.
        let mut out_srcs = self.out_cob.values().flat_map(|f| f.iter().map(|(c, _)| c.reconst_src()));
        let Some(out_src) = out_srcs.next() else { return };
        debug_assert!(out_srcs.all(|s| s == out_src), "verify: out_cob terms have differing sources");

        // 2. every closed circle in out_cob's source is one of in_cob's target circles (open arcs,
        // which may later close up, are not checked).
        let in_tgt = self.in_cob.reconst_tgt();
        debug_assert!(out_src.comps().filter(|c| c.is_circle()).all(|c| in_tgt.contains(c)), "verify: out_cob source circle not in in_cob target");
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
                    let dot = if col.is_a() == o { Dot::X } else { Dot::Y };
                    CobComp::cup(t).add_dot(dot)
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