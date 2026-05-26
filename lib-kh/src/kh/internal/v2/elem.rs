use std::collections::HashMap;
use std::fmt::Display;

use itertools::Itertools;
use maplit::hashmap;
use yui_core::bitseq::Bit;
use yui_core::{Ring, RingOps};
use yui_link::{Edge, Link, Node};

use crate::kh::internal::v2::cob::CobComp;
use crate::kh::{KhAlgGen, KhChain};

use super::cob::{Bottom, Dot, Cob, LcCobTrait, LcCob};
use super::tng::{Tng, TngComp};
use crate::ext::LinkExt;
use super::complex::TngComplexKey;

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

    pub fn append_node(&mut self, x: &Node) { 
        if x.is_crossing() { 
            self.append_crossing(x)
        } else { 
            self.append_arcs(x)
        }
    }

    fn append_crossing(&mut self, x: &Node) {
        assert!(x.is_crossing());

        let r = self.state[x];
        let a = x.resolve(r);
        let t = Tng::from_resolved(&a);

        self.connect_id_cob(&t, Some(r));
    }

    fn append_arcs(&mut self, x: &Node) {
        assert!(x.is_resolved());

        let t = Tng::from_resolved(x);
        self.connect_id_cob(&t, None);
    }

    pub fn insert_loop(&mut self, c: Edge) { 
        let t = Tng::from(TngComp::circ([c]));
        self.connect_id_cob(&t, None);
    }

    fn connect_id_cob(&mut self, t: &Tng, r: Option<Bit>) { 
        let id = Cob::id(&t);

        let mors = std::mem::take(&mut self.out_cob);
        self.out_cob = mors.into_iter().map(|(mut k, f)| {
            if let Some(r) = r { 
                k.state.push(r);
            }
            let f = f.connect(&id);
            (k, f)
        }).collect();
    }

    pub fn deloop(&mut self, k: &TngComplexKey, c: &TngComp) {
        let Some(f) = self.out_cob.remove(k) else { return };
        let marked = self.base_pt.map(|e| c.contains(e)).unwrap_or(false);

        let k0 = k + KhAlgGen::X;
        let f0 = f.clone().cap_off(Bottom::Tgt, c, Dot::None);
        self.out_cob.insert(k0, f0);

        if !marked { 
            let k1 = k + KhAlgGen::I;
            let f1 = f.cap_off(Bottom::Tgt, c, Dot::Y);
            self.out_cob.insert(k1, f1);    
        }
    }

    pub fn insert_cob(&mut self, k: TngComplexKey, v: LcCob<R>) {
        self.out_cob.insert(k, v);
    }

    pub fn remove_cob(&mut self, k: &TngComplexKey) -> Option<LcCob<R>> { 
        self.out_cob.remove(k)
    }

    pub fn modify<F>(&mut self, f: F)
    where F: Fn(TngComplexKey, LcCob<R>) -> (TngComplexKey, LcCob<R>) { 
        let retr_cob = std::mem::take(&mut self.out_cob);
        self.out_cob = retr_cob.into_iter().map(|(k, cob)|
            f(k, cob)
        ).collect();
    }

    pub fn is_evalable(&self) -> bool { 
        let init = LcCob::from(self.in_cob.clone());
        self.out_cob.values().all(|c| init.is_stackable(c)) && 
        self.out_cob.values().all(|c| c.iter().all(|(c, _)| c.tgt().is_empty()))
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
                    let t = TngComp::from(circ.clone());
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