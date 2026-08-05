use std::ops::RangeInclusive;

use num_traits::Zero;
use yui_core::bitseq::Bit;
use yui_core::{Ring, RingOps};
use yui_link::{Edge, Node, Path};

use itertools::Itertools;

use crate::kh::{KhAlgGen, KhChain, KhGen};
use crate::tng::{Tng, TngComp, Cob, End, Dot, LcCob, LcCobTrait, TngComplex, TngComplexElem, TngComplexKey, circles_of, label_assignments, cap_circles, expanded_key};

// Owns the canonical-cycle elements and transforms them in lockstep with the `TngComplex`: each
// element's `out_cob` (cap cobordism) tracks the complex via the same append/deloop/eliminate.
pub(crate) struct TngElemBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    elements: Vec<TngComplexElem<R>>,
}

impl<R> TngElemBuilder<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub(crate) fn new() -> Self {
        Self { elements: vec![] }
    }

    pub fn set<I>(&mut self, elements: I)
    where I: IntoIterator<Item = TngComplexElem<R>> {
        self.elements = elements.into_iter().collect();
    }

    pub fn take(&mut self) -> Vec<TngComplexElem<R>> {
        std::mem::take(&mut self.elements)
    }

    pub fn content(&self) -> &[TngComplexElem<R>] {
        &self.elements
    }

    pub(crate) fn retain<F>(&mut self, f: F)
    where F: FnMut(&TngComplexElem<R>) -> bool {
        self.elements.retain(f);
    }

    // Bilinearly merge `other` (a parallel chunk's elements) into self: each pair shares `state`
    // and `in_cob`, and their `out_cob`s combine over the product of keys via horizontal composition.
    pub(crate) fn merge(&mut self, other: Vec<TngComplexElem<R>>) {
        if other.is_empty() {
            return; // nothing to merge (e.g. single-crossing append, or no canon cycles)
        }
        assert_eq!(self.elements.len(), other.len());
        for (e, o) in self.elements.iter_mut().zip(other) {
            Self::merge_into(e, &o);
        }
    }

    fn merge_into(e: &mut TngComplexElem<R>, o: &TngComplexElem<R>) {
        assert_eq!(e.state(), o.state());
        assert_eq!(e.in_cob(), o.in_cob());

        e.modify_out_cob(|lhs| lhs.iter().flat_map(|(k, f)| {
            o.out_cob().iter().map(move |(l, g)| {
                (k + l, f.apply_bilin(g, |c1, c2| c1.connect(c2)))
            })
        }).collect());
    }

    pub(crate) fn append_node(&mut self, x: &Node) {
        for e in self.elements.iter_mut() {
            Self::append_node_to(e, x);
        }
    }

    fn append_node_to(e: &mut TngComplexElem<R>, x: &Node) {
        let (t, r) = if x.is_crossing() {
            let r = e.state()[x];
            (Tng::from_resolved(&x.resolve(r), e.base_pt()), Some(r))
        } else {
            assert!(x.is_resolved());
            (Tng::from_resolved(x, e.base_pt()), None)
        };
        Self::connect_id(e, &t, r);
    }

    pub(crate) fn insert_loop(&mut self, c: Edge) {
        for e in self.elements.iter_mut() {
            Self::insert_loop_to(e, c);
        }
    }

    fn insert_loop_to(e: &mut TngComplexElem<R>, c: Edge) {
        let marked = e.base_pt() == Some(c);
        let t = Tng::from(TngComp::from_path(Path::circ([c]), marked));
        Self::connect_id(e, &t, None);
    }

    // Append the identity cobordism on `t` to the out_cob, pushing the resolution bit `r` onto each key.
    fn connect_id(e: &mut TngComplexElem<R>, t: &Tng, r: Option<Bit>) {
        let id = Cob::id(t);
        e.modify_out_cob(|out| out.into_iter().map(|(mut k, f)| {
            if let Some(r) = r {
                k.state.push(r);
            }
            (k, f.connect(&id))
        }).collect());
    }

    pub(crate) fn deloop(&mut self, k: &TngComplexKey, c: &TngComp) {
        for e in self.elements.iter_mut() {
            Self::deloop_from(e, k, c);
        }
    }

    // note: capping can annihilate the cobordism (e.g. a Y-dotted cap over 𝔽₂); a zero entry must
    // not survive in `out_cob` — its phantom key made `d_sym` non-deterministic.
    fn deloop_from(e: &mut TngComplexElem<R>, k: &TngComplexKey, c: &TngComp) {
        let marked = e.base_pt().map(|b| c.contains(b)).unwrap_or(false);
        e.modify_out_cob(|mut cob| {
            let Some(f) = cob.remove(k) else { return cob };

            let k0 = k + KhAlgGen::X;
            let f0 = f.clone().cap_off(End::Tgt, c, None);
            if !f0.is_zero() {
                cob.insert(k0, f0);
            }

            if !marked {
                let k1 = k + KhAlgGen::I;
                let f1 = f.cap_off(End::Tgt, c, Some(Dot::Y));
                if !f1.is_zero() {
                    cob.insert(k1, f1);
                }
            }
            cob
        });
    }

    pub(crate) fn eliminate(&mut self, complex: &TngComplex<R>, i: &TngComplexKey, j: &TngComplexKey) {
        for e in self.elements.iter_mut() {
            Self::eliminate_from(e, complex, i, j);
        }
    }

    //  Gaussian Elimination
    //
    //       a
    //  v0 - - -> v1         .             .
    //     \   / b
    //       /         ==>
    //     /   \ c              d - ca⁻¹b
    //  w0 -----> w1         w0 ---------> w1
    //       d
    fn eliminate_from(e: &mut TngComplexElem<R>, complex: &TngComplex<R>, i: &TngComplexKey, j: &TngComplexKey) {
        debug_assert!(complex.has_edge(i, j));

        e.modify_out_cob(|mut cob| {
            // mors into i can be simply dropped.
            cob.remove(i);

            // mors into j must be redirected by -ca⁻¹.
            let Some(b) = cob.remove(j) else { return cob };

            let (h, t) = complex.ht();
            let ainv = complex.edge(i, j).inv().unwrap();
            let ainv_b = b.stack(&ainv);

            for k in complex.vertex(i).out_edges() {
                if k == j { continue }

                let c = complex.edge(i, k);
                let c_ainv_b = ainv_b.stack(c).reduce(h, t);
                let s = if let Some(d) = cob.remove(k) {
                    d - c_ainv_b
                } else {
                    -c_ainv_b
                };

                if !s.is_zero() {
                    cob.insert(*k, s);
                }
            }
            cob
        });
    }

    pub(crate) fn eval(&self, h: &R, t: &R) -> Vec<KhChain<R>> {
        self.elements.iter().map(|e| e.eval(h, t)).collect()
    }

    // Like `eval`, but expands each vertex's remaining circles over their label assignments (the same
    // pairing as `into_raw_complex`) — so it also works on an un-delooped complex (`no_full_deloop`).
    // Keys/generators dropped by the q-filter are skipped, matching `into_raw_complex_filtered`.
    pub(crate) fn eval_with(&self, c: &TngComplex<R>, h: &R, t: &R, q_range: Option<RangeInclusive<isize>>) -> Vec<KhChain<R>> {
        let q_shift = c.deg_shift().1;
        let in_window = |g: &KhGen|
            q_range.as_ref().is_none_or(|r| r.contains(&(q_shift + g.rel_q_deg())));

        self.elements.iter().map(|e| {
            let init = LcCob::from(e.in_cob().clone());
            e.out_cob().iter().filter(|(k, _)| c.contains_key(k)).flat_map(|(k, retr)| {
                let circles = circles_of(c.vertex(k).tng());
                label_assignments(&circles).into_iter().filter_map(|b| {
                    let g = expanded_key(k, &b).as_gen();
                    if !in_window(&g) {
                        return None;
                    }
                    let gc = cap_circles(retr.clone(), End::Tgt, &circles, &b, h, t);
                    let x = init.stack(&gc).eval(h, t);
                    Some((g, x))
                }).collect_vec()
            }).collect()
        }).collect()
    }
}
