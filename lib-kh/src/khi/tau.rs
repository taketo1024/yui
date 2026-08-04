//! The chain-level involution `τ : CKh(D) → CKh(D)` induced by the symmetry
//! of a strongly invertible link diagram (Section 2.1 of the reference).
//! Used to build the involutive Khovanov complex `CKhI = Cone(1 + τ)`.
//!
//! Reference:
//! - T. Sano, "Involutive Khovanov homology and equivariant knots",
//!   Algebr. Geom. Topol. 25 (2025), 5059–5111.
//!   <https://doi.org/10.2140/agt.2025.25.5059>, <https://arxiv.org/abs/2404.08568>

use std::collections::HashMap;

use yui_link::{InvLink, State};

use crate::kh::{KhAlgGen, KhGen, KhTensor};

/// Build the chain-level involution `τ : KhGen → KhGen` induced by the link
/// involution of `l`.
///
/// The returned closure precomputes per-state circle permutations from the
/// link's edge involution, then applies them to each generator.
pub(crate) fn tau_map(l: &InvLink) -> impl Fn(&KhGen) -> KhGen + Send + Sync + 'static {
    assert!(l.nodes().all(|x| x.is_crossing()));

    let n = l.n_crossings();
    let inner = l.inner();

    let x_index = (0..n).map(|i| (l.node(i), i)).collect::<HashMap<_, _>>();

    let state_map = State::generate(n).map(|s| {
        let t = State::from_iter((0..n).map(|i| {
            let x = l.node(i);
            let tx = l.inv_node(x);
            let ti = x_index[tx];
            s[ti]
        }));
        (s, t)
    }).collect::<HashMap<_, _>>();

    // For each state, the circles sorted by `min_edge` (matches `KhCubeVertex::new`).
    let circles_at = |s: State| {
        let mut cs = inner.resolve_by(&s).comps();
        cs.sort_by_key(|c| c.min_edge());
        cs
    };

    let label_map = State::generate(n).map(|s| {
        let t = state_map[&s];
        let v_circles = circles_at(s);
        let w_circles = circles_at(t);

        debug_assert_eq!(v_circles.len(), w_circles.len());

        let map = (0..v_circles.len()).map(|i| {
            let c0 = &v_circles[i];
            let e = l.inv_edge(c0.min_edge());
            let j = w_circles.iter().position(|c| c.contains(e)).unwrap();
            (i, j)
        }).collect::<HashMap<_, _>>();

        (s, map)
    }).collect::<HashMap<_, _>>();

    move |x: &KhGen| {
        let s = *x.state();
        let t = state_map[&s];
        let map = &label_map[&s];
        let label = x.tensor();

        let mut seq = vec![KhAlgGen::I; label.len()];
        for (i, e) in label.iter().enumerate() {
            if e.is_X() {
                let j = map[&i];
                seq[j] = KhAlgGen::X;
            }
        }

        KhGen::new(t, KhTensor::from_iter(seq))
    }
}

#[cfg(test)]
mod tests {
    use yui_link::InvLink;
    use crate::kh::{KhGen, KhTensor, KhAlgGen::{I, X}};
    use yui_link::State;
    use super::*;

    #[test]
    fn tau_trefoil_states() {
        let l = InvLink::from_symmetric_pd_code([[5,3,6,2],[1,5,2,4],[3,1,4,6]]);
        let tau = tau_map(&l);

        // empty tensor decoration to focus on state.
        let empty = KhTensor::empty();
        let probe = |s: State| {
            let r = l.inner().resolve_by(&s).comps().len();
            let mut tensor = KhTensor::empty();
            for _ in 0..r { tensor.push(I); }
            tau(&KhGen::new(s, tensor)).state().clone()
        };
        let _ = empty; // silence unused

        assert_eq!(probe(State::from([0,0,0])), State::from([0,0,0]));
        assert_eq!(probe(State::from([1,0,0])), State::from([1,0,0]));
        assert_eq!(probe(State::from([0,1,0])), State::from([0,0,1]));
        assert_eq!(probe(State::from([0,0,1])), State::from([0,1,0]));
        assert_eq!(probe(State::from([1,1,0])), State::from([1,0,1]));
        assert_eq!(probe(State::from([1,0,1])), State::from([1,1,0]));
        assert_eq!(probe(State::from([0,1,1])), State::from([0,1,1]));
        assert_eq!(probe(State::from([1,1,1])), State::from([1,1,1]));
    }

    #[test]
    fn tau_trefoil_labels() {
        let l = InvLink::from_symmetric_pd_code([[5,3,6,2],[1,5,2,4],[3,1,4,6]]);
        let tau = tau_map(&l);

        let apply = |s: State, t: KhTensor| {
            let g = tau(&KhGen::new(s, t));
            *g.tensor()
        };

        assert_eq!(apply(State::from([0,0,0]), KhTensor::from([I, X])), KhTensor::from([I, X]));
        assert_eq!(apply(State::from([0,1,0]), KhTensor::from([X])),    KhTensor::from([X]));
        assert_eq!(apply(State::from([1,1,0]), KhTensor::from([I, X])), KhTensor::from([I, X]));
        assert_eq!(apply(State::from([1,1,1]), KhTensor::from([I, I, X])), KhTensor::from([I, X, I]));
    }

    #[test]
    fn tau_is_involution() {
        let l = InvLink::from_symmetric_pd_code([[5,3,6,2],[1,5,2,4],[3,1,4,6]]);
        let tau = tau_map(&l);

        for s in State::generate(3) {
            let r = l.inner().resolve_by(&s).comps().len();
            for label in KhTensor::generate(r) {
                let x = KhGen::new(s, label);
                let txx = tau(&tau(&x));
                assert_eq!(txx, x, "τ² should be identity, but τ²({x:?}) = {txx:?}");
            }
        }
    }
}
