//! Greedy block-pivot finder for the selective build — collects an independent
//! (diagonal) block of invertible edges directly from candidate keys, with no
//! SpMat / PivotFinder (the selective deloop already minted the pivots).

use ahash::AHashSet;
use itertools::Itertools;
use yui_core::{Ring, RingOps};

use crate::tng::{LcCobTrait, TngComplex, TngComplexKey};

/// A diagonal (independent) block of invertible edges among `candidates`, as
/// `(source, target)` pairs whose eliminations are mutually independent.
pub(crate) fn find_block<R>(complex: &TngComplex<R>, candidates: &[TngComplexKey]) -> Vec<(TngComplexKey, TngComplexKey)>
where R: Ring, for<'x> &'x R: RingOps<R> {
    // cap the block so eliminate_block's contribution buffer can't blow up memory.
    const MAX_PIVOT: usize = 4096;

    let ordered = candidates.iter()
        .filter(|k| complex.contains_key(k)) // skip keys removed by a prior block this round.
        .sorted_by_key(|k| (complex.vertex(k).c_weight(), **k));

    let mut used = AHashSet::new();
    let mut block = vec![];

    for &k in ordered {
        if block.len() >= MAX_PIVOT { break }
        if used.contains(&k) { continue }
        let Some((s, t)) = pick_pivot(complex, &k, &used) else { continue };

        // keep the block diagonal: reject a pivot crossing any member.
        let crosses = block.iter().any(|(s2, t2)|
            complex.has_edge(&s, t2) || complex.has_edge(s2, &t)
        );
        if crosses { continue }

        used.insert(s);
        used.insert(t);
        block.push((s, t));
    }
    block
}

// `k`'s lightest invertible incident edge with a free other endpoint, oriented
// `(source, target)`: an out-edge `(k, l)` or an in-edge `(j, k)`.
fn pick_pivot<R>(complex: &TngComplex<R>, k: &TngComplexKey, used: &AHashSet<TngComplexKey>) -> Option<(TngComplexKey, TngComplexKey)>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let v = complex.vertex(k);
    let out = v.out_edges()
        .filter(|&l| !used.contains(l) && complex.edge(k, l).is_invertible())
        .map(|l| (*k, *l));
    let inc = v.in_edges()
        .filter(|&j| !used.contains(j) && complex.edge(j, k).is_invertible())
        .map(|j| (*j, *k));
    out.chain(inc).min_by_key(|(s, t)| (edge_weight(complex, s, t), *s, *t))
}

// fill estimate of eliminating `(s, t)`: (col nnz − 1)·(row nnz − 1).
fn edge_weight<R>(complex: &TngComplex<R>, s: &TngComplexKey, t: &TngComplexKey) -> usize
where R: Ring, for<'x> &'x R: RingOps<R> {
    let ns = complex.vertex(s).out_edges().count();
    let nt = complex.vertex(t).in_edges().count();
    (ns - 1) * (nt - 1)
}

#[cfg(test)]
mod tests {
    use yui_link::Link;
    use crate::tng::builder::{TngComplexBuilder, BuildConfig, DeloopMode, ElimMode};
    use crate::tng::LcCob;
    use super::*;

    // Deloop the trefoil but keep all invertible edges (no elim) — the fixture.
    fn delooped_trefoil() -> TngComplex<i64> {
        let l = Link::test_data("3_1");
        let config = BuildConfig { deloop_mode: DeloopMode::Greedy, elim_mode: ElimMode::None, h_range: None };
        TngComplexBuilder::from_link(&l, &0, &0, false)
            .with_config(config)
            .run()
            .into_tng_complex()
    }

    // The block is a diagonal of invertible edges: members share no endpoint and
    // have no cross edge.
    #[test]
    fn find_block_diagonal_independent() {
        let c = delooped_trefoil();
        let candidates = c.h_range().flat_map(|i| c.keys_of_deg(i).copied().collect_vec()).collect_vec();
        let block = find_block(&c, &candidates);

        assert!(!block.is_empty(), "expected a non-empty block");
        for (s, t) in &block {
            assert!(c.edge(s, t).is_invertible());
        }
        for (p, (sp, tp)) in block.iter().enumerate() {
            for (q, (sq, tq)) in block.iter().enumerate() {
                if p == q { continue }
                assert!(sp != sq && tp != tq && sp != tq && tp != sq, "shared endpoint");
                assert!(!c.has_edge(sp, tq), "cross edge ({sp},{tq}) breaks diagonality");
            }
        }
    }

    // Surviving keys with their out-edges, canonicalized for equality comparison.
    fn snapshot(c: &TngComplex<i64>) -> Vec<(TngComplexKey, Vec<(TngComplexKey, LcCob<i64>)>)> {
        let mut keys = c.h_range().flat_map(|i| c.keys_of_deg(i).copied().collect_vec()).collect_vec();
        keys.sort();
        keys.into_iter().map(|k| {
            let mut es = c.vertex(&k).out_edges().map(|l| (*l, c.edge(&k, l).clone())).collect_vec();
            es.sort_by_key(|(l, _)| *l);
            (k, es)
        }).collect()
    }

    // Eliminating an independent block in parallel == eliminating its pivots serially.
    #[test]
    fn eliminate_block_matches_serial() {
        let c = delooped_trefoil();
        let candidates = c.h_range().flat_map(|i| c.keys_of_deg(i).copied().collect_vec()).collect_vec();
        let block = find_block(&c, &candidates);
        assert!(!block.is_empty());

        let mut par = delooped_trefoil();
        par.eliminate_block(&block);

        let mut ser = delooped_trefoil();
        for (k0, k1) in &block {
            ser.eliminate(k0, k1);
        }

        assert_eq!(snapshot(&par), snapshot(&ser));
    }
}
