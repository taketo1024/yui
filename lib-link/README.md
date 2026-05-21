## yui-link

[![crates.io](https://img.shields.io/crates/v/yui-link.svg)](https://crates.io/crates/yui-link)
[![docs.rs](https://docs.rs/yui-link/badge.svg)](https://docs.rs/yui-link)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Knots and links for the [`yui`](https://github.com/taketo1024/yui) workspace: planar diagrams, braid words, and a couple of derived invariants. Used by [`yui-kh`](https://crates.io/crates/yui-kh) (Khovanov homology) and other knot-homology layers.

## Layout

```text
src/
├── link/
│   ├── link.rs       — `Link`, `Edge`, `State`, `XCode`
│   ├── node.rs       — `Node`, `NodeType` {XL, XR, V, H}, `NodeOri`
│   ├── path.rs       — `Path`: an arc or a circle of edges
│   ├── graph.rs      — `seifert_graph`
│   └── inv_link.rs   — `InvLink`: involutive link
├── braid.rs          — `Braid`, `Generator`
├── misc/
│   └── jones.rs      — `jones_polynomial`
└── test_data.rs      — hard-coded PD codes / braid words (test-only)
```

## Key types

### `Link`

A knot or link given by a planar diagram. Internally a `Vec<Node>` plus a `HashSet<Edge>`, where `Edge = usize` is an opaque edge id. Construct from:

- a PD code: `Link::from_pd_code([[1,4,2,5], ...])`,
- explicit nodes: `Link::from_nodes([...])`,
- a JSON file under `~/.yui/data/links/`: `Link::load("3_1")`.

Provides `n_crossings`, `n_comps`, `comps`, `writhe`, `mirror`, `resolve_at` / `resolve_by` (Khovanov-style 0/1-smoothings), `seifert_state`, `seifert_circles`, and traversal helpers.

### `Node`

A vertex in a planar diagram — either a crossing or a smoothing — together with an orientation:

```text
   3   2         3   2         3   2         3   2
    \ /           \ /           \ /           \_/
     \    = XL,    /    = XR,   | |   = V,     _    = H,
    / \           / \           / \           / \
   0   1         0   1         0   1         0   1
```

Carries a `NodeType` (XL / XR / V / H), a `NodeOri` (`↑ ↓ ← →` or none), and `[Edge; 4]` of incident edges in the order shown. `resolve(bit)` turns a crossing (XL/XR) into a smoothing (V/H); `mirror()` swaps XL ↔ XR; `sign()` returns the crossing sign as `Option<Sign>` (from (type, orientation)).

### `Path`

An ordered list of `Edge`s with a `closed: bool` flag — either an *arc* (open path) or a *circle* (closed loop). Constructed via `Path::arc(edges)` / `Path::circ(edges)` / `Path::new(edges, closed)`. Returned as the component type by `Link::comps()` and `Link::seifert_circles()`. Methods include `is_arc`, `is_circle`, `ends`, `contains`, `connect` (concatenate two arcs sharing an endpoint), `reduce` (canonicalize the edge sequence), and `unori_eq` (equality up to orientation).

### `Braid`

A word in the Artin generators of the braid group: a strand count plus a `Vec<Generator>`, where `Generator(i32)` encodes the signed half-twist index (positive = right-handed). `closure()` produces the corresponding `Link`. Other operations: `inv`, `Mul` for concatenation, `display()` for ASCII rendering, `Braid::load(name)` for JSON-backed braid words.

### `InvLink`

An *involutive link*: a `Link` together with an involution on it — an edge bijection `e ↦ e'`, the induced node bijection, and an optional axis base point. Build directly via `InvLink::new(link, edge_map_fn, base_pt)`. `InvLink::from_symmetric_pd_code` is a convenience constructor for *strongly invertible* knots, where the involution acts on edges as `e ↦ (n+1-e) mod n + 1`; `InvLink::load(name)` returns hard-coded entries for a small table of such knots (`3_1`, `4_1`, `5_1`, `5_2a/b`, …, `7_7a/b`). Most read-only `Link` methods are delegated.

### Derived invariants

- `seifert_graph(&Link) -> petgraph::Graph<Path, usize>` — vertices are Seifert circles, edges are crossings.
- `jones_polynomial(&Link) -> LPoly<'q', i32>` — Kauffman-bracket computation summing over all `2^n` resolutions.

## Conventions

- **PD code.** Each crossing is `[a, b, c, d]` ordered counter-clockwise from the lower-left, with `a → c` the incoming under-strand. See:
  - KnotAtlas — [katlas.org/wiki/Planar_Diagrams](https://katlas.org/wiki/Planar_Diagrams)
  - KnotInfo — [knotinfo.org/descriptions/pd_notation.html](https://knotinfo.org/descriptions/pd_notation.html)
- **Edge ids.** `Edge = usize`. Only used to identify endpoint coincidence; ids don't have to be contiguous or start at 0, but `InvLink::from_symmetric_pd_code` expects `1..=n_edges` so the involution formula is well-defined.

## Quick example

```rust
use yui_link::Link;

// Trefoil knot from a PD code.
let l = Link::from_pd_code([[1,4,2,5], [3,6,4,1], [5,2,6,3]]);
assert_eq!(l.n_crossings(), 3);
assert_eq!(l.writhe(),     -3);
assert_eq!(l.n_comps(),     1);
```

## Feature flags

- `test-utils` — exposes the hard-coded `Link::test_data(name)` / `Braid::test_data(name)` constructors to downstream crates (behind `#[cfg(test)]` by default).

## License

This library is licensed under the [MIT License](https://opensource.org/licenses/MIT).

---

*This README was generated by [Claude](https://www.anthropic.com/claude).*
