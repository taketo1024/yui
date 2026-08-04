## yui-link

[![crates.io](https://img.shields.io/crates/v/yui-link.svg)](https://crates.io/crates/yui-link)
[![docs.rs](https://docs.rs/yui-link/badge.svg)](https://docs.rs/yui-link)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Knots and links for the [`yui`](https://github.com/taketo1024/yui) workspace: planar diagrams, braid words, and a couple of derived invariants. Used by [`yui-kh`](https://crates.io/crates/yui-kh) (Khovanov homology) and other knot-homology layers.

## Layout

```text
src/
├── link/
│   ├── link.rs       — `Link`: the type itself — accessors, components, traversal
│   ├── link_ops.rs   — operations *on* a link: mirror, connected sum, resolutions, Seifert's algorithm
│   ├── construct.rs  — links built *from patterns*: twist knots, pretzels, cables, Whitehead doubles
│   ├── pd_code.rs    — PD-code conversion and `Link::load`
│   ├── builder.rs    — `LinkBuilder`: assemble a diagram port by port
│   ├── node.rs       — `Node`: vertex of a diagram (crossing or smoothing)
│   └── path.rs       — `Path`: an arc or a circle of edges
├── inv_link/
│   ├── inv_link.rs   — `InvLink`: involutive link, with mirror and connected sum
│   └── construct.rs  — the equivariant constructions
├── braid/
│   ├── braid.rs      — `Braid`: word in the Artin generators
│   └── braid_gen.rs  — `BraidGen`: signed Artin generator
├── misc/
│   └── jones.rs      — `jones_polynomial`, `det`
└── test_data.rs      — hard-coded PD codes / braid words (test-only)
resources/
└── inv_link/         — bundled symmetric PD codes for `InvLink::load`
    └── *.json        — 3_1, 4_1, 5_1, 5_2a/b, …, 7_7a/b
```

## Key types

### `Link`

A knot or link given by a planar diagram. Internally a `Vec<Node>` plus a `Vec<Edge>` of free loops and an optional base point. `Edge = usize` is an opaque edge id. Construct from:

- a PD code: `Link::from_pd_code([[1,4,2,5], ...])`,
- explicit nodes and loops: `Link::new(nodes, loops)`, or the wrapper `Link::from_nodes(nodes)` for no loops,
- a JSON file under `<DATA_DIR>/links/`: `Link::load("3_1")`,
- predefined shapes: `Link::empty()`, `Link::unknot()`, `Link::unlink(n)`.

Free loops are closed components without crossings. They participate in `comps()`, `n_comps()`, `n_edges()`, and `edges()` as expected — each loop shows up as a one-element closed `Path` in `comps()`.

Each link carries a `base_pt: Option<Edge>`, defaulting to the minimum edge of the diagram (or `None` for the empty link). Set explicitly via `with_base_pt(e)` (consuming builder; asserts `e` is a real edge of the diagram).

Accessors: `n_crossings`, `n_comps`, `comps`, `writhe`, `loops`, `n_loops`, `base_pt`, `with_base_pt`, `reindexed`, traversal helpers.

Operations on a diagram: `mirror`, `conn_sum` / `conn_sum_at`, `cc_at` (crossing change), `resolve_at` / `resolve_by` (Khovanov-style 0/1-smoothings), `seifert_state`, `seifert_circles`, `seifert_graph`.

Constructions are associated functions, taking the companion link (if any) as an argument: `Link::twist_knot(n)`, `Link::pretzel(a, b, c)`, `Link::cable2(&l)` (blackboard-framed 2-cable), `Link::whitehead_double(&l, positive, tw)` (`tw` from the Seifert framing) and `Link::whitehead_double_bbf(&l, positive, tw)` (from the blackboard framing).

For anything else, `LinkBuilder` assembles a diagram port by port; `build()` rejects unconnected or duplicated ports and non-planar wirings.

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

An ordered list of `Edge`s with a `closed: bool` flag — either an *arc* (open path) or a *circle* (closed loop). Constructed via `Path::arc(edges)` / `Path::circ(edges)` / `Path::new(edges, closed)`. Returned as the component type by `Link::comps()` and `Link::seifert_circles()`; free loops become one-element circles. Methods: `is_arc`, `is_circle`, `ends`, `contains`, `connect` (concatenate two arcs sharing an endpoint), `reduce` (canonicalize the edge sequence), `unori_eq` (equality up to orientation).

### `Braid`

A word in the Artin generators of the braid group: a strand count plus a `Vec<BraidGen>`. Each `BraidGen` wraps a non-zero `i8` whose absolute value is the (1-based) strand index and whose sign distinguishes σ_i from σ_i⁻¹. The `i8` storage keeps individual generators small; it caps strand counts at `i8::MAX = 127`, which is well above typical use.

`From<T>` is implemented for `BraidGen` for `T ∈ {i8, i16, i32, i64}` (out-of-range values panic). Correspondingly, `Braid::from(<[T; N]>)` and `Braid::from_iter::<T>` accept any of those element types, so literals like `Braid::from([1, 1, -2])` (default `i32`) work directly.

Methods: `closure() -> Link` (free strands become free loops in the resulting link), `inv`, `reduced()` (collapses adjacent σ σ⁻¹ pairs in a single stack pass; does **not** apply the braid relations), `extend(by)` (appends `by` free strands on the right), `display()` (ASCII rendering), `is_id`, `Mul` for concatenation, `Braid::load(name)` for JSON-backed braid words.

### `InvLink`

An *involutive link*: a `Link` together with an involution on it — an edge bijection `e ↦ e'`, the induced node bijection, and (via the inner `Link`) an optional axis base point. Constructors:

- `InvLink::new(inner, e_map)` — `e_map: IntoIterator<Item = (Edge, Edge)>`. The constructor asserts that `e_map` covers every link edge, has image within the edge set, and is involutive.
- `InvLink::from_symmetric_pd_code(pd_code)` — convenience for *strongly invertible* knots whose PD code uses edges `1..=n` with the standard involution `e ↦ (n+1-e) mod n + 1`.
- `InvLink::load(name)` — reads `<DATA_DIR>/inv_link/<name>.json`. Bundled entries (3_1, 4_1, 5_1, 5_2a/b, …, 7_7a/b) ship in `lib-link/resources/inv_link/` and are copied into the data dir by `scripts/fetch-knot-data.py`.

`with_base_pt(e)` sets the base point; it asserts that `e` is on-axis (`inv_edge(e) == e`). `inv_edge(e)` and `inv_node(x)` look up the involution. Most read-only `Link` methods are delegated, including `base_pt()`.

`mirror` and `conn_sum` / `conn_sum_at` are the equivariant counterparts of the `Link` operations — the connected sum splices along on-axis edges. `InvLink::sym_pretzel(a, b, a)` (all-odd) and `InvLink::whitehead_double(&k, positive, tw)` (even `tw`) are the equivariant constructions; each recovers the strong inversion by reindexing to the standard involution.

### Derived invariants

- `jones_polynomial(&Link) -> LPoly<'q', i32>` — Kauffman-bracket computation summing over all `2^n` resolutions.

## Conventions

- **PD code.** Each crossing is `[a, b, c, d]` ordered counter-clockwise from the lower-left, with `a → c` the incoming under-strand. See:
  - KnotAtlas — [katlas.org/wiki/Planar_Diagrams](https://katlas.org/wiki/Planar_Diagrams)
  - KnotInfo — [knotinfo.org/descriptions/pd_notation.html](https://knotinfo.org/descriptions/pd_notation.html)
- **Edge ids.** `Edge = usize`. Only used to identify endpoint coincidence; ids don't have to be contiguous or start at 0. `InvLink::from_symmetric_pd_code` expects `1..=n_edges` so the involution formula is well-defined.
- **Base point.** Defaults to the minimum edge of the diagram. For `InvLink`, the base point must be fixed by the involution (use `with_base_pt(e)` to set a specific on-axis edge).

## Data directory

The `*::load(name)` constructors read JSON files from a user-data directory — `$YUI_DATA_DIR` if set, otherwise the platform default (`~/Library/Application Support/yui/` on macOS, `${XDG_DATA_HOME:-~/.local/share}/yui/` on Linux, `%APPDATA%\yui\` on Windows). To populate it from the [KnotInfo](https://knotinfo.org/) database:

```bash
python3 scripts/fetch-knot-data.py            # default data dir
python3 scripts/fetch-knot-data.py --out DIR  # or a custom directory
```

The script writes per-knot PD codes to `<DATA_DIR>/links/`, braid words to `<DATA_DIR>/braid/`, and copies this crate's bundled symmetric PD codes (`resources/inv_link/*.json`) to `<DATA_DIR>/inv_link/`.

## Quick example

```rust
use yui_link::Link;

// Trefoil knot from a PD code.
let l = Link::from_pd_code([[1,4,2,5], [3,6,4,1], [5,2,6,3]]);
assert_eq!(l.n_crossings(), 3);
assert_eq!(l.writhe(),     -3);
assert_eq!(l.n_comps(),     1);

// Free loops via `unlink(n)`.
let u = Link::unlink(3);
assert_eq!(u.n_loops(), 3);
assert_eq!(u.n_comps(), 3);
```

## Feature flags

- `test-utils` — exposes the hard-coded `Link::test_data(name)` / `Braid::test_data(name)` / `InvLink::test_data(name)` constructors to downstream crates (behind `#[cfg(test)]` by default).

## License

This library is licensed under the [MIT License](https://opensource.org/licenses/MIT).

---

*This README was generated by [Claude](https://www.anthropic.com/claude).*
