# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.
It is written and maintained by Claude, not by hand.

## Overview

`yui` (結) is a Rust workspace for homology computations, focused on knot homology theories — specifically Khovanov homology. The codebase is mathematical in nature: the types model algebraic structures (rings, modules, chain complexes), and algorithms operate over them.

## Commands

```bash
# Build
cargo build
cargo build --release

# Test all crates
cargo test

# Test a single crate
cargo test -p yui-core
cargo test -p yui-matrix
cargo test -p yui-homology
cargo test -p yui-link
cargo test -p yui-kh

# Run a specific test
cargo test -p yui-homology grid

# Run the CLI binary
cargo run -p ykh -- --help
cargo run --release -p ykh -- kh "3_1"

# Doctests are disabled in every crate (`[lib] doctest = false`), so the README
# examples are not compiled by `cargo test`. Run them explicitly:
cargo test --workspace --doc

# Multithread is on by default; disable it:
cargo test -p yui-matrix --no-default-features

# Logging: the ykh CLI gates env_logger by its own flag (--log 1/2/3 = info/debug/trace);
# RUST_LOG is IGNORED by ykh, but respected by the profile examples (from_default_env).
cargo run --release -p ykh -- kh "3_1" --log 2
RUST_LOG=info cargo run --release --example profile_kh_br44
```

## Workspace layout and dependency order

```
lib-core      — algebraic trait hierarchy, concrete numeric/poly/lc types, algorithms, extensions
    ↓
lib-matrix    — sparse/dense matrix types and decompositions (PLUQ, SNF)
    ↓
lib-homology  — chain complexes, chain maps, homology computation
    ↓
lib-link      — knot/link data structures (planar diagrams, braids)
    ↓
lib-kh        — Khovanov homology (KhComplex, KhHomology, KhI invariant, ss)
    ↓
bin-ykh       — CLI frontend (clap-based, commands: kh, khi, ckh, ckhi, cc, sl2)
```

## Key abstractions

### `lib-core`

Directory layout:

```
src/
├── abst/    — algebraic trait hierarchy (Ring, EucRing, Field, ...)
├── conc/    — concrete types (Ratio, FF, Sign, Poly, Lc, BitSeq, ...)
├── ext/     — extension traits (IteratorExt, RangeExt, IntoDigits, TeX, ...)
├── algo/    — algorithms (TopSort, UnionFind, KeyedUnionFind, rep_comb)
└── util/    — formatting, data-dir resolution, thread-safe counter
```

- **Trait hierarchy** (`src/abst/`):
  - `MathType` (basic.rs) — math object identified by `math_symbol()`; bounds: `Default + Eq + Clone + Send + Sync + Display + Debug + 'static`.
  - `IndexType` (basic.rs) — blanket-impl shape trait for keys/indices: `MathType` bounds + `Hash + Ord`.
  - `AddMon → AddGrp`, `Mon → Ring → EucRing → Field`, plus `RMod` for modules over a ring.
  - Integer-specific traits in `conc/num/int.rs`: `IntType` (signed, totally-ordered `EucRing`) and `IntOps` helper.
- **Concrete numeric types** (`src/conc/num/`): `Ratio<T>` (rationals), `FF2` (𝔽₂), `FF<P>` (𝔽_p), `QuadInt<I, D>` (with aliases `GaussInt`, `EisenInt`), integer impls for `i32/i64/i128/BigInt`, and `Sign` (the multiplicative `{±1}`, with `GetSign`) — reachable as `yui_core::num::Sign`, not at the crate root.
- **Polynomials** (`src/conc/poly/`): `Poly<X, R>`, `MultiVar`, typed variable wrappers `Var`, `Var2`, `Var3`.
- **`Lc<X, R>`** (`src/conc/lc/`): a formal linear combination over keys `X: LcKey` with coefficients in `R: Ring`. This is the fundamental "element" type for chain modules. `LcKey: MathType + IndexType` — any hashable, displayable math object can be a key. Wrap a non-math type with `AsKey<T>` if needed.
- **Misc concrete types** (`src/conc/misc/`): `BitSeq<I>` (packed sequence of bits over a word `I`, length up to `I::BITS`; aliases `BitSeq8` … `BitSeq128`), `BitMap`, `U256`.
- **Algorithms** (`src/algo/`):
  - `TopSort` — incremental builder for topological sort. Backed by `petgraph::algo::toposort`.
  - `UnionFind` / `KeyedUnionFind` — disjoint-set structure. Inner uses `petgraph::unionfind::UnionFind` (path compression + union by rank). `KeyedUnionFind` wraps with `indexmap::IndexSet<X>` for hashable keys.
  - `rep_comb` — stars-and-bars combinatorics helper.
- **Extension traits** (`src/ext/`):
  - `CloneAnd::clone_and(f)` — clone, mutate, return.
  - `IntoDigits` — decompose into base-10 digits.
  - `DivRound` — rounded integer division.
  - `IteratorExt::range()` — `Option<RangeInclusive<_>>` of an iterator's min/max.
  - `RangeExt` — shift range endpoints.
  - `TeX` — LaTeX rendering.

### `lib-matrix`

See `lib-matrix/README.md` (embedded in the crate-level docs) for the user-facing overview. Highlights for working in this crate:

Directory layout:

```
src/
├── perm.rs       — `Perm`: permutation of 0..n
├── dense/
│   ├── mat.rs    — `Mat<R>` (over `nalgebra::DMatrix`) + `MatTrait`
│   ├── pluq.rs   — dense PLUQ + solver
│   ├── snf.rs    — Smith normal form
│   └── lll.rs    — LLL / Hermite normal form
└── sparse/
    ├── sp_mat.rs — `SpMat<R>` (over `nalgebra_sparse::CscMatrix`)
    ├── sp_vec.rs — `SpVec<R>`
    ├── trans.rs  — `Trans<R>`: forward/backward basis-change tracking
    ├── pivot.rs  — pivot finder (Bouillaguet–Delaplace–Voge)
    ├── pluq.rs   — sparse PLUQ + incremental solver
    ├── schur.rs  — Schur-complement reduction
    ├── snf.rs    — sparse Smith normal form
    └── triang.rs — triangular solves
```

- **`SpMat<R>`** — sparse matrix (CSC). Core type for all differential maps.
  - `iter()` walks stored triplets and may yield explicit zeros; `iter_nz()` filters them.
  - Block ops: `block_split` / `block_combine` (2×2 quadrants), `h_split` / `v_split` (1-D), `h_stack` / `v_stack` (1-D combine).
  - Permutation: `permute(p, q)`, `permute_rows`, `permute_cols`, `permute_and_split(p, q, r)`.
- **`SpVec<R>`** — sparse column vector (CSC with one column).
- **`Mat<R>`** — dense matrix.
- **`Perm`** — permutation of `0..n`, stored as `Either<usize, Vec<usize>>` so `Perm::id(n)` is zero-cost.
  - Composition is right-to-left math convention: `(p * q)(i) = p(q(i))`.
  - `apply_to(y)` is the left action `result[p(i)] = y[i]`; `apply_inv_to(y)` is the inverse-left action without allocating an inverse perm.
  - `forward_indices(n, prefix)` builds a perm that pulls a given prefix to the front.
- **`Trans<R>`** — composable forward/backward transformation; tracks basis changes through reductions.
- **`MatTrait`** — minimal shape interface (`shape`, `n_rows`, `n_cols`, `is_square`); both `Mat` and `SpMat` implement it.

Algorithms (intentionally narrow — targeted at the homology pipeline):

| Algorithm | Dense | Sparse |
|---|---|---|
| Heuristic pivot finder | — | `sparse::pivot::find_pivots` |
| Schur-complement reduction | — | `sparse::schur::Schur` |
| Triangular solve | (used internally by dense PLUQ) | `sparse::triang` |
| PLUQ decomposition over a ring | `dense::pluq::pluq` | `sparse::pluq::pluq` / `pre_pluq` |
| Linear solve over a field | `dense::pluq::solve_pluq` | `sparse::pluq::solve_pluq` (+ incremental `solve_pluq_incr`) |
| Smith normal form | `dense::snf` | `sparse::snf::sp_snf` |
| LLL / Hermite normal form | `dense::lll::{lll, lll_hnf}` | — |

### `lib-homology`

See `lib-homology/README.md` for the user-facing overview. Highlights for working in this crate:

Directory layout:

```
src/
├── conc/
│   ├── add_ind.rs     — `AddInd`: grading-index trait; `isize2`, `isize3` types
│   ├── summand.rs     — `Summand<X, R>`: free / fg R-module summand
│   ├── gr_mod.rs      — `GrMod<I, X, R>`: I-graded R-module (aliases `GrMod1`/`GrMod2`/`GrMod3`)
│   ├── complex.rs     — `ChainComplex<I, X, R>` + homology methods (aliases `ChainComplex1`/`2`/`3`)
│   ├── chain_map.rs   — `ChainMap<'a, 'c, ...>`, mapping cone
│   └── generic.rs     — `GenericKey<I>`, `GenericChainComplex<I, R>`, `GenericGrMod<I, R>`
├── algo/
│   ├── homology_calc.rs — `HomologyCalc`: SNF-based homology at one index
│   └── chain_reducer.rs — `ChainReducer`: pivot/Schur reduction of a complex
└── utils/
    ├── grid.rs        — `Grid<K, V>`: sparse K-indexed table with defaulting `Index`
    ├── format.rs      — `rmod_str` (e.g. "Z² ⊕ Z/2 ⊕ Z/2")
    ├── to_string.rs   — `ToSeqString`, `ToTableString`
    ├── matrix.rs      — `make_matrix`: the matrix of a map between two generator sets
    └── tex.rs         — `ToTexSeq`, `ToTexTable`, `tex_rmod_str`
```

Naming convention (Option B): the unsuffixed name is the **generic struct** parameterized over the grading `I`; numeric suffixes are aliases for fixed gradings.

- `ChainComplex<I, X, R>` is the struct; `ChainComplex1<X, R> = ChainComplex<isize, X, R>` (singly-graded), `ChainComplex2<X, R>`, `ChainComplex3<X, R>`. Same shape for `GrMod`, `Grid`, `GenericChainComplex`, `GenericGrMod`. (`yui-core`'s `Poly` uses the older convention `PolyBase` / `Poly` / `Poly2` / ... — not yet aligned.)

Key types:

- **`AddInd`**: `IndexType + Copy + Zero + Add + Sub + Neg`. Implementors: `isize`, `isize2`, `isize3` only — `usize`-based variants were dropped (not closed under subtraction).
- **`Summand<X, R>`**: stores raw generators (`indexmap::IndexSet<X, ahash::RandomState>`), free rank, torsion coefficients, and a `Trans<R>` recording the basis change from raw to "SNF basis." `vectorize`/`devectorize` move between `Lc<X, R>` and coordinate `SpVec<R>`.
- **`ChainComplex<I, X, R>`**: stores a `GrMod<I, X, R>` of summands, a degree-shift `d_deg: I`, and an `Arc<dyn Fn(I, &Lc<X, R>) -> Lc<X, R> + Send + Sync>` differential closure. Optionally caches per-index `SpMat<R>`s (populated by `reduced()` / `from_d_matrices`) so `d_matrix(i)` becomes a clone instead of re-vectorizing each generator.
- **`ChainMap<'a, 'c, I, X, Y, R>`**: holds `&'a` source/target complexes and a closure of lifetime `'c`. `cone(support, target_based)` produces a new `ChainComplex` and so requires `'c: 'static`.
- **`GenericChainComplex<I, R>`** (alias for `ChainComplex<I, GenericKey<I>, R>`): built by `from_d_matrices(d_deg, [(i, M_i), ...])` when only the matrices are known.

Computation flow:

- `ChainComplex::homology()` iterates over the support, calling `HomologyCalc::calculate(d_{i-d_deg}, d_i, true)` per index — runs SNF, returns `(rank, torsion, Trans)`. Each resulting summand merges the original `c.trans()` with the homology `Trans` so cycles can be pulled back.
- `ChainComplex::reduced()` → `ChainReducer::reduce(complex, true)`: repeatedly finds a unit pivot (`yui-matrix`'s `find_pivots`), applies the Schur complement to shrink the complex, and accumulates the basis change. The reducer's matrices are passed to the resulting `ChainComplex` via `with_d_matrices`, so the reduced complex's `d_matrix(i)` is a cache hit.
- `ChainReducer::into_generic_complex()` exits to a `GenericChainComplex` when symbolic generators are no longer needed.

### `lib-kh`

- **`kh/`**: `KhComplex<R>` — the Khovanov chain complex for a `Link`. Built via `TngComplexBuilder` (v2, default) or `KhCube` (v1, simpler). `KhHomology` = homology of `KhComplex`. The bigraded homology is a `Grid2<KhComplexSummand<R>>`.
  - `KhAlg<R>` encodes the Frobenius algebra structure (parameters `h`, `t`).
  - `ss` module computes the Rasmussen s-invariant.
  - `kh/ext/sl2.rs`: sl(2) representation extensions.
- **`khi/`**: Khovanov homology of involutive links (`KhIComplex`, `KhIHomology`, `ssi` invariant).
- **`TngComplexBuilder`** (`tng/builder/builder.rs`): the efficient construction algorithm — processes crossings one by one, maintaining a `TngComplex` (cobordism-based) and auto-delooping/eliminating where possible. Tuned via `BuildConfig { node_order: NodeOrder, strategy: Strategy, cut: CutOption, h_range, q_range, max_elim_cost, no_full_deloop }`:
  - `NodeOrder` (crossing order): `MinCut` (default; bounds boundary cutwidth, wins on wide knots), `Given` (PD order, no reordering).
  - `Strategy` (simplification): `Greedy` (default; deloop every circle, eliminate immediately), `MinFill` (whole-degree deloop + Markowitz elimination), `NoElim` (deloop, no eliminate), `None` (raw merge; finalize still deloops to a valid complex).
- **`SymTngBuilder`** (`tng/builder/sym_builder.rs`): the equivariant variant for `InvLink` — processes axis-symmetric crossings singly, off-axis crossings in `(x, τx)` pairs, and maintains a `key_map` so `tau_map()` yields the chain-level `τ` for `KhIComplex`. Its `SymBuildConfig` adds `preprocess` and equivariant chunking controls. Forwards to its `inner: TngComplexBuilder` via `delegate!`.

### `lib-link`

See `lib-link/README.md` for the user-facing overview. Highlights for working in this crate:

Directory layout:

```
src/
├── link/
│   ├── link.rs       — `Link`, `Edge`, `State`: the type, accessors, `comps`, traversal, `reindexed`
│   ├── link_ops.rs   — operations *on* a link: `mirror`, `conn_sum`, `cc_at`, `resolve_*`, `seifert_*`
│   ├── construct.rs  — links from *patterns*: `twist_knot`, `pretzel`, `cable2`, `whitehead_double`
│   ├── pd_code.rs    — `PDCodeX`, `from_pd_code`, `pd_code`, `Link::load`
│   ├── builder.rs    — `LinkBuilder`, `Port`, `LinkError`
│   ├── node.rs       — `Node`, `NodeType` {XL, XR, V, H}, `Slot` {SW, SE, NE, NW}
│   └── path.rs       — `Path`: arc or circle of edges
├── inv_link/
│   ├── inv_link.rs   — `InvLink`: involutive link, plus `mirror` and `conn_sum`
│   └── construct.rs  — `sym_pretzel`, `whitehead_double` (equivariant)
├── braid/
│   ├── braid.rs      — `Braid`
│   └── braid_gen.rs  — `BraidGen` (signed Artin generator, stored as `i8`)
├── misc/jones.rs     — `jones_polynomial`
└── test_data.rs      — hard-coded PD codes / braid words (`#[cfg(any(test, feature = "test-utils"))]`)
resources/inv_link/   — bundled symmetric PD codes (3_1 … 7_7b) for `InvLink::load`
```

Key types:

- **`Link`** — planar diagram. Fields: `nodes: Vec<Node>`, `loops: Vec<Edge>` (free loops — closed components without crossings), `base_pt: Option<Edge>` (defaults to the minimum edge in the diagram). No internal `edges` field — `edges()` recomputes a sorted+deduped `Vec<Edge>` on demand; `n_edges()` is O(1) (`nodes.len() * 2 + loops.len()`, relying on the "each node-edge appears exactly twice" invariant). Primary constructor: `Link::new(nodes, loops)`; `from_nodes(nodes)` is a wrapper for the no-loops case. Predefined shapes: `empty`, `unknot`, `unlink(n)`. `with_base_pt(e)` is a consuming builder (asserts `e` is a real edge).

- **`Node`** — vertex of a diagram: `NodeType ∈ {XL, XR, V, H}` (crossings vs. smoothings), `incoming: Option<(Slot, Slot)>`, and `[Edge; 4]` for incident edges. `is_crossing`, `is_resolved`, `resolve(bit)`, `mirror`, `sign()`.
  - **`Slot`** — one of a node's four ends, `{SW, SE, NE, NW}` counter-clockwise from the lower left. `Slot::ALL`, `index()`, `shift(k)`, `From<usize>`. A link-level port is `(usize, Slot)` — node index plus slot.
  - **Orientation** is the pair of slots the two strands *enter* by, sorted; `None` when the node is not coherently oriented. Valid exactly when the two lie on different strands (`Node::orientable`), the strand pairing being `NodeType::paired_slot` (`XL/XR: SW<->NE`, `V: SW<->NW`, `H: SW<->SE`). This replaced the old `NodeOri` enum, which could not express a `V`/`H` node whose two strands run opposite ways (the diagonal pairs). A crossing's orientation survives exactly one of its two smoothings; `resolve` clears it for the other.

- **`Path`** — `Vec<Edge>` + `closed: bool`. Used as the `comps()` element type; free loops become one-element circle paths.

- **`LinkBuilder`** — assembles a diagram from crossings and `Port = (NodeIndex, slot)` pairs: `add_crossing`, `add_node`, `add_loop`, `add_v_twist` / `add_h_twist` (a chain of `k ≥ 1` crossings, returning its four corner ports in CCW order SW/SE/NE/NW), `add_link` (absorb an existing link, returning its node-index → vertex map), `connect` / `disconnect`, `edge_at`. `build()` returns `Result<Link, LinkError>` — it rejects unconnected ports, duplicated ports and non-planar wirings. Every construction in `construct.rs` goes through it.

- **`BraidGen`** — signed Artin generator, single `i8` field (so individual generators fit in one byte and strand counts cap at 127). `From<T>` is implemented for `T ∈ {i8, i16, i32, i64}` via a `macro_rules!` (out-of-range values panic via `i8::try_from`). The `pub(super) fn from_raw(value: i8)` is a no-check inner constructor used by `Braid::from_iter` to skip the nonzero assertion (the same panic happens later in `closure()` anyway).

- **`Braid`** — `(strands: usize, elements: Vec<BraidGen>)`. `Braid::from_iter<T>` and `Braid::from(<[T; N]>)` are macro-generated for `T ∈ {i8, i16, i32, i64}`, so default-`i32` literals work. `closure() -> Link` allows free strands — strands untouched by any crossing become free loops in the result. `reduced()` does single-pass stack-based cancellation of σ σ⁻¹ pairs only (does **not** apply the braid relations). `extend(by)` appends `by` free strands on the right.

- **`InvLink`** — `Link` + `e_map: HashMap<Edge, Edge>` + `x_map: HashMap<Node, Node>`. Base point lives in the inner `Link` (delegated). `InvLink::new(inner, e_map)` is the only constructor: it asserts no free loops (unsupported), that `e_map` is an involution covering every edge, and that it carries each node to a node of the same type.
  - **Involutive vs strongly invertible**: an `InvLink` is only an *involutive* link. `is_strongly_invertible()` (τ reverses the orientation) and `is_2periodic()` (τ preserves it) decide the two cases, both from the orientation via `preserves_dir_at` — never from index arithmetic, which is degenerate when the diagram has 2 edges. `on_axis_edges()` lists the τ-fixed edges. `conn_sum` and `whitehead_double` assert knot + strongly invertible.
  - **`si_knot_from(inner)`** (public) builds the involution from a diagram *based on its axis*: τ reverses the traversal, so walking both ways from the base point pairs each edge with its image — `τ(e_{k+i}) = e_{k-i}`. No search, no relabelling, and the old `e ↦ (n+1-e) mod n + 1` formula is gone. `from_symmetric_pd_code` is just this applied to `Link::from_pd_code` (the least edge is on-axis by that convention).
  - **Determinism**: `Link::conn_sum_at` and `whitehead_double_impl` carry a base point through `LinkBuilder`'s renumbering (via `edge_at` before `build`), so `InvLink::conn_sum` (splice `self`'s *other* on-axis edge to `other`'s base point) and `whitehead_double` hand `si_knot_from` a diagram already based on the axis. Base points are no longer normalised to edge 1.
  - `InvLink::load(name)` reads `<DATA_DIR>/inv_link/<name>.json`; bundled JSON files at `lib-link/resources/inv_link/` are copied into the data dir by `scripts/fetch-knot-data.py`.
  - `inv_edge(e)` / `inv_node(x)` look up the involution (renamed from `inv_e` / `inv_x`).
  - `with_base_pt(e)` asserts `e` is on-axis (`inv_edge(e) == e`), then delegates to the inner `Link::with_base_pt`.
  - `reversed()` on `Node` / `Link` / `InvLink` reverses the orientation — each strand enters by the slot it used to leave. Crossing signs and the writhe are unchanged, and τ is untouched.

Conventions specific to this crate:

- `Link::edges()` returns `Vec<Edge>` (sorted+deduped), not an iterator — the old `impl Iterator<Item = &Edge>` was changed when the `edges` field was dropped.
- **`link_ops.rs` vs `construct.rs`**: an *operation* takes a link you already have and derives something from it (so it is a method, `&self`); a *construction* builds a diagram out of a pattern (so it is an associated function — `Link::whitehead_double(&l, …)`, not `l.whitehead_double(…)`), taking the companion link as an argument when it needs one. `InvLink` mirrors the same split. Keep `link.rs` / `inv_link.rs` to the type, its accessors and traversal.
- `Link::pretzel(a, b, c)` numbers its edges from the top spanning arc, so its least edge is on-axis — which is what `InvLink::sym_pretzel` needs.
- `serde` / `serde_json` are hard dependencies (not feature-gated): all three `*::load` paths parse JSON files from the data dir.

## Style conventions

(Personal coding and communication preferences are kept outside this repository; below is what is specific to this workspace.)

- Keep the code readable line-by-line against an explicit picture of the math: picture-true names, no construction-history artifacts. Prefer keeping `lib-kh` types close to the Bar-Natan / Khovanov vocabulary; when speed would obscure the math, surface the trade-off explicitly instead of silently optimizing.
- Macro-generated methods (e.g. `delegate!`-produced ones) are invisible to the LSP — fall back to text search on the underlying delegated target.
- **Dependent projects**: `yui-hfk` and `yui-kr` depend on this workspace — test them after API changes here.
- Knot-name references: **KnotInfo** is the only database site (there is no "LinkInfo"); **KnotAtlas** writes `K11a84` (no underscore).

## Commit conventions

- A commit message is a **single subject line** — no body, no `Co-Authored-By` trailer (this overrides any global/default trailer rule).

## Versioning & releases

- All crates currently share one version. Patch releases are **independent** — bump only the crate
  that changed, e.g. `yui-core` to `0.5.1`, and publish just that one. Minor and major bumps are
  **lockstep**: move every crate together, so the versions re-converge.
- Versions are literal in each crate's manifest, not inherited from `[workspace.package]`, precisely
  so a single crate can diverge for a patch. `authors` / `license` / `edition` / `rust-version` /
  `repository` are inherited.
- Publish order is forced by the dependency chain:
  `yui-core → yui-matrix → yui-homology → yui-link → yui-kh → ykh`, each after the previous
  appears in the registry index.
- `Cargo.lock` is tracked, so the `ykh` binary has a reproducible build. Library consumers ignore it.

## Benchmarking & profiling

- Profile with `RUST_LOG=debug` and capture stderr to a log — INFO hides the per-degree deloop/eliminate detail, and the summary alone hides scheduling non-determinism. lib-kh, lib-matrix and lib-homology all log, so a narrow filter like `yui_kh::tng::builder=debug` drops the matrix-reduction lines; narrow only when the output is genuinely too noisy.
- Benches are only compiled by `cargo bench` / `cargo build --benches`; renames can silently break them — run `cargo build --benches` after API changes.
- Criterion baselines live in `target/criterion/`, which `cargo clean` wipes — copy them aside first if you need bench deltas across a clean.
- If `cargo bench` hangs at 0% CPU, it's Criterion's broken gnuplot probe — put a fast-failing `gnuplot` stub on `PATH`.
- When cross-checking against another implementation, feed both sides the **identical valid PD code** over the same ring. A `braid.closure()` re-export is not a valid PD code, and comparing against one has produced spurious disagreements.

## Generic patterns to know

- Nearly all types are generic over a ring `R: Ring, for<'x> &'x R: RingOps<R>`. The `for<'x>` HRTB on the reference impl is required everywhere.
- The `I: AddInd` type parameter is the grading — `isize` for singly-graded, `isize2` for bigraded.
- `auto_impl_ops` is used heavily to derive arithmetic operator impls from a single `+` or `*` impl.
- `multithread` is a Cargo feature (default on) that switches `iter()` to `par_iter()` via `cfg_if!` blocks inside `ChainComplex::d_matrix`.
- `hashmap!` literals come from the `maplit` crate (`use maplit::hashmap;`).
- Iterator helpers: prefer `IteratorExt::range(self) -> Option<RangeInclusive<_>>` over the older free-function form for computing min..=max of an iterator.
- Logging uses the `log` crate; the ykh CLI inits `env_logger` from its `--log 1/2/3` flag (`RUST_LOG` is ignored there); the profile examples init `from_default_env` (`RUST_LOG` works).
