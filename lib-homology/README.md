## yui-homology

[![crates.io](https://img.shields.io/crates/v/yui-homology.svg)](https://crates.io/crates/yui-homology)
[![docs.rs](https://docs.rs/yui-homology/badge.svg)](https://docs.rs/yui-homology)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Chain complexes and homology computation for the [`yui`](https://github.com/taketo1024/yui) workspace. Built on top of [`yui-matrix`](https://crates.io/crates/yui-matrix); used by [`yui-kh`](https://crates.io/crates/yui-kh) (Khovanov homology) and other knot-homology layers.

## Layout

```text
src/
├── conc/
│   ├── add_ind.rs     — `AddInd`: trait for grading types; `isize2`, `isize3`
│   ├── summand.rs     — `Summand<X, R>`: free / finitely-generated R-module summand
│   ├── gr_mod.rs      — `GrMod<I, X, R>`: I-graded R-module
│   ├── complex.rs     — `ChainComplex<I, X, R>` + homology methods
│   ├── chain_map.rs   — `ChainMap`, mapping cone
│   └── generic.rs     — `GenericKey<I>`, `GenericChainComplex<I, R>`
├── algo/
│   ├── homology_calc.rs — `HomologyCalc`: SNF-based homology at a single index
│   └── chain_reducer.rs — `ChainReducer`: pivot/Schur reduction of a complex
└── utils/
    ├── grid.rs        — `Grid<K, V>`: sparse K-indexed table with defaulting `Index`
    ├── format.rs      — `rmod_str` (e.g. "Z² ⊕ Z/2 ⊕ Z/2")
    ├── to_string.rs   — `ToSeqString`, `ToTableString`
    └── tex.rs         — `ToTexSeq`, `ToTexTable`, `tex_rmod_str`
```

## Key types

### Grading

- **`AddInd`** — marker trait for grading types. Implemented for `isize`, `isize2`, `isize3`. Requires `Copy + Zero + Add + Sub + Neg`.

### Modules

- **`Summand<X, R>`** — a free or finitely-generated `R`-module. Stores the raw generators (`indexmap::IndexSet<X>`, which `from_raw_generators` asserts are distinct rather than silently deduplicating), the free rank, torsion coefficients, and a `Trans<R>` recording the basis change from raw generators to the "SNF basis." Use `generator(i)`, `vectorize(z)`, `devectorize(v)` to move between the abstract `Lc<X, R>` form and the coordinate vector.

- **`GrMod<I, X, R>`** — an `I`-graded R-module: a sparse `Grid` of `Summand<X, R>` keyed by grading index. Missing indices are treated as the zero summand. Aliases: `GrMod1<X, R>`, `GrMod2<X, R>`, `GrMod3<X, R>`.

### Chain complexes

- **`ChainComplex<I, X, R>`** — a chain complex: a `GrMod` of summands plus a degree-shift `d_deg: I` and a differential closure `Fn(I, &Lc<X, R>) -> Lc<X, R>`. Aliases: `ChainComplex1<X, R>`, `ChainComplex2<X, R>`, `ChainComplex3<X, R>`. Key methods: `d`, `d_matrix(i)`, `homology()`, `reduced()`, `as_generic()`, plus `check_d_all()` / `check_d_at(i)` / `check_d_for(i, x)` asserting `d² = 0`.

- **`GenericChainComplex<I, R>`** — a chain complex whose generators are anonymous `GenericKey<I>` placeholders, so the complex is fully described by its differential matrices. Built by `from_d_matrices(d_deg, [(i, M_i), ...])`. Useful when you only have matrices, not symbolic generators.

- **`ChainMap<'a, 'c, I, X, Y, R>`** — a chain map between two complexes. Holds references to source / target (`'a`) and a closure (`'c`). Provides `apply`, `make_matrix(i)`, `cone(support, target_based)`, plus `check_all()` / `check_at(i)` / `check_for(i, x)` for verifying `df = fd`.

### Storage

- **`Grid<K, V>`** — sparse `K`-indexed table backed by a `HashMap` plus a stored `V::default()`. `grid[k]` never panics: missing keys return the default. Aliases `Grid1<V>`, `Grid2<V>`, `Grid3<V>` for `isize`, `isize2`, `isize3`.

## Algorithms

- **`HomologyCalc::calculate(d_prev, d_next, with_trans)`** — computes the homology at one index from two consecutive differentials via Smith normal form. Returns `(rank, torsion, Option<Trans>)`. The `with_trans` flag controls whether to record a basis change for cycle representatives.

- **`ChainReducer::reduce(&complex, with_trans)`** — iteratively reduces a chain complex by Gaussian elimination. Each step uses [`yui-matrix`'s pivot finder](https://crates.io/crates/yui-matrix) (Bouillaguet–Delaplace–Voge) to locate a unit pivot and applies a Schur-complement reduction. When `with_trans` is set, the reduction `Trans<R>` is accumulated per index so that cycles in the reduced complex pull back to the original basis. `complex.reduced()` is the typical caller; for matrices-only output use `into_generic_complex()`.

## Computation flow

```text
ChainComplex<I, X, R>
    │ .reduced()            (optional: shrink before SNF)
    ↓
ChainComplex<I, X, R>'      (smaller, same homology)
    │ .homology()
    ↓  per index i:
       d_prev = d_matrix(i - d_deg)
       d_next = d_matrix(i)
       HomologyCalc::calculate(d_prev, d_next, true)
       → Summand built from (rank, tors, Trans)
    ↓
GrMod<I, X, R>              (the homology)
```

## Quick example

```rust
use maplit::hashmap;
use yui_homology::GenericChainComplex1;
use yui_matrix::sparse::SpMat;

// 0 → Z⁴ →d₁ Z⁶ →d₂ Z⁴ → 0  (the 2-sphere as a triangulated complex)
let c = GenericChainComplex1::<i32>::from_d_matrices(-1, hashmap! {
    0 => SpMat::from_row_major((0, 4), []),
    1 => SpMat::from_row_major((4, 6), [
        -1,-1, 0,-1, 0, 0,
         1, 0,-1, 0,-1, 0,
         0, 1, 1, 0, 0,-1,
         0, 0, 0, 1, 1, 1,
    ]),
    2 => SpMat::from_row_major((6, 4), [
         1, 1, 0, 0,
        -1, 0, 1, 0,
         1, 0, 0, 1,
         0,-1,-1, 0,
         0, 1, 0,-1,
         0, 0, 1, 1,
    ]),
});

let h = c.homology();
assert_eq!(h[0].rank(), 1);   // H₀(S²) = Z
assert_eq!(h[1].rank(), 0);   // H₁(S²) = 0
assert_eq!(h[2].rank(), 1);   // H₂(S²) = Z
```

## Feature flags

- `multithread` (default) — parallelizes `d_matrix` column construction via [`rayon`].

## License

This library is licensed under the [MIT License](https://opensource.org/licenses/MIT).

---

*This README was generated by [Claude](https://www.anthropic.com/claude).*
