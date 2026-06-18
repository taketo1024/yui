## yui-kh

[![crates.io](https://img.shields.io/crates/v/yui-kh.svg)](https://crates.io/crates/yui-kh)
[![docs.rs](https://docs.rs/yui-kh/badge.svg)](https://docs.rs/yui-kh)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Khovanov homology and its involutive variant for knots and links, in the [`yui`](https://github.com/taketo1024/yui) workspace. Built on top of [`yui-link`](https://crates.io/crates/yui-link) for diagrams, [`yui-homology`](https://crates.io/crates/yui-homology) for chain complexes, and [`yui-matrix`](https://crates.io/crates/yui-matrix) for linear algebra. Powered by Bar-Natan's tangle-cobordism formalism for fast computations.

## Layout

```text
src/
├── kh/                     — Khovanov homology
│   ├── alg.rs              — `KhAlg`: Frobenius algebra A = R[X]/(X² − hX − t)
│   ├── kh_gen.rs           — `KhGen`: generator (state + tensor label)
│   ├── cube.rs             — `KhCube`: the cube of resolutions (textbook construction)
│   ├── complex.rs          — `KhComplex`: Khovanov chain complex
│   ├── homology.rs         — `KhHomology`: bigraded homology
│   ├── canon_cycle.rs      — Lee/Bar-Natan canonical cycles
│   ├── ss.rs               — `ss_invariant`: family of slice-torus invariants
│   └── ext/                — chain-level extras (crossing change, sl₂ action)
├── khi/                    — Involutive Khovanov homology
│   ├── khi_gen.rs          — `KhIGen = EitherKey<KhGen, KhGen>` (B/Q sides)
│   ├── tau.rs              — chain-level τ for a strongly invertible link
│   ├── complex.rs          — `KhIComplex`: built as `Cone(1 + τ)`
│   ├── homology.rs         — `KhIHomology`
│   └── ssi.rs              — `ssi_invariants`: equivariant Rasmussen pair (s̲, s̄)
├── tng/                    — Bar-Natan's tangle / cobordism category
│   ├── tng.rs              — `Tng`, `TngComp`: Temperley–Lieb diagrams
│   ├── cob.rs              — `Cob`, `LcCob`: cobordism morphisms with dots
│   ├── complex.rs          — `TngComplex`: tangle chain complex
│   ├── elem.rs             — `TngComplexElem`: cycle/element tracker
│   └── builder/
│       ├── builder.rs      — `TngComplexBuilder`: scan-and-cancel algorithm
│       └── sym_builder.rs  — `SymTngBuilder`: equivariant variant for `InvLink`
└── util/                   — small helpers (divisibility, bigraded view)
```

## Key types

### Frobenius algebra

- **`KhAlg<R>`** — the rank-two Frobenius algebra `A = R[X]/(X² − hX − t)` parametrized by `(h, t) ∈ R²`. Specialisations: `(0, 0)` original Khovanov, `(0, 1)` Lee, `(H, 0)` Bar-Natan with `R = F[H]`. Provides `mul`, `comul`, and the dual element `Y = X − h`.
- **`KhAlgGen`** — the two basis elements `{1, X}`.
- **`KhTensor`** — a tensor product of `KhAlgGen`s, packed into a `BitSeq` (length ≤ 64).
- **`KhGen`** — a generator of `CKh`: a resolution `State` plus a `KhTensor` label. Intrinsic gradings via `rel_h_deg()` and `rel_q_deg()` (the complex applies the shift).

### Chain complex and homology

- **`KhComplex<R>`** — Khovanov chain complex of a `Link`, parameterised by `(h, t)`. Constructed via the fast tangle builder (`new`) or the textbook cube (`new_no_simplify`). Carries `KhAlg`, `deg_shift`, a `reduced` flag, and a list of canonical cycles.
- **`KhHomology<R>`** — bigraded homology of `KhComplex`. Computed via [`yui-homology`'s `ChainReducer` + Smith normal form](https://crates.io/crates/yui-homology).
- **`ss_invariant(&Link, c, reduced) -> i32`** — the slice-torus invariant `ss̃_c(K) = 2·d_c(D) + w(D) − r(D) + 1` of a knot.

### Involutive variant (KhI)

- **`KhIComplex<R>`** — the involutive Khovanov complex `CKhI = Cone(CKh →^{Q(1+τ)} Q·CKh)` of a strongly invertible link `InvLink`.
- **`KhIHomology<R>`**, **`KhIGen`** (= `EitherKey<KhGen, KhGen>`), **`KhIGenExt`** trait.
- **`ssi_invariants(&InvLink, c, reduced) -> (i32, i32)`** — the equivariant Rasmussen pair `(s̲, s̄)`.

### Tangle complex (the fast backend)

- **`Tng`**, **`TngComp`** — Temperley–Lieb diagrams; the "objects" of Bar-Natan's category.
- **`Cob`**, **`CobComp`**, **`LcCob<R>`** — dotted cobordism morphisms, modulo `S`, `T`, `4Tu`, and the dotted skein `X² = h·X + t`.
- **`TngComplex<R>`** — chain complex over `Cob_{/l}`; the differential is `LcCob<R>`. Supports `merge` (Bar-Natan tensor product), `deloop` (loop ≅ ∅_X ⊕ ∅_1 in the `(h, t)`-generalised form), and `eliminate` (Gauss-eliminate an invertible edge).
- **`TngComplexBuilder<R>`** — assembles the complex incrementally, choosing crossings by a connection heuristic and deloop+eliminating after each step. The crossing order (`NodeOrder`) and the simplification strategy (`BuildMode`) are tunable via `BuildConfig`.
- **`SymTngBuilder<R>`** — the equivariant variant for `InvLink`: processes axis-symmetric crossings singly, off-axis crossings in `(x, τ·x)` pairs, and maintains a `key_map` so `tau_map()` can return the chain-level `τ` used to build `KhIComplex`.

### Extras (`kh::ext`)

- **`cc_map0`, `cc_map1`** — the two crossing-change chain maps `CKh(D) → CKh(D')` coming from the homology generators of the Hopf-link complex.
- **`KhSl2Map`** — the `e`-operator (raising op of an `sl₂`-action) on `CKh`, with `StringDecomp` for reading off the resulting `R[x]/(x^l)`-summand structure.

## References

The math in this crate follows these papers (cited at the head of each file too):

- **K00** — M. Khovanov, *A categorification of the Jones polynomial*, Duke Math. J. **101** (2000), 359–426. [doi](https://doi.org/10.1215/S0012-7094-00-10131-7), [arXiv](https://arxiv.org/abs/math/9908171).
- **BN02** — D. Bar-Natan, *On Khovanov's categorification of the Jones polynomial*, Algebr. Geom. Topol. **2** (2002), 337–370. [doi](https://doi.org/10.2140/agt.2002.2.337), [arXiv](https://arxiv.org/abs/math/0201043).
- **BN05** — D. Bar-Natan, *Khovanov's homology for tangles and cobordisms*, Geom. Topol. **9** (2005), 1443–1499. [doi](https://doi.org/10.2140/gt.2005.9.1443), [arXiv](https://arxiv.org/abs/math/0410495).
- **K06** — M. Khovanov, *Link homology and Frobenius extensions*, Fund. Math. **190** (2006), 179–190. [doi](https://doi.org/10.4064/fm190-0-6), [arXiv](https://arxiv.org/abs/math/0411447).
- **BN07** — D. Bar-Natan, *Fast Khovanov homology computations*, J. Knot Theory Ramif. **16** (2007), 243–255. [doi](https://doi.org/10.1142/S0218216507005294), [arXiv](https://arxiv.org/abs/math/0606318).
- **SS24** — T. Sano and K. Sato, *A family of slice-torus invariants from the divisibility of Lee classes*, Topol. Appl. **357** (2024), 109059. [doi](https://doi.org/10.1016/j.topol.2024.109059), [arXiv](https://arxiv.org/abs/2211.02494).
- **Sano25** — T. Sano, *Involutive Khovanov homology and equivariant knots*, Algebr. Geom. Topol. **25** (2025), 5059–5111. [doi](https://doi.org/10.2140/agt.2025.25.5059), [arXiv](https://arxiv.org/abs/2404.08568).
- **ISST26** — H. Imori, T. Sano, K. Sato, M. Taniguchi, *Cobordism maps in Khovanov homology and singular instanton homology II* (preprint). [arXiv](https://arxiv.org/abs/2510.09399).
- **Sano26y** — T. Sano, *A y-ification of Khovanov homology* (preprint). [arXiv](https://arxiv.org/abs/2602.17435).

## Quick example

```rust
use yui_link::Link;
use yui_kh::kh::{KhHomology, ss_invariant};

// Khovanov homology of the trefoil over Z.
let l = Link::from_pd_code([[1,4,2,5], [3,6,4,1], [5,2,6,3]]);
let h = KhHomology::<i32>::new(&l, &0, &0, false);
assert_eq!(h[ 0].rank(), 2);                    // Kh⁰  = Z²
assert_eq!(h[-2].tors(), &vec![2]);             // Kh⁻² = Z/2
assert_eq!(h[-3].rank(), 1);                    // Kh⁻³ = Z

// Rasmussen-type s-invariant over (R, c) = (Z, 2).
assert_eq!(ss_invariant(&l, &2, false), -2);
```

## CLI

A command-line front-end ships in the sibling crate [`ykh`](../bin-ykh) — `ykh kh "3_1"`, `ykh khi "3_1"`, etc.

## License

This library is licensed under the [MIT License](https://opensource.org/licenses/MIT).

---

*This README was generated by [Claude](https://www.anthropic.com/claude).*
