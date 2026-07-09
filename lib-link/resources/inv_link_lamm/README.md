# Symmetric diagrams of strongly invertible knots (Lamm's templates)

PD codes of transvergent (axis-symmetric) diagrams for strongly invertible
knots up to 10 crossings, built from

> C. Lamm, *Symmetric diagrams for all strongly invertible knots up to 10
> crossings*, https://arxiv.org/abs/2210.13198

All files load via `InvLink::from_symmetric_pd_code` (edges numbered so the
strong inversion is the standard `e ↦ (n+1-e) mod n + 1`). Every diagram is
verified by comparing its Khovanov homology over `ℤ` with the KnotInfo
reference diagram of the named knot, and is **mirror-normalized to KnotInfo's
chirality**. Note this differs from the convention of our earlier
`inv_link/5_2a.json` etc. (read from Lobb–Watson, arXiv:1908.00082 §6.5),
which match these files only after mirroring; e.g. KnotInfo's `9_46` has
writhe −3, the mirror of the positive-writhe `9_46` used in our ssi papers.

Generation: `scripts/parse_si_template.py` parses the arXiv source (vector
figures `t_*.pdf` + the tuple tables in the tex) into combinatorial template
data; a `LinkBuilder`-based generator (session tooling, not yet committed)
assembles and verifies the diagrams.

## Reconstruction model

A *template* is the upper half `U` of a transvergent diagram: an open arc with
its two endpoints on the axis (the two τ-fixed points of the knot). Each twist
slot is a **finger**: a marked point of `U` is pulled down to the axis, where
`|tᵢ|` crossings are inserted between the finger and its mirror image
(`tᵢ > 0`: the `NW–SE` strand passes over, i.e. `XL` in lib-link's convention;
verified on all 36 knots of template `C₁`). The lower half is the mirror copy
of `U` — reflection plus crossing switch, which preserves each crossing's
`NodeType` while mirroring its ports (`s ↦ 3-s`). Crossing count:
`2·c(U) + Σ|tᵢ|`, matching the diagram sizes in Lamm's Table 1.

For the doubly transvergent template `X` (Appendix B) the same model applies
twice: the quarter `Q` (an open arc with ONE crossing, one end on each axis)
is reflected into four copies, with twist chains on both axes; x-chains use
the convention above, y-chains the 90°-rotated one (`t > 0` = `XR`). Each
diagram carries two strong inversions `τx`, `τy`.

## Wiring freedoms and identification-driven selection

Two aspects of the reconstruction are not determined by the figures alone.

**Entry side of a finger.** Splicing a twist chain into an arc-edge, the
strand may enter the chain at its left or right corner. For an isolated finger
the two choices give the same tangle (they differ by a π-rotation of the twist
region about a vertical axis — a rigid motion), and in context only one choice
is planar; `LinkBuilder::build`'s planarity validation selects it. This
freedom never changes the resulting knot.

**Nesting order on shared edges.** When several fingers attach to the *same*
arc-edge (template `X`: four slots on edge `e0`, going to *different* axes),
the diagram depends on the order in which they are spliced along the edge.
Unlike the entry side, different orders can **both be planar and give
different knots**. The drawn dashed lines fix one order, but the figure is a
single picture serving all tuples of the template; for a particular tuple the
isotoped diagram may effectively slide one finger's base past another's.

Where the figure under-determines the wiring, the generator enumerates the
finitely many structurally valid variants — planar, correct crossing count,
valid strong inversion(s) — and accepts the one whose Khovanov homology equals
that of the *named* knot. The paper's knot name is the ground truth, not the
dashed lines. Safeguards against accepting a wrong-but-Kh-equal diagram:

- the crossing count `2c + Σ|tᵢ|` and the validity of the involution(s) are
  pinned before Kh is consulted;
- Kh over `ℤ` is injective on knots ≤ 10 crossings up to mirror, with the
  single exception of the pair `10_22` / `10_35`;
- for the doubly transvergent diagrams, Lamm's Prop. 11.3 gives slack: for a
  prime hyperbolic knot with cyclic period 2, *any* doubly transvergent
  diagram realizes the two inversion classes on its two axes — so
  knot-correctness plus double-transvergence already yields the intended
  `(K, τ₁), (K, τ₂)` pair, independent of which valid wiring was chosen.

Cases where the drawn attachment provably disagrees with the paper's own knot
tables (found by this method, each documented in the generator):

- **`B₄`, slot 2**: the middle dashed line touches the arch strand (arc-edge
  3), but only attaching at the 8-curl's neck (edge 4; edge 6 is equivalent)
  reproduces `10_56 = (-1,2,-1)` and `10_57 = (-1,2,2)`. The four `B₄` knots
  with `t₂ = 0` are insensitive to this and verify either way.
- **`E₂`**: the two slots' roles are swapped relative to their axis anchors'
  left-to-right order (calibrated on `10_97 = (2,1)`; `9_29 = (1,1)` is
  order-insensitive).

## The `9_35` entry of Fig. 16 (probable erratum)

The printed tuple `9_35 = (1, 0, -1 | 0, 0, 2)` does not produce `9_35` under
*any* wiring variant: all `2⁶` entry-side combinations, both nesting orders of
the shared-edge fingers, and every sign flip of the printed support yield a
12-crossing diagram of `8_3`. An exhaustive search over all weight-4 tuples
(the weight is fixed by the crossing count) finds `9_35` only with a
**different support**: `(1, 0, 0 | -1, 0, 2)` — slots `x₁, y₁, y₃` instead of
`x₁, x₃, y₃`. A support difference cannot be explained by wiring freedom, so
either the paper's entry is a typo, or Lamm's reading of template `X`'s
attachments differs from the vector figure in a way only this entry exposes.
The other 33 entries — including all others using the same slots — verify
under one fixed convention. This bundle uses the corrected tuple.

Caveat inherited from the paper: for `9_35`, `9_40`, `9_48`, `10_75` the
inequivalence of the two axes' inversion classes is stated as unverified
(their symmetry groups may exceed `D₂`); the corrected `9_35` diagram carries
the same caveat.

## Naming

- `<knot>.json` — a single transvergent diagram, inversion class unlabeled
  (Appendix A templates give one strong inversion per knot; torus knots have
  only one class).
- `<knot>a.json` / `<knot>b.json` — the two strong-inversion classes:
  - 3-bridge knots with cyclic period 2 (template `X`): `a` = `τx`, `b` = `τy`
    of the same doubly transvergent diagram. Exceptions: `8_18` (period 4) and
    `8_19` (torus), whose two axes give equivalent inversions (§11.4) — these
    keep only their plain Appendix-A files;
  - 2-bridge knots `b(p,q)`: `a` = the class built from the smaller even
    representative among `{q, q⁻¹ mod p}`, `b` = the other. For amphichiral
    2-bridge knots (`q² ≡ -1`), `b` is the mirror image of `a`'s diagram
    (mirroring swaps the two classes). For the six knots with `q² ≡ 1`
    (`7_4, 7_7, 9_10, 9_17, 9_23, 9_31`) the single fraction gives only one
    class; `b` comes from Lamm's vertical (Fig. 6) construction — the 4-plat
    of the palindromic odd-length positive continued fraction, drawn
    vertically with the axis through the middle box.

These `a`/`b` labels follow the conventions above; they are **not** aligned
with the labels of Lobb–Watson §6.5 (used in `resources/inv_link/`).

Not yet included: `10_22` and `10_35` (2-bridge, both `b(49,·)`) — the pair is
Kh-identical, so assigning the right name to each fraction needs KnotInfo's
2-bridge notation.
