# Symmetric diagrams of strongly invertible knots

PD codes of transvergent diagrams for the strongly invertible knots up to 10
crossings, built from

> C. Lamm, *Symmetric diagrams for all strongly invertible knots up to 10
> crossings*, https://arxiv.org/abs/2210.13198

A *transvergent* diagram has its axis of π-rotation in the diagram plane; that
rotation is the strong inversion. Each file is numbered so its least edge lies
on the axis, which is what `InvLink::from_symmetric_pd_code` expects — it
recovers τ by walking the traversal both ways from the base point. Beyond that
the numbering is arbitrary: it is not canonical, and neither is the traversal
direction. A canonical choice, if one is wanted, is the least
`Link::reindexed(e, 1)` over the two on-axis edges — edge 1 stays on the axis,
so the result still loads.

The bundle is **complete**: every prime strongly invertible knot with at most 10
crossings is here, with both equivalence classes of strong inversions wherever
there are two. Lamm's own counts close the arithmetic — 118 knots of bridge
number 3 (§5), and among two-bridge knots, all of which are strongly invertible,
75 + 10 + 6 = 91 non-torus (§11) plus the 4 torus knots — so 213 knots in 368
diagrams. Every knot in this range has bridge number 2 or 3, so nothing else can
occur.

## Layout and loading

`scripts/fetch-knot-data.py` flattens the three directories into
`<DATA_DIR>/inv_link/`, so the load name is the file stem:
`InvLink::load("9_46")`, `InvLink::load("9_46a")`.

| dir | files | source |
|---|---|---|
| `A/` | 118 | Appendix A, templates `B₁`–`B₆`, `C₁`–`C₆`, `D₁`, `E₁`–`E₃` |
| `B/` | 64 | Appendix B, template `X` — doubly transvergent, 32 knots × 2 axes |
| `two_bridge/` | 186 | §4, two-bridge plat diagrams — 91 knots × 2 classes, + 4 torus knots |

`A/` is exactly Lamm's 118 three-bridge knots. The four (2,n) torus knots
`3_1, 5_1, 7_1, 9_1` sit in `two_bridge/` — they come from the same plat
construction, but have only one equivalence class of strong inversions (Lamm's
two-class statement excludes torus knots), so they are unlabeled.

`10_35` comes from the fractions 49/20 and 49/22 at 10 crossings, `10_22` from
49/34 and 49/36 at **14** — Lamm's condition for `c_t(K) = c(K)` fails for it.
Kh cannot tell these two knots apart, so their fractions are named by the
oriented isometry instead (`scripts/name_49.py`).

The 32 knots of `B/` also appear in `A/` with fewer crossings (`10_141`: 10 vs
14, `10_144`: 10 vs 12, `9_46`: 9 vs 10) — prefer `A/` when one inversion is
enough; which axis it realizes is not recorded.

## Naming

`<knot>.json` leaves the equivalence class unlabeled; `<knot>a/b.json` are the
two classes.

- **Template `X`** — `a` is the horizontal axis, `b` the vertical one of the
  same doubly transvergent diagram. `8_18` (period 4) and `8_19` (torus) are
  excluded: their axes give equivalent inversions, so only their `A/` file
  exists.
- **Two-bridge** — for `b(p,q)`, `a` comes from the smaller even representative
  of `{q, q⁻¹ mod p}`, `b` from the other. When `q² ≡ −1` the knot is
  amphichiral and `b` is the mirror of `a`, since mirroring swaps the classes.
  The six knots with `q² ≡ 1` (`7_4, 7_7, 9_10, 9_17, 9_23, 9_31`) get one
  class from the fraction; their `b` uses Lamm's vertical construction — the
  plat of the palindromic odd-length positive continued fraction, axis through
  the middle box.

These labels are **not** aligned with Lobb–Watson (arXiv:1908.00082 §6.5).

## Verification

Every diagram is checked against **KnotInfo's own diagram of the knot it is named
after**, and is stored mirrored when that is what matches — so the whole bundle
carries KnotInfo's chirality. All of this is reproducible from `scripts/`; each
script's header comment gives its usage.

| check | tool | result |
|---|---|---|
| strongly invertible knot | `is_knot() && is_strongly_invertible()` | 368/368 |
| Kh over ℤ vs KnotInfo's diagram | `ykh kh` | 368/368 |
| same knot, same chirality | SnapPy complement isometry | 362/362 |
| crossing count `2·c(U) + Σ\|tᵢ\|` | template vs diagram | 118/118 |

Kh alone cannot certify this bundle. Over ℤ it fails to separate 13 pairs of
knots ≤ 10 crossings, and 24 knots share their mirror's homology
(`scripts/kh_collisions.py`):

```
10_40=10_103   10_71=10_104   10_59=10_106   10_81=10_109   8_8=10_129
10_137=10_155  8_16=10_156    10_22=10_35    10_25=10_56    10_41=10_94
10_43=10_91    10_60=10_86    10_73=10_83
```

The isometry check closes the gap. Two knot exteriors admit an
orientation-*preserving* isometry exactly when the knots agree including
chirality (Gordon–Luecke), and SnapPy reports the orientation as the
determinant of the isometry's cusp map. Reversing the knot's direction sends
`(μ, λ) ↦ (−μ, −λ)`, determinant `+1`, so the test ignores direction — Lamm's
diagrams carry none anyway. It needs hyperbolicity, so it skips the 6 torus
knots `3_1, 5_1, 7_1, 8_19, 9_1, 10_124`; each of those has a Kh table unique
among all 498 tables of knots ≤ 10 crossings and their mirrors, so Kh settles
them. Together the two checks are complete over this range.

## How the diagrams were built

A **template** shows the upper half of a transvergent diagram plus integer
markings on the axis. Each marking carries a **twist number** `tᵢ` and a dotted
line to the strand it acts on: that strand is pulled to the axis, where `|tᵢ|`
crossings go between it and its mirror image (`tᵢ > 0` puts the `NW–SE` strand
over, i.e. `XL` here — checked on all 36 knots of `C₁`). The lower half is the
mirror of the upper, preserving each crossing's `NodeType` and sending its
ports `s ↦ 3−s`. So `c = 2·c(U) + Σ|tᵢ|`, matching Lamm's Table 1. For template
`X` a **quarter** plays the same role, with markings on both axes: the tuple
before `|` runs along the horizontal axis left to right, the one after it down
the vertical axis, where the sign convention is rotated 90° (`tᵢ > 0` = `XR`).

Two things the figures do not fix. **Which side the strand enters a twist
region** — the two choices differ by a π-rotation of that region, only one is
planar, and `LinkBuilder::build` picks it; the knot never changes. **The order
of markings sharing one strand** — in template `X` four markings sit on one
strand but belong to different axes, and different orders can both be planar
and give different knots.

Where a figure under-determines the reading, the generator enumerated the valid
variants — planar, right crossing count, valid inversions — and kept the one
whose Khovanov homology matched the *named* knot; the name is ground truth, not
the dotted lines. Two readings disagree with Lamm's own tables:

- **`B₄`, marking 2** — the dotted line touches the arch strand, but only
  attaching at the 8-curl's neck gives `10_56 = B₄(−1,2,−1)` and
  `10_57 = B₄(−1,2,2)`. The four `B₄` knots with `t₂ = 0` verify either way.
- **`E₂`** — the two markings' roles are swapped relative to their anchors'
  left-to-right order (calibrated on `10_97 = E₂(2,1)`).

**`9_35` in template `X` is probably an erratum.** The printed
`(1, 0, −1 | 0, 0, 2)` gives `8_3` under every reading — all `2⁶` entry sides,
both orders on the shared strand, every sign flip. Searching all weight-4
tuples (the weight is forced by the crossing count) finds `9_35` only at a
different support, `(1, 0, 0 | −1, 0, 2)`, which no ambiguity above can
explain. The other 33 entries verify under one fixed convention, so this bundle
uses the corrected tuple.

Caveat from the paper: for `9_35`, `9_40`, `9_48` and `10_75` it is unclear
whether the two axes belong to different equivalence classes (their symmetry
groups may exceed `D₂`).
