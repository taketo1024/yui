#!/usr/bin/env python3
"""
Parse the arXiv source of

    C. Lamm, "Symmetric diagrams for all strongly invertible knots up to
    10 crossings", https://arxiv.org/abs/2210.13198

into combinatorial template data (Rust `Template` literals) from which the
symmetric (transvergent) diagrams of all 118 strongly invertible 3-bridge
knots with <= 10 crossings are built (see lib-link/resources/inv_link_lamm/).

Usage:
    python3 scripts/parse_si_template.py [ARXIV_SRC_DIR] > templates_gen.rs

ARXIV_SRC_DIR is the unpacked arXiv source (default: references/arXiv-2210.13198v2
under the repo root); it must contain the template figures t_*.pdf and
Strongly_invertible.tex (Appendix A holds the twist-tuple tables).

Geometry model: each template figure draws the UPPER HALF of a transvergent
diagram as a single stroked path whose moveto-separated pieces follow the knot
arc in traversal order; each gap between consecutive pieces is an under-pass.
Dashed subpaths are twist slots: one end on the axis (x-position = tuple
order), the other end touching the strand that is pulled down to the axis
where |t| crossings are inserted (t > 0: NW-SE strand over = XL).

Canonical crossing frame: under-in at slot 0 (SW), under-out at slot 2 (NE);
the over strand occupies slots 1/3, entering at 1 iff cross(u, v) > 0 in the
math frame (y up), where u = under direction, v = over direction.
Every crossing is XL in this frame.
"""

import re
import sys
import zlib
from pathlib import Path

SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    Path(__file__).resolve().parent.parent / "references" / "arXiv-2210.13198v2"

# ---------------------------------------------------------------- pdf parsing

def content_stream(pdf_path):
    data = pdf_path.read_bytes()
    out = []
    for m in re.finditer(rb"stream\r?\n(.*?)endstream", data, re.S):
        try:
            out.append(zlib.decompress(m.group(1)).decode("latin1"))
        except zlib.error:
            pass
    return "\n".join(out)

def parse_paths(stream):
    """Returns list of (width, dashed, subpaths); subpath = list of sampled points (math frame, y up)."""
    toks = re.findall(r"\[|\]|[-+]?\d*\.?\d+(?:e[-+]?\d+)?|[A-Za-z']+", stream)
    paths, stack, width, dashed = [], [], 1.0, False
    cur, subpaths, pt = [], [], None
    start = None

    def flushsub():
        nonlocal cur
        if cur:
            subpaths.append(cur)
            cur = []

    i = 0
    nums = []
    in_dash_array = False
    dash_nonempty = False
    while i < len(toks):
        t = toks[i]
        i += 1
        if t == "[":
            in_dash_array = True
            dash_nonempty = False
            continue
        if t == "]":
            in_dash_array = False
            continue
        if re.match(r"^[-+]?[\d.]", t):
            if in_dash_array:
                dash_nonempty = dash_nonempty or float(t) != 0
            else:
                nums.append(float(t))
            continue
        # operator
        if t == "w":
            width = nums[-1]
        elif t == "d":
            dashed = dash_nonempty
        elif t == "m":
            flushsub()
            pt = (nums[-2], -nums[-1])
            cur = [pt]
        elif t == "l":
            q = (nums[-2], -nums[-1])
            cur += sample_line(pt, q)
            pt = q
        elif t == "c":
            x1, y1, x2, y2, x3, y3 = nums[-6:]
            p1, p2, p3 = (x1, -y1), (x2, -y2), (x3, -y3)
            cur += sample_cubic(pt, p1, p2, p3)
            pt = p3
        elif t in ("S", "s"):
            flushsub()
            if subpaths:
                paths.append((width, dashed, subpaths))
            subpaths = []
        nums = [] if t not in ("w",) else nums
        if t in ("m", "l", "c", "w", "d"):
            nums = []
    return paths

def sample_line(p, q, n=32):
    return [(p[0] + (q[0] - p[0]) * k / n, p[1] + (q[1] - p[1]) * k / n) for k in range(1, n + 1)]

def sample_cubic(p0, p1, p2, p3, n=64):
    pts = []
    for k in range(1, n + 1):
        t = k / n
        s = 1 - t
        x = s**3 * p0[0] + 3 * s**2 * t * p1[0] + 3 * s * t**2 * p2[0] + t**3 * p3[0]
        y = s**3 * p0[1] + 3 * s**2 * t * p1[1] + 3 * s * t**2 * p2[1] + t**3 * p3[1]
        pts.append((x, y))
    return pts

# ------------------------------------------------------------- reconstruction

def dist(p, q):
    return ((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2) ** 0.5

def extract_template(pdf_path):
    paths = parse_paths(content_stream(pdf_path))

    solid = [p for p in paths if not p[1]]
    dashed = [p for p in paths if p[1]]

    # the axis: solid subpath spanning the largest x-range at the lowest y
    # (math frame: y up, so axis y is the minimum). The curve: thickest stroke.
    curve_pieces = []
    axis_y = None
    for w, _, subs in solid:
        for sp in subs:
            xs = [p[0] for p in sp]
            ys = [p[1] for p in sp]
            if max(ys) - min(ys) < 1.0 and max(xs) - min(xs) > 50:
                axis_y = sum(ys) / len(ys)
            else:
                curve_pieces.append(sp)
    assert axis_y is not None, f"{pdf_path.name}: no axis found"

    # traversal = concatenation of pieces in draw order; ends must lie on the axis
    first, last = curve_pieces[0][0], curve_pieces[-1][-1]
    assert abs(first[1] - axis_y) < 2, f"{pdf_path.name}: curve does not start on axis ({first})"
    assert abs(last[1] - axis_y) < 2, f"{pdf_path.name}: curve does not end on axis ({last})"

    # flatten with piece index; arc position = (piece, sample index)
    # gaps: between piece k end and piece k+1 start
    gaps = []
    for k in range(len(curve_pieces) - 1):
        a, b = curve_pieces[k][-1], curve_pieces[k + 1][0]
        d = dist(a, b)
        assert d < 15, f"{pdf_path.name}: gap {k} too wide ({d:.1f}) — pieces not in traversal order?"
        gaps.append((k, a, b))

    def seg_intersect(p, q, a, b):
        """Does segment p-q intersect segment a-b?"""
        def orient(u, v, w):
            return (v[0] - u[0]) * (w[1] - u[1]) - (v[1] - u[1]) * (w[0] - u[0])
        d1, d2 = orient(a, b, p), orient(a, b, q)
        d3, d4 = orient(p, q, a), orient(p, q, b)
        return (d1 * d2 < 0) and (d3 * d4 < 0)

    # find the over-passage for each gap: the piece whose chord crosses the
    # gap segment (extended slightly for stroke thickness), excluding the
    # under strand's own neighborhood
    def over_at(gap):
        k, a, b = gap
        # extend the gap segment by 20% on both ends
        ext = 0.2
        ea = (a[0] - ext * (b[0] - a[0]), a[1] - ext * (b[1] - a[1]))
        eb = (b[0] + ext * (b[0] - a[0]), b[1] + ext * (b[1] - a[1]))
        hits = []
        for pi, sp in enumerate(curve_pieces):
            for si in range(len(sp) - 1):
                if pi == k and si > len(curve_pieces[k]) - 20:
                    continue
                if pi == k + 1 and si < 20:
                    continue
                if seg_intersect(sp[si], sp[si + 1], ea, eb):
                    hits.append((pi, si))
        assert len(hits) == 1, f"{pdf_path.name}: gap {k}: {len(hits)} over candidates {hits}"
        return hits[0]

    # crossings: one per gap
    crossings = []  # (under_gap_idx, over_piece, over_sample)
    for g in gaps:
        pi, si = over_at(g)
        crossings.append((g, pi, si))

    # traversal events: under-visits at gaps (arc pos = end of piece k),
    # over-visits at (over_piece, over_sample)
    events = []  # (piece, sample_idx, kind, crossing_id)
    for cid, ((k, a, b), pi, si) in enumerate(crossings):
        events.append((k, len(curve_pieces[k]) - 1, "under", cid))
        events.append((pi, si, "over", cid))
    events.sort(key=lambda e: (e[0], e[1]))

    # entry slots
    def direction(pi, si):
        sp = curve_pieces[pi]
        j0, j1 = max(0, si - 3), min(len(sp) - 1, si + 3)
        dx = sp[j1][0] - sp[j0][0]
        dy = sp[j1][1] - sp[j0][1]
        n = (dx * dx + dy * dy) ** 0.5
        return (dx / n, dy / n)

    visits = []
    for pi, si, kind, cid in events:
        (k, a, b), opi, osi = crossings[cid]
        u = ((b[0] - a[0]), (b[1] - a[1]))
        n = (u[0] ** 2 + u[1] ** 2) ** 0.5
        u = (u[0] / n, u[1] / n)
        v = direction(opi, osi)
        cross = u[0] * v[1] - u[1] * v[0]
        over_in_slot = 1 if cross > 0 else 3
        if kind == "under":
            visits.append((cid, 0))
        else:
            visits.append((cid, over_in_slot))

    # twist slots: dashed subpaths; axis end -> tuple order, other end -> arc edge
    slots = []
    for _, _, subs in dashed:
        for sp in subs:
            e0, e1 = sp[0], sp[-1]
            if abs(e0[1] - axis_y) < 2:
                anchor, touch, tangent_from = e0, e1, sp[-6]
            elif abs(e1[1] - axis_y) < 2:
                anchor, touch, tangent_from = e1, e0, sp[5]
            else:
                raise AssertionError(f"{pdf_path.name}: dashed line not anchored on axis")
            # the attached strand: extend the dashed line beyond the touch end along
            # its LOCAL direction and take the first curve chord it crosses (robust
            # when another strand passes nearby sideways)
            n = dist(tangent_from, touch)
            u = ((touch[0] - tangent_from[0]) / n, (touch[1] - tangent_from[1]) / n)
            hit = None
            for reach in range(0, 25):
                probe = (touch[0] + u[0] * reach * 0.25, touch[1] + u[1] * reach * 0.25)
                cands = []
                for pi, spc in enumerate(curve_pieces):
                    d, si = min((dist(p, probe), si) for si, p in enumerate(spc))
                    if d < 1.6:
                        cands.append((d, pi, si))
                if cands:
                    hit = min(cands)[1:]
                    break
            assert hit, f"{pdf_path.name}: dashed line at x={anchor[0]:.1f} points at no strand"
            pi, si = hit
            # arc edge index = number of visit-events strictly before this arc position
            ei = sum(1 for (vp, vs, _, _) in events if (vp, vs) < (pi, si))
            slots.append((anchor[0], ei, (pi, si)))
    slots.sort(key=lambda s: s[0])  # tuple order = left-to-right on the axis

    n_cross = len(crossings)
    assert len(visits) == 2 * n_cross
    return n_cross, visits, [e for _, e, _ in slots], slots

# ----------------------------------------------------------------- tex tuples

def parse_tex():
    tex = (SRC / "Strongly_invertible.tex").read_text()
    a = tex.index("\\section{Appendix A}")
    b = tex.index("\\section{Appendix B}") if "\\section{Appendix B}" in tex else len(tex)
    body = tex[a:b]

    knots = []  # (template, knot, tuple)
    cur_tpl = None
    for line in body.splitlines():
        m = re.search(r"includegraphics\[.*?\]\{t_(\w+)\}", line)
        if m:
            cur_tpl = m.group(1)
            continue
        m = re.match(
            r"\$(?:\\phantom\{1\})?(\d+)_\{(\d+)(?:\\phantom\{0\})?\}\$\s*&=&\((.*)\)", line.strip()
        )
        if m:
            cr, idx, rest = m.group(1), m.group(2), m.group(3)
            ts = [int(x) for x in re.findall(r"-?\d+", rest)]
            knots.append((cur_tpl, f"{cr}_{idx}", ts))
    return knots

# ----------------------------------------------------------------------- main

def main():
    knots = parse_tex()
    by_tpl = {}
    for tpl, k, ts in knots:
        by_tpl.setdefault(tpl, []).append((k, ts))

    names = ["B1", "B2", "B3", "B4", "B5", "B6",
             "C1", "C2a", "C2b", "C3a", "C3b", "C4", "C5", "C6",
             "D1", "E1", "E2", "E3"]

    # Two figures' dashed lines point at strands that do NOT reproduce Lamm's
    # knot tables; the corrected attachments below are fixed empirically by the
    # Kh-homology verification (see notes):
    #   B4: slot 2 belongs on the 8-curl's neck (edge 4; edge 6 equivalent),
    #       not the arch strand (edge 3) its dashed line touches.
    #   E2: the two slots' roles are swapped relative to the anchors' x-order.
    OVERRIDES = {"B4": [5, 4, 2], "E2": [2, 6]}

    print("// generated by parse_lamm.py — do not edit")
    print()
    for name in names:
        n_cross, visits, twist_edges, slots = extract_template(SRC / f"t_{name}.pdf")
        if name in OVERRIDES:
            twist_edges = OVERRIDES[name]
        klist = by_tpl.get(name, [])
        # sanity: tuple lengths match slot count
        for k, ts in klist:
            assert len(ts) == len(twist_edges), f"{name}: {k} tuple len {len(ts)} != {len(twist_edges)} slots"
        vis = ", ".join(f"({c}, {s})" for c, s in visits)
        print(f"pub static {name.upper()}: Template = Template {{")
        print(f"    name: \"{name}\",")
        print(f"    types: &[{', '.join(['XL'] * n_cross)}],")
        print(f"    visits: &[{vis}],")
        print(f"    twist_edges: &[{', '.join(map(str, twist_edges))}],")
        print("};")
        print()

    print("pub fn all() -> Vec<(&'static Template, Vec<KnotSpec>)> {")
    print("    vec![")
    for name in names:
        klist = by_tpl.get(name, [])
        specs = ", ".join(
            f"spec(\"{k}\", &[{', '.join(map(str, ts))}])" for k, ts in klist
        )
        print(f"        (&{name.upper()}, vec![{specs}]),")
    print("    ]")
    print("}")

    total = sum(len(v) for v in by_tpl.values())
    print(f"// total knots: {total}", file=sys.stderr)
    for name in names:
        print(f"//   {name}: {len(by_tpl.get(name, []))} knots", file=sys.stderr)

if __name__ == "__main__":
    main()
