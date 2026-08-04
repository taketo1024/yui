#!/usr/bin/env python3
"""
Download KnotInfo's database dump and convert it into per-knot JSON files
for yui's runtime data dir.

Usage:
    python3 scripts/fetch-knot-data.py [--out DIR] [--in FILE] [--force] [--keep-xls]

    --out DIR      output directory (default: $YUI_DATA_DIR, else the platform
                   user-data dir for yui).
    --in FILE      use this .xls file instead of downloading.
    --force        overwrite existing JSON files.
    --keep-xls     keep the downloaded .xls file (printed at the end).

Output layout:
    <out>/links/<name>.json     PD code  : JSON array of [a,b,c,d] crossing quads
    <out>/braid/<name>.json     braid    : JSON array of signed generator indices
    <out>/inv_link/<name>.json  involutive-link PD codes bundled in lib-link/resources/

The PD / braid columns are looked up by header name (row 0 of the sheet),
so column reordering upstream does not break the script. The inv_link/
directory is always copied from the repo's `lib-link/resources/inv_link/`
(the Lamm bundle's lamm/A, lamm/B, lamm/two_bridge are flattened into it),
regardless of whether the .xls download happens.

Requires: Python 3.8+ and xlrd<2.0  (`pip install 'xlrd<2.0'`).
"""

import argparse
import json
import os
import shutil
import sys
import tempfile
import urllib.request
from pathlib import Path

try:
    import xlrd
except ImportError:
    sys.exit(
        "error: xlrd is not installed. Install with:\n"
        "    pip install 'xlrd<2.0'\n"
        "\n"
        "If pip refuses with `externally-managed-environment` (PEP 668),\n"
        "use a venv. From the workspace root:\n"
        "    python3 -m venv .venv\n"
        "    .venv/bin/pip install 'xlrd<2.0'\n"
        "    .venv/bin/python scripts/fetch-knot-data.py"
    )


URL = "https://knotinfo.org/knotinfo_data_complete.xls"

COL_NAME  = "name"
COL_PD    = "pd_notation"
COL_BRAID = "braid_notation"


def platform_data_dir() -> Path:
    """Mirror the Rust `directories` crate's ProjectDirs(_, _, "yui").data_dir()."""
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "yui"
    if sys.platform.startswith("win"):
        appdata = os.environ.get("APPDATA") or str(Path.home() / "AppData" / "Roaming")
        return Path(appdata) / "yui" / "data"
    # Linux & other Unix
    xdg = os.environ.get("XDG_DATA_HOME") or str(Path.home() / ".local" / "share")
    return Path(xdg) / "yui"


def resolve_out_dir(cli_value: str | None) -> Path:
    if cli_value:
        return Path(cli_value)
    if env := os.environ.get("YUI_DATA_DIR"):
        return Path(env)
    return platform_data_dir()


def find_col(header, key):
    for i, h in enumerate(header):
        if h == key:
            return i
    raise SystemExit(f"error: column '{key}' not found in header row")


def parse_array(s: str, kind: str, name: str):
    """KnotInfo PD/braid notation is already a JSON-shaped array."""
    s = s.strip()
    if not s:
        return None
    try:
        return json.loads(s)
    except json.JSONDecodeError as e:
        print(f"  skip {name}: malformed {kind}: {e}", file=sys.stderr)
        return None


def copy_bundled_resources(out_dir: Path, force: bool, clean: bool = False) -> None:
    """Copy bundled JSON resources (inv_link/*.json) from the repo into out_dir.

    Copying only ever adds, so a resource deleted from the repo would linger in an
    existing data dir and keep loading. `clean` first removes the directories this
    function owns (the top-level names under resources/), leaving the downloaded
    links/ and braid/ untouched.
    """
    repo_root = Path(__file__).resolve().parent.parent
    src_root = repo_root / "lib-link" / "resources"
    if not src_root.exists():
        print(f"  skip resource copy: {src_root} not found", file=sys.stderr)
        return

    if clean:
        for owned in {p.parts[0] for p in (s.relative_to(src_root) for s in src_root.rglob("*.json"))}:
            target = out_dir / owned
            if target.is_dir():
                shutil.rmtree(target)
                print(f"  cleaned {target}")

    copied = skipped = 0
    for src in src_root.rglob("*.json"):
        rel = src.relative_to(src_root)
        # the Lamm bundle is grouped by construction in the repo (lamm/A, lamm/B,
        # lamm/two_bridge); flatten it into inv_link/ so names load directly.
        if rel.parts[:2] == ("inv_link", "lamm"):
            rel = Path("inv_link") / rel.name
        dst = out_dir / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists() and not force:
            skipped += 1
            continue
        shutil.copyfile(src, dst)
        copied += 1
    print(f"==> bundled resources: copied {copied}, skipped {skipped} (under {src_root})")


def download_xls(dest: Path) -> None:
    print(f"==> fetch knot data from KnotInfo ({URL})")
    with urllib.request.urlopen(URL) as r, open(dest, "wb") as f:
        while chunk := r.read(1 << 16):
            f.write(chunk)


def convert(xls_path: Path, out_dir: Path, force: bool) -> None:
    print(f"==> extract and save data to {out_dir}")
    wb = xlrd.open_workbook(xls_path)
    sh = wb.sheet_by_index(0)
    header = sh.row_values(0)

    i_name  = find_col(header, COL_NAME)
    i_pd    = find_col(header, COL_PD)
    i_braid = find_col(header, COL_BRAID)

    links_dir = out_dir / "links"
    braid_dir = out_dir / "braid"
    links_dir.mkdir(parents=True, exist_ok=True)
    braid_dir.mkdir(parents=True, exist_ok=True)

    n_links = n_braid = skipped = 0

    # row 0 = machine header, row 1 = human header, row 2+ = data
    for r in range(2, sh.nrows):
        row = sh.row_values(r)
        name = (row[i_name] or "").strip()
        if not name:
            skipped += 1
            continue

        pd    = parse_array(row[i_pd],    "pd_notation",    name)
        braid = parse_array(row[i_braid], "braid_notation", name)

        if pd is not None:
            path = links_dir / f"{name}.json"
            if force or not path.exists():
                path.write_text(json.dumps(pd))
                n_links += 1
        if braid is not None:
            path = braid_dir / f"{name}.json"
            if force or not path.exists():
                path.write_text(json.dumps(braid))
                n_braid += 1
        if pd is None and braid is None:
            skipped += 1

    print(f"==> wrote {n_links} links, {n_braid} braids, skipped {skipped}")
    print(f"    output dir: {out_dir}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--out", dest="out", default=None,
                    help="destination directory (overrides $YUI_DATA_DIR and the platform default)")
    ap.add_argument("--in", dest="src", default=None, type=Path,
                    help="use this .xls file instead of downloading")
    ap.add_argument("--force", action="store_true",
                    help="overwrite existing JSON files")
    ap.add_argument("--clean", action="store_true",
                    help="remove bundled-resource dirs (inv_link/) before copying, so resources "
                         "deleted from the repo do not linger")
    ap.add_argument("--keep-xls", action="store_true",
                    help="keep the downloaded .xls (only relevant without --in)")
    args = ap.parse_args()

    out_dir = resolve_out_dir(args.out)

    if args.src:
        if not args.src.exists():
            sys.exit(f"error: input file not found: {args.src}")
        convert(args.src, out_dir, args.force)
        copy_bundled_resources(out_dir, args.force, args.clean)
        return

    # Download to a temp file
    tmp_dir = Path(tempfile.mkdtemp(prefix="knotinfo-"))
    xls = tmp_dir / "knotinfo_data_complete.xls"
    try:
        download_xls(xls)
        convert(xls, out_dir, args.force)
        copy_bundled_resources(out_dir, args.force, args.clean)
    finally:
        if args.keep_xls:
            print(f"    kept .xls at: {xls}")
        else:
            try:
                xls.unlink()
                tmp_dir.rmdir()
            except OSError:
                pass


if __name__ == "__main__":
    main()
