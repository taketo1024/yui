#!/bin/sh
# Usage: khi-ssi.sh <name> <symmetric-pd-json> [cut] [h-range]
#
# Computes ssi of the strongly-invertible knot via the cobordism-level cone
# (`ykh khi -s --cob-cone`), logging to profiles/<date>_<name>.{out,log}.
# `cut` defaults to "auto(2)"; pass manual edge-cuts as "e,e,e;e,e,e".
# `h-range` (e.g. "0..1") overrides the ssi window; it must contain 0..=1.
# The log starts with the git branch/commit, so runs are self-documenting.
set -e

name="$1"
pd="$2"
cut="${3:-auto(2)}"
hrange="${4:-..1}"
date=$(date +%Y%m%d)
out="profiles/${date}_${name}.out"
log="profiles/${date}_${name}.log"

mkdir -p profiles
cargo build --release -p ykh

{
    echo "# branch: $(git branch --show-current), commit: $(git rev-parse --short HEAD)"
    echo "# started: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$log"

./target/release/ykh khi "$pd" -c H -s --mode min-fill --cob-cone --cut "$cut" --h-range "$hrange" --log 2 \
    > "$out" 2>> "$log"

echo "done: $(cat "$out")"
