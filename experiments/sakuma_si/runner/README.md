# Sakuma-SI runner

Standalone crate (not a workspace member). Requires the `yui` repo at this
branch (`big-link` feature + the `h_range` pass-through).

    cd experiments/sakuma_si/runner

    # ssi of Wh+(K_a # -K_b), K = a Sakuma 2-bridge knot with two inversions:
    RUST_LOG=debug cargo run --release -- wh 8_21 pos  > out.txt 2> log.txt

    # full sweep over all two-inversion pairs (ssi(a), ssi(b), ssi(a # -b)):
    RUST_LOG=info  cargo run --release            > sweep.csv 2> sweep.log

`ssi_heavy` (used by `wh`) is min-fill + auto(2) cut + window 0..=1.
Edit it in `src/main.rs` for a manual `--cut` if a diagram needs one.
