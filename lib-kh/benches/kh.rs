//! Kh (non-involutive) microbenchmarks for `lib-kh`.
//!
//! Run with:
//! ```
//! cargo bench -p yui-kh --bench kh                   # all kh cases
//! cargo bench -p yui-kh --bench kh -- kh_modes       # mode comparison
//! cargo bench -p yui-kh --bench kh -- 14n_19265      # one case
//! ```
//!
//! Heavy cases (44/60-crossing) live as `examples/profile_*` and are
//! intentionally outside criterion. KhI benchmarks live in `benches/khi.rs`.

use criterion::{criterion_group, criterion_main, Criterion};
use yui_link::Link;
use yui_kh::kh::KhComplex;
use yui_kh::tng::builder::{BuildConfig, BuildMode, NodeOrder};

const NODE_ORDERS: [(&str, NodeOrder); 1] =
    [("mincut", NodeOrder::MinCut)];

const MODES: [(&str, BuildMode); 3] =
    [("greedy", BuildMode::Greedy), ("selective", BuildMode::Selective), ("minfill", BuildMode::MinFill)];

const KNOTS: [&str; 8] = ["3_1", "4_1", "5_1", "5_2", "6_2", "6_3", "7_3", "8_19"];

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn build(l: &Link, cfg: BuildConfig) -> KhComplex<i32> {
    KhComplex::<i32>::new_with_config(l, &0, &0, false, cfg)
}

/// `KhComplex::new` across node orders (greedy mode).
fn bench_kh_node(c: &mut Criterion) {
    let mut group = c.benchmark_group("kh_node");

    for name in KNOTS {
        let l = Link::test_data(name);
        for (sn, node) in NODE_ORDERS {
            group.bench_function(format!("{name}/{sn}"), |b| {
                b.iter(|| build(&l, BuildConfig { node, ..Default::default() }))
            });
        }
    }

    // 14-crossing knot — big enough that node-ordering choices matter.
    let l = Link::test_data("14n_19265");
    group.sample_size(10);
    for (sn, node) in NODE_ORDERS {
        group.bench_function(format!("14n_19265/{sn}"), |b| {
            b.iter(|| build(&l, BuildConfig { node, ..Default::default() }))
        });
    }
    group.finish();
}

/// Mode comparison (greedy / selective / min-fill), min-cut node order.
fn bench_kh_modes(c: &mut Criterion) {
    let mut group = c.benchmark_group("kh_modes");

    for name in KNOTS {
        let l = Link::test_data(name);
        for (mn, mode) in MODES {
            group.bench_function(format!("{name}/{mn}"), |b| {
                b.iter(|| build(&l, BuildConfig { mode, ..Default::default() }))
            });
        }
    }

    let l = Link::test_data("14n_19265");
    group.sample_size(10);
    for (mn, mode) in MODES {
        group.bench_function(format!("14n_19265/{mn}"), |b| {
            b.iter(|| build(&l, BuildConfig { mode, ..Default::default() }))
        });
    }
    group.finish();
}

criterion_group!(kh, bench_kh_node, bench_kh_modes);
criterion_main!(kh);
