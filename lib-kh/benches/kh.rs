//! Microbenchmarks for `lib-kh`.
//!
//! Use `samply` (or similar) on a single representative run for true
//! profiling; these are *signal* benchmarks for regression detection.
//!
//! Run with:
//! ```
//! cargo bench -p yui-kh -- kh                        # kh-side group
//! cargo bench -p yui-kh -- khi                       # khi-side group
//! cargo bench -p yui-kh -- khi_complex_k18           # one group
//! cargo bench -p yui-kh -- 14n_19265                 # one case
//! ```
//!
//! Heavy cases (44-crossing kh, 60-crossing khi) live as `examples/profile_*`
//! and are intentionally outside the bench harness — too slow for criterion's
//! statistical model.

use criterion::{criterion_group, criterion_main, Criterion};
use yui_core::num::FF2;
use yui_link::{InvLink, Link};
use yui_kh::kh::KhComplex;
use yui_kh::khi::KhIComplex;
use yui_kh::tng::builder::{BuildConfig, SymBuildConfig, NodeOrder};

const NODE_ORDERS: [(&str, NodeOrder); 2] =
    [("loop", NodeOrder::LoopGreedy), ("mincut", NodeOrder::MinCut)];

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// `KhComplex::new` via `TngComplexBuilder` (v2).
fn bench_kh_complex_new(c: &mut Criterion) {
    let mut group = c.benchmark_group("kh_complex_new");

    for name in ["3_1", "5_1", "6_3", "8_19"] {
        let l = Link::test_data(name);
        for (sn, node) in NODE_ORDERS {
            group.bench_function(format!("{name}/{sn}"), |b| {
                b.iter(|| KhComplex::<i32>::new_with_config(&l, &0, &0, false, BuildConfig { node, ..Default::default() }))
            });
        }
    }

    // 14-crossing knot — bigger enough that node-ordering choices matter.
    {
        let l = Link::test_data("14n_19265");
        group.sample_size(10);
        for (sn, node) in NODE_ORDERS {
            group.bench_function(format!("14n_19265/{sn}"), |b| {
                b.iter(|| KhComplex::<i32>::new_with_config(&l, &0, &0, false, BuildConfig { node, ..Default::default() }))
            });
        }
    }

    group.finish();
}

/// `KhIComplex::new` via `SymTngBuilder` over `FF2` (the supported coefficient
/// ring for KhI: `char(R) = 2` is asserted by `new_no_simplify`).
fn bench_khi_complex_new(c: &mut Criterion) {
    let mut group = c.benchmark_group("khi_complex_new");
    group.sample_size(20);

    let zero = FF2::default();

    for name in ["3_1", "4_1", "6_3"] {
        let l = InvLink::test_data(name);
        for (sn, node) in NODE_ORDERS {
            group.bench_function(format!("{name}/{sn}"), |b| {
                b.iter(|| KhIComplex::<FF2>::new_with_config(&l, &zero, &zero, false, SymBuildConfig { node, ..Default::default() }))
            });
        }
    }

    group.finish();
}

/// 18-crossing strongly invertible knot supplied by the user — a realistic
/// scale-up test. Set its own group so we can dial sample size independently.
fn bench_khi_complex_k18(c: &mut Criterion) {
    let mut group = c.benchmark_group("khi_complex_k18");
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(60));

    let pd: &[[u8; 4]] = &[
        [1,27,2,26],[5,16,6,17],[6,32,7,31],[10,27,11,28],[11,1,12,36],
        [13,8,14,9],[14,20,15,19],[17,4,18,5],[18,24,19,23],[21,32,22,33],
        [22,16,23,15],[25,3,26,2],[28,9,29,10],[29,24,30,25],[30,4,31,3],
        [33,20,34,21],[34,8,35,7],[35,13,36,12],
    ];
    let l = InvLink::from_symmetric_pd_code(pd.iter().copied());

    let zero = FF2::default();

    group.bench_function("kh_complex", |b| {
        b.iter(|| KhComplex::<FF2>::new(l.inner(), &zero, &zero, false))
    });

    group.bench_function("khi_complex", |b| {
        b.iter(|| KhIComplex::<FF2>::new(&l, &zero, &zero, false))
    });

    group.finish();
}

/// Selective-deloop variant — the path that exercises `is_dotted_cup`.
fn bench_khi_selective(c: &mut Criterion) {
    use yui_kh::tng::builder::{SymBuildConfig, BuildMode};

    let mut group = c.benchmark_group("khi_selective");
    group.sample_size(10);

    let zero = FF2::default();
    let sel = || SymBuildConfig { mode: BuildMode::Selective, ..Default::default() };

    for name in ["3_1", "4_1", "6_3"] {
        let l = InvLink::test_data(name);
        group.bench_function(name, |b| {
            b.iter(|| KhIComplex::<FF2>::new_with_config(&l, &zero, &zero, false, sel()))
        });
    }

    // K18 — realistic scale where dotted caps actually arise during selective deloop.
    let pd: &[[u8; 4]] = &[
        [1,27,2,26],[5,16,6,17],[6,32,7,31],[10,27,11,28],[11,1,12,36],
        [13,8,14,9],[14,20,15,19],[17,4,18,5],[18,24,19,23],[21,32,22,33],
        [22,16,23,15],[25,3,26,2],[28,9,29,10],[29,24,30,25],[30,4,31,3],
        [33,20,34,21],[34,8,35,7],[35,13,36,12],
    ];
    let l = InvLink::from_symmetric_pd_code(pd.iter().copied());
    group.measurement_time(std::time::Duration::from_secs(60));
    group.bench_function("k18", |b| {
        b.iter(|| KhIComplex::<FF2>::new_with_config(&l, &zero, &zero, false, sel()))
    });

    group.finish();
}

criterion_group!(kh,  bench_kh_complex_new);
criterion_group!(khi, bench_khi_complex_new, bench_khi_complex_k18);
criterion_group!(sel, bench_khi_selective);
criterion_main!(kh, khi, sel);
