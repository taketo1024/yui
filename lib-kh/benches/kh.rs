//! Microbenchmarks for `lib-kh`.
//!
//! These are *signal* benchmarks — they confirm regression on small/medium
//! knots, they do not predict behavior on the 60+ crossing targets. Use
//! `samply` (or similar) on a single representative run for true profiling.
//!
//! Run with:
//! ```
//! cargo bench -p yui-kh                              # all
//! cargo bench -p yui-kh -- kh_complex_new            # filter by group
//! cargo bench -p yui-kh -- kh_complex_new/v2/8_19    # one case
//! ```

use criterion::{criterion_group, criterion_main, Criterion};
use yui_core::num::FF2;
use yui_link::{InvLink, Link};
use yui_kh::kh::{KhComplex, KhHomology};
use yui_kh::khi::KhIComplex;

/// `KhComplex::new` via `TngComplexBuilder` (v2).
fn bench_kh_complex_new(c: &mut Criterion) {
    let mut group = c.benchmark_group("kh_complex_new");

    for name in ["3_1", "5_1", "6_3", "8_19"] {
        let l = Link::test_data(name);
        group.bench_function(name, |b| {
            b.iter(|| KhComplex::<i32>::new(&l, &0, &0, false))
        });
    }

    group.finish();
}

/// Full Khovanov homology (complex + reduction + SNF).
fn bench_kh_homology(c: &mut Criterion) {
    let mut group = c.benchmark_group("kh_homology");
    group.sample_size(20);

    for name in ["3_1", "5_1", "6_3", "8_19"] {
        let l = Link::test_data(name);
        group.bench_function(name, |b| {
            b.iter(|| KhHomology::<i32>::new(&l, &0, &0, false))
        });
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
        group.bench_function(name, |b| {
            b.iter(|| KhIComplex::<FF2>::new(&l, &zero, &zero, false))
        });
    }

    group.finish();
}

/// 18-crossing strongly invertible knot supplied by the user — a realistic
/// scale-up test. Set its own group so we can dial sample size independently.
fn bench_khi_complex_k18(c: &mut Criterion) {
    let mut group = c.benchmark_group("khi_complex_k18");
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(60));

    let pd: &[[usize; 4]] = &[
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

criterion_group!(
    benches,
    bench_kh_complex_new,
    bench_kh_homology,
    bench_khi_complex_new,
    bench_khi_complex_k18,
);
criterion_main!(benches);
