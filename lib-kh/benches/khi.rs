//! KhI (involutive) microbenchmarks for `lib-kh`, over `FF2` (the supported
//! coefficient ring: `char(R) = 2`).
//!
//! Run with:
//! ```
//! cargo bench -p yui-kh --bench khi                  # all khi cases
//! cargo bench -p yui-kh --bench khi -- khi_modes     # mode comparison
//! cargo bench -p yui-kh --bench khi -- k18           # 18-crossing scale test
//! ```
//!
//! Heavy cases (44/60-crossing) live as `examples/profile_*`, outside criterion.

use criterion::{criterion_group, criterion_main, Criterion};
use yui_core::num::FF2;
use yui_link::InvLink;
use yui_kh::kh::KhComplex;
use yui_kh::khi::KhIComplex;
use yui_kh::tng::builder::{SymBuildConfig, BuildMode, NodeOrder};

const NODE_ORDERS: [(&str, NodeOrder); 1] =
    [("mincut", NodeOrder::MinCut)];

const MODES: [(&str, BuildMode); 2] =
    [("greedy", BuildMode::Greedy), ("minfill", BuildMode::MinFill)];

// 18-crossing strongly invertible knot — realistic scale-up.
const K18_PD: &[[u8; 4]] = &[
    [1,27,2,26],[5,16,6,17],[6,32,7,31],[10,27,11,28],[11,1,12,36],
    [13,8,14,9],[14,20,15,19],[17,4,18,5],[18,24,19,23],[21,32,22,33],
    [22,16,23,15],[25,3,26,2],[28,9,29,10],[29,24,30,25],[30,4,31,3],
    [33,20,34,21],[34,8,35,7],[35,13,36,12],
];

/// `KhIComplex::new` across node orders (greedy mode), on small knots.
fn bench_khi_node(c: &mut Criterion) {
    let mut group = c.benchmark_group("khi_node");
    group.sample_size(20);
    let zero = FF2::default();

    for name in ["3_1", "4_1", "6_3"] {
        let l = InvLink::test_data(name);
        for (sn, node) in NODE_ORDERS {
            group.bench_function(format!("{name}/{sn}"), |b| {
                b.iter(|| KhIComplex::<FF2>::new_with_config(&l, &zero, &zero, false, SymBuildConfig { node_order: node, ..Default::default() }))
            });
        }
    }
    group.finish();
}

/// Mode comparison (greedy / selective / min-fill) on small knots — the crossover
/// regime where selective tends to win. (Replaces the old `compare_selective` example.)
fn bench_khi_modes(c: &mut Criterion) {
    let mut group = c.benchmark_group("khi_modes");
    group.sample_size(20);
    let zero = FF2::default();

    // SI knots via their symmetric PD codes (resources/inv_link); test_data covers 3_1/4_1/6_3.
    let knots: Vec<(&str, InvLink)> = vec![
        ("3_1",  InvLink::test_data("3_1")),
        ("4_1",  InvLink::test_data("4_1")),
        ("5_1",  InvLink::from_symmetric_pd_code([[1,7,2,6],[3,9,4,8],[5,1,6,10],[7,3,8,2],[9,5,10,4]])),
        ("5_2",  InvLink::from_symmetric_pd_code([[3,11,4,10],[5,9,6,8],[6,2,7,1],[9,5,10,4],[11,3,12,2],[12,8,1,7]])),
        ("6_3",  InvLink::test_data("6_3")),
        ("7_4",  InvLink::from_symmetric_pd_code([[2,8,3,7],[3,15,4,14],[5,13,6,12],[8,2,9,1],[10,16,11,15],[11,7,12,6],[13,5,14,4],[16,10,1,9]])),
        ("9_46", InvLink::from_symmetric_pd_code([[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]])),
    ];

    for (name, l) in &knots {
        for (mn, mode) in MODES {
            group.bench_function(format!("{name}/{mn}"), |b| {
                b.iter(|| KhIComplex::<FF2>::new_with_config(l, &zero, &zero, false, SymBuildConfig { mode, ..Default::default() }))
            });
        }
    }
    group.finish();
}

/// Matrix cone vs. cobordism-level cone (`--cob-cone`), greedy mode, on small/medium knots.
fn bench_khi_cone(c: &mut Criterion) {
    let mut group = c.benchmark_group("khi_cone");
    group.sample_size(20);
    let zero = FF2::default();

    let knots: Vec<(&str, InvLink)> = vec![
        ("3_1",  InvLink::test_data("3_1")),
        ("4_1",  InvLink::test_data("4_1")),
        ("6_3",  InvLink::test_data("6_3")),
        ("7_4",  InvLink::from_symmetric_pd_code([[2,8,3,7],[3,15,4,14],[5,13,6,12],[8,2,9,1],[10,16,11,15],[11,7,12,6],[13,5,14,4],[16,10,1,9]])),
        ("9_46", InvLink::from_symmetric_pd_code([[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]])),
    ];

    for (name, l) in &knots {
        group.bench_function(*name, |b| {
            b.iter(|| KhIComplex::<FF2>::new_with_config(l, &zero, &zero, false, SymBuildConfig::default()))
        });
    }
    group.finish();
}

/// 18-crossing scale: kh vs khi (greedy), plus the khi mode comparison at scale
/// (the selective path exercises `is_dotted_cup`).
fn bench_khi_k18(c: &mut Criterion) {
    let mut group = c.benchmark_group("k18");
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(60));

    let l = InvLink::from_symmetric_pd_code(K18_PD.iter().copied());
    let zero = FF2::default();

    group.bench_function("kh", |b| {
        b.iter(|| KhComplex::<FF2>::new(l.inner(), &zero, &zero, false))
    });
    for (mn, mode) in MODES {
        group.bench_function(format!("khi/{mn}"), |b| {
            b.iter(|| KhIComplex::<FF2>::new_with_config(&l, &zero, &zero, false, SymBuildConfig { mode, ..Default::default() }))
        });
    }
    group.finish();
}

criterion_group!(khi, bench_khi_node, bench_khi_modes, bench_khi_cone);
criterion_group!(k18, bench_khi_k18);
criterion_main!(khi, k18);
