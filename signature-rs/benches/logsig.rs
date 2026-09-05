use std::time::Duration;

use criterion::{
    BatchSize, BenchmarkId, Criterion, SamplingMode, Throughput, criterion_group, criterion_main,
};
use ndarray::Array2;
use ordered_float::NotNan;
use rand::{RngExt, SeedableRng, rngs::StdRng};
use signature_rs::LogSignatureBuilder;

type S = NotNan<f64>;

/// Kernel calls per timed iteration: amortizes the once-per-install
/// park/unpark handshake (`pool.install` from a non-member thread) below
/// measurement noise. Reported times are per call via `Throughput`.
const KERNEL_CALLS_PER_ITER: usize = 64;

fn grid() -> Vec<(usize, usize)> {
    std::env::var("BENCH_GRID")
        .unwrap_or_else(|_| "2x12,3x8,4x6,5x5,6x4,8x3,12x2".to_owned())
        .split(',')
        .map(|tok| {
            let (d, m) = tok.split_once('x').expect("grid entries are `dxm`");
            (d.parse().unwrap(), m.parse().unwrap())
        })
        .collect()
}

fn e2e_grid() -> Vec<(usize, usize)> {
    std::env::var("BENCH_E2E_GRID")
        .unwrap_or_else(|_| "2x12,3x8,4x6,5x5,6x4,8x3,12x2".to_owned())
        .split(',')
        .map(|tok| {
            let (d, m) = tok.split_once('x').expect("grid entries are `dxm`");
            (d.parse().unwrap(), m.parse().unwrap())
        })
        .collect()
}

fn threads() -> Vec<usize> {
    std::env::var("BENCH_THREADS")
        .unwrap_or_else(|_| "1,2,4,8,16,32".to_owned())
        .split(',')
        .map(|t| {
            t.parse()
                .expect("thread counts are comma-separated integers")
        })
        .collect()
}

fn e2e_path_length() -> usize {
    std::env::var("BENCH_N")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(1000)
}

fn synthetic_path(d: usize, n: usize) -> Array2<S> {
    let mut rng = StdRng::seed_from_u64(0x9E37_79B9_7F4A_7C15);
    let mut path = Array2::<f64>::zeros((n + 1, d));
    for k in 1..=n {
        for c in 0..d {
            path[[k, c]] = path[[k - 1, c]] + rng.random_range(-1.0..1.0);
        }
    }
    path.mapv(|v| NotNan::new(v).expect("path contains NaN"))
}

fn builder(d: usize, m: usize) -> LogSignatureBuilder<u8> {
    LogSignatureBuilder::<u8>::new()
        .with_num_dimensions(d)
        .with_max_degree(m)
}

const WARM_SEGMENTS: usize = 16;

fn bench_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("build");
    for (d, m) in grid() {
        let b = builder(d, m);
        group.bench_with_input(
            BenchmarkId::new("dxm", format!("{d}x{m}")),
            &(d, m),
            |bencher, _| {
                bencher.iter(|| std::hint::black_box(b.build::<S>()));
            },
        );
    }
    group.finish();
}

fn bench_e2e(c: &mut Criterion) {
    let mut group = c.benchmark_group("logsig_e2e");
    group.sampling_mode(SamplingMode::Flat);
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(5));
    let n = e2e_path_length();
    for (d, m) in e2e_grid() {
        let b = builder(d, m);
        let path = synthetic_path(d, n);
        for t in threads() {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(t)
                .build()
                .expect("failed to build rayon thread pool");
            group.bench_with_input(
                BenchmarkId::new("dxmxn", format!("{d}x{m}x{n}/{t}t")),
                &(d, m, n),
                |bencher, _| {
                    let view = path.view();
                    // Per-iteration install: a full path build makes thousands
                    // of pool dispatches, so the single install handshake is
                    // negligible.
                    bencher
                        .iter(|| pool.install(|| std::hint::black_box(b.build_from_path(&view))));
                },
            );
        }
    }
    group.finish();
}

fn bench_concat_single(c: &mut Criterion) {
    let mut group = c.benchmark_group("concat_single");
    for (d, m) in grid() {
        let b = builder(d, m);
        let path = synthetic_path(d, 200);
        let acc = b.build_from_path(&path.slice(ndarray::s![..=WARM_SEGMENTS, ..]).view());
        let seg = b.build_from_path(
            &path
                .slice(ndarray::s![WARM_SEGMENTS..=WARM_SEGMENTS + 1, ..])
                .view(),
        );
        // The batch's displacement set: the same segment folded K times
        // (the steady-state concat workload).
        let segs: Vec<_> = std::iter::repeat(seg.clone())
            .take(KERNEL_CALLS_PER_ITER)
            .collect();
        for t in threads() {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(t)
                .build()
                .expect("failed to build rayon thread pool");
            group.throughput(Throughput::Elements(KERNEL_CALLS_PER_ITER as u64));
            group.bench_with_input(
                BenchmarkId::new("dxm", format!("{d}x{m}/{t}t")),
                &(d, m),
                |bencher, _| {
                    bencher.iter_batched(
                        || acc.clone(),
                        |mut acc| {
                            // The production fold shape: one continuous
                            // batch dispatch over all K displacements
                            // (same batch driver as build_from_path).
                            pool.install(|| {
                                acc.concatenate_batch(&segs);
                            });
                            std::hint::black_box(acc)
                        },
                        BatchSize::SmallInput,
                    );
                },
            );
        }
    }
    group.finish();
}

criterion_group!(benches, bench_build, bench_e2e, bench_concat_single);
criterion_main!(benches);
