use commutator_rs::Commutator;
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use lie_rs::{LieSeries, IDENTITY_SHIFTS};
use lyndon_rs::LyndonBasis;
use ordered_float::NotNan;

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

/// Grids for the batch-kernel group. `4x10` is included because it is the
/// smallest default-shape grid whose sweep dominates the per-job prologue —
/// the regime where the parallel sweep actually approaches linear scaling
/// (see `bench_commutation_batch`).
fn batch_grid() -> Vec<(usize, usize)> {
    std::env::var("BENCH_BATCH_GRID")
        .unwrap_or_else(|_| "2x12,3x8,4x6,5x5,6x4,8x3,12x2,4x10".to_owned())
        .split(',')
        .map(|tok| {
            let (d, m) = tok.split_once('x').expect("grid entries are `dxm`");
            (
                d.parse().expect("d is an integer"),
                m.parse().expect("m is an integer"),
            )
        })
        .collect()
}

/// Jobs per batch call in `bench_commutation_batch`.
const BATCH_JOBS: usize = 64;

fn dense_coefficients(len: usize) -> Vec<S> {
    (0..len)
        .map(|i| NotNan::new((i % 17 + 1) as f64 / 19.0 - 0.4).unwrap())
        .collect()
}

fn make_series(d: usize, m: usize) -> (usize, LieSeries<u8, S>, LieSeries<u8, S>, LyndonBasis<u8>) {
    let basis = LyndonBasis::<u8>::new(d, lyndon_rs::Sort::Lexicographical);
    let words = basis.generate_basis(m);
    let len = words.len();
    let a = LieSeries::new(words.clone(), dense_coefficients(len));
    let b = LieSeries::new(words, dense_coefficients(len));
    (len, a, b, basis)
}

/// Install a global tracing subscriber so instrumented spans (built with
/// `--features tracing`) report their timings; `RUST_LOG` filters output.
#[cfg(feature = "tracing")]
fn init_tracing() {
    use tracing_subscriber::{EnvFilter, fmt::format::FmtSpan};
    let _ = tracing_subscriber::fmt()
        .with_span_events(FmtSpan::CLOSE)
        .with_thread_ids(true)
        .with_env_filter(
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("lie_rs=debug")),
        )
        .with_writer(std::io::stderr)
        .try_init();
}

#[cfg(not(feature = "tracing"))]
fn init_tracing() {}

fn bench_commutation_kernel(c: &mut Criterion) {
    init_tracing();
    let mut group = c.benchmark_group("commutator_kernel");
    for (d, m) in grid() {
        let (len, a, b, _) = make_series(d, m);
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
                    bencher.iter(|| {
                        // One install per K calls: the closure runs on a pool
                        // worker, so the park/unpark handshake is amortized
                        // away and per-call times reflect the kernel itself.
                        pool.install(|| {
                            for _ in 0..KERNEL_CALLS_PER_ITER {
                                let mut out = vec![S::default(); len];
                                LieSeries::commutator_coefficients(
                                    &a,
                                    &a.coefficients,
                                    &b.coefficients,
                                    &mut out,
                                );
                                std::hint::black_box(&out);
                            }
                        })
                    });
                },
            );
        }
    }
    group.finish();
}

fn bench_commutator(c: &mut Criterion) {
    let mut group = c.benchmark_group("commutator");
    for (d, m) in grid() {
        let (_, a, b, _) = make_series(d, m);
        group.bench_with_input(
            BenchmarkId::new("dxm", format!("{d}x{m}")),
            &(d, m),
            |bencher, _| {
                bencher.iter(|| std::hint::black_box(a.commutator(&b)));
            },
        );
    }
    group.finish();
}

fn bench_lie_series_new(c: &mut Criterion) {
    let mut group = c.benchmark_group("lie_series_new");
    for (d, m) in grid() {
        let (len, _, _, basis) = make_series(d, m);
        group.bench_with_input(
            BenchmarkId::new("dxm", format!("{d}x{m}")),
            &(d, m),
            |bencher, _| {
                bencher.iter(|| LieSeries::new(basis.generate_basis(m), vec![S::default(); len]));
            },
        );
    }
    group.finish();
}

fn bench_commutation_batch(c: &mut Criterion) {
    use lie_rs::lie_series::{GatingCache, KernelJob, commutator_coefficients_batch_with_cache};

    let mut group = c.benchmark_group("commutator_kernel_batch");
    for (d, m) in batch_grid() {
        let (len, a, b, _) = make_series(d, m);
        for t in threads() {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(t)
                .build()
                .expect("failed to build rayon thread pool");
            group.throughput(Throughput::Elements(BATCH_JOBS as u64));
            group.bench_with_input(
                BenchmarkId::new("dxm", format!("{d}x{m}/{t}t")),
                &(d, m),
                |bencher, _| {
                    // Full-support jobs (the fold's big-node pattern): the
                    // kernel's own full-support lists, reused verbatim, so
                    // the batch measures sweep throughput with the prologue
                    // at its production steady state (memoized/full-support).
                    let dense_nz = a.nonzero_coefficient_indices(&a.coefficients);
                    bencher.iter(|| {
                        let mut outs: Vec<Vec<S>> =
                            (0..BATCH_JOBS).map(|_| vec![S::default(); len]).collect();
                        let mut cache = GatingCache::default();
                        let mut jobs: Vec<KernelJob<S>> = outs
                            .iter_mut()
                            .map(|o| KernelJob {
                                a: a.coefficients.as_ptr(),
                                a_len: len,
                                a_nonzero: &dense_nz,
                                b: b.coefficients.as_ptr(),
                                b_len: len,
                                b_nonzero: &dense_nz,
                                result: o.as_mut_ptr(),
                                result_len: len,
                                a_shift: IDENTITY_SHIFTS.as_ptr(),
                                b_shift: IDENTITY_SHIFTS.as_ptr(),
                                r_shift: IDENTITY_SHIFTS.as_ptr(),
                            })
                            .collect();
                        pool.install(|| {
                            commutator_coefficients_batch_with_cache(&a, &mut jobs, &mut cache)
                        });
                        std::hint::black_box(outs.iter().map(|o| o[0]).sum::<S>());
                    });
                },
            );
        }
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_commutation_kernel,
    bench_commutation_batch,
    bench_commutator,
    bench_lie_series_new
);
criterion_main!(benches);
