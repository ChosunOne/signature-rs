//! SCRATCH profiling driver (never commit): mimics the criterion benches
//! with explicit phase counters. Usage:
//!   prof_concat <d> <m> <threads> <concat|e2e> [iters]
//! Prints lie_rs::lie_series::prof counters divided by iterations.

use ndarray::Array2;
use ordered_float::NotNan;
use rand::{RngExt, SeedableRng, rngs::StdRng};
use signature_rs::LogSignatureBuilder;

type S = NotNan<f64>;

const KERNEL_CALLS_PER_ITER: usize = 64;

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

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 {
        eprintln!("usage: prof_concat <d> <m> <threads> <concat|e2e> [iters]");
        std::process::exit(2);
    }
    let d: usize = args[1].parse().unwrap();
    let m: usize = args[2].parse().unwrap();
    let threads: usize = args[3].parse().unwrap();
    let mode = args[4].clone();
    let iters: usize = args.get(5).map(|v| v.parse().unwrap()).unwrap_or(5);

    let b = LogSignatureBuilder::<u8>::new()
        .with_num_dimensions(d)
        .with_max_degree(m);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("pool");

    let mut wall_ns: u128 = 0;
    match mode.as_str() {
        "concat" => {
            let path = synthetic_path(d, 200);
            let warm = 16;
            let acc = pool
                .install(|| b.build_from_path(&path.slice(ndarray::s![..=warm, ..]).view()));
            let seg = pool.install(|| {
                b.build_from_path(&path.slice(ndarray::s![warm..=warm + 1, ..]).view())
            });
            let segs: Vec<_> = std::iter::repeat(seg.clone())
                .take(KERNEL_CALLS_PER_ITER)
                .collect();
            for _ in 0..iters {
                let mut acc = acc.clone();
                let t0 = std::time::Instant::now();
                pool.install(|| acc.concatenate_batch(&segs));
                wall_ns += t0.elapsed().as_nanos();
                std::hint::black_box(&acc);
            }
        }
        "e2e" => {
            let n = 1000;
            let path = synthetic_path(d, n);
            let view = path.view();
            // one warmup build (plan/table caches), then timed iters
            let _ = pool.install(|| b.build_from_path(&view));
            for _ in 0..iters {
                let t0 = std::time::Instant::now();
                let sig = pool.install(|| b.build_from_path(&view));
                wall_ns += t0.elapsed().as_nanos();
                std::hint::black_box(&sig);
            }
        }
        other => panic!("unknown mode {other}"),
    }

    println!(
        "config d={d} m={m} threads={threads} mode={mode} iters={iters} wall_ms_per_iter={:.3}",
        wall_ns as f64 / iters as f64 / 1e6
    );
}
