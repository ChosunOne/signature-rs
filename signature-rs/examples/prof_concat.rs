//! Parallel concatenation and end-to-end build demo + profiling entry point.
//!
//! Demonstrates the two hot paths of `signature-rs` on a synthetic random
//! walk:
//!   * `LogSignature::concatenate_batch` — folding a batch of segment
//!     signatures into one result (the parallel tournament reduction with
//!     the SIMD-across-folds cohort engine),
//!   * `LogSignatureBuilder::build_from_path` — the end-to-end path →
//!     log-signature computation.
//!
//! Profiling: build with frame pointers and capture with `cargo flamegraph`
//! (requires perf + inferno):
//!
//!   CARGO_PROFILE_BENCH_DEBUG=true RUSTFLAGS="-Cforce-frame-pointers=yes" \
//!   cargo flamegraph --release -p signature-rs --example prof_concat \
//!     -o /tmp/flame-concat.svg --deterministic \
//!     -c "record -F 997 --call-graph fp -g" -- 3 8 8 concat
//!
//! Frame-pointer unwinding is required: the default dwarf call graphs fail
//! to unwind the rayon worker stacks in this workload.
//!
//! Usage: prof_concat <d> <m> <threads> <concat|e2e> [iters]

use ndarray::Array2;
use ordered_float::NotNan;
use rand::{RngExt, SeedableRng, rngs::StdRng};
use signature_rs::LogSignatureBuilder;

type S = NotNan<f64>;

/// Segment signatures concatenated per timed `concat` iteration.
const BATCH_SEGMENTS: usize = 64;

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

fn die(msg: &str) -> ! {
    eprintln!("error: {msg}");
    std::process::exit(2);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 {
        die("usage: prof_concat <d> <m> <threads> <concat|e2e> [iters]\n       example: prof_concat 3 8 8 concat");
    }
    let d: usize = args[1].parse().unwrap_or_else(|_| die("d must be a usize"));
    let m: usize = args[2].parse().unwrap_or_else(|_| die("m must be a usize"));
    let threads: usize = args[3].parse().unwrap_or_else(|_| die("threads must be a usize"));
    let mode = args[4].clone();
    if mode != "concat" && mode != "e2e" {
        die("mode must be 'concat' or 'e2e'");
    }
    let iters: usize = args
        .get(5)
        .map(|v| v.parse().unwrap_or_else(|_| die("iters must be a usize")))
        .unwrap_or(5);

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
                .take(BATCH_SEGMENTS)
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
        other => die(&format!(
            "unknown mode {other} (expected 'concat' or 'e2e')"
        )),
    }

    println!(
        "config d={d} m={m} threads={threads} mode={mode} iters={iters} wall_ms_per_iter={:.3}",
        wall_ns as f64 / iters as f64 / 1e6
    );
}
