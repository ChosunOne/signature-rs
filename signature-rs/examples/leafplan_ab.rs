//! A/B driver for the leaf steady-plan cache (`leaf_steady_cache`):
//! measures the steady-state e2e `build_from_path` loop at a given
//! (d, m, threads) — the tournament path whose per-group plan derivation
//! the cache shares.
//!
//! Deterministic (LCG displacement path) — safe to diff across binaries.
//!
//! Usage: leafplan_ab [d] [m] [threads] [e2e_iters] [path_len]

use ndarray::Array2;
use ordered_float::NotNan;
use signature_rs::LogSignatureBuilder;
use std::hint::black_box;
use std::time::Instant;

type F = NotNan<f64>;

fn lcg_path(d: usize, n: usize, seed: u64) -> Array2<F> {
    let mut s = seed;
    let mut lcg = || {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
        ((s >> 33) % 11) as i64 - 5
    };
    Array2::from_shape_fn((n, d), |_| NotNan::new(lcg() as f64).unwrap())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let d: usize = args.get(1).map(|s| s.parse().unwrap()).unwrap_or(3);
    let m: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(8);
    let threads: usize = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(1);
    let iters: usize = args.get(4).map(|s| s.parse().unwrap()).unwrap_or(12);
    let n: usize = args.get(5).map(|s| s.parse().unwrap()).unwrap_or(1000);

    let b = LogSignatureBuilder::<u8>::new()
        .with_num_dimensions(d)
        .with_max_degree(m);

    // Steady-state e2e loop (3 excluded warm-up iterations: the first
    // calls pay the plan-cache cold derivations, mirroring production's
    // first build on a configuration).
    let path = lcg_path(d, n, 0xC0FFEE);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("thread pool");
    for _ in 0..3 {
        black_box(pool.install(|| b.build_from_path(&path.view())));
    }
    let t0 = Instant::now();
    for _ in 0..iters {
        black_box(pool.install(|| b.build_from_path(&path.view())));
    }
    let e2e_ms = t0.elapsed().as_secs_f64() * 1e3 / iters as f64;

    println!(
        "d={d} m={m} t={threads} n={n}: e2e_avg={e2e_ms:.3}ms (iters={iters})"
    );
}
