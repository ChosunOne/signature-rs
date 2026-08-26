//! Scaling probe for structure-constant table construction in `LieSeries::new`.
//!
//! Motivation: at (d=3, m=8) with D=1318 Lyndon words (1.7M ordered basis
//! pairs), table construction was observed to run single-threaded for
//! >4 minutes, while iisignature performs the analogous offline computation
//! in <1s. At (3,6) with D=196 it takes ~6ms. These tests exist to make the
//! anomaly reproducible and instrumentable (enable the `progress` feature to
//! watch per-row progress).
//!
//! Run:
//!   cargo test -p lie-rs --release --test struct_table_scaling -- --ignored --nocapture
//! Run with progress bars:
//!   cargo test -p lie-rs --release --features progress --test struct_table_scaling -- --ignored --nocapture
//! Run with tracing output (span timings):
//!   cargo test -p lie-rs --release --features tracing --test struct_table_scaling -- --ignored --nocapture

use std::time::Instant;

use lie_rs::LieSeries;

use lyndon_rs::lyndon::{LyndonBasis, Sort};
use num_rational::Ratio;

/// Install a scoped tracing subscriber (for the current thread/test) that
/// prints a `close` event with the elapsed time for every instrumented span.
/// Filter from `RUST_LOG`, defaulting to `lie_rs=debug,commutator_rs=debug`.
fn tracing_timing_guard() -> tracing::subscriber::DefaultGuard {
    use tracing_subscriber::{EnvFilter, fmt::format::FmtSpan, util::SubscriberInitExt};
    tracing_subscriber::fmt()
        .with_span_events(FmtSpan::CLOSE)
        .with_env_filter(
            EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| EnvFilter::new("lie_rs=debug,commutator_rs=debug")),
        )
        .with_test_writer()
        .set_default()
}

fn time_struct_table(d: usize, m: usize, label: &str) {
    let _tracing_guard = tracing_timing_guard();
    let basis = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
    let t = Instant::now();
    let series =
        LieSeries::<u8, Ratio<i128>>::new(basis.clone(), vec![Ratio::from_integer(0); basis.len()]);
    let elapsed = t.elapsed();
    println!(
        "{label} d={d} m={m} D={} pairs={} struct_table={:.3}s",
        series.coefficients.len(),
        series.coefficients.len() * series.coefficients.len(),
        elapsed.as_secs_f64(),
    );
}

/// Small sizes: must complete; guards against regressions in the normal suite.
#[test]
fn struct_table_small_sizes_complete() {
    time_struct_table(3, 4, "small");
    time_struct_table(3, 6, "small");
    time_struct_table(3, 7, "small");
}

/// The case where construction time was observed to explode.
#[test]
#[ignore = "minutes-long scaling probe; run explicitly"]
fn struct_table_d3_m8_completes() {
    time_struct_table(3, 8, "probe");
}

/// Ladder for localizing where per-pair cost blows up.
#[test]
#[ignore = "minutes-long ladder probe; run explicitly"]
fn struct_table_ladder() {
    for m in [4, 5, 6, 7, 8] {
        time_struct_table(3, m, "ladder");
    }
    for (d, m) in [(2, 12), (4, 7), (8, 4)] {
        time_struct_table(d, m, "ladder");
    }
}
