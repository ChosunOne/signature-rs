//! # lie-rs
//!
//! A Rust library for Lie series and Baker-Campbell-Hausdorff computations.
//!
//! ## Quick Start
//!
//! ```rust
//! use lie_rs::prelude::*;
//! use lyndon_rs::prelude::*;
//! use num_rational::Ratio;
//!
//! let basis = LyndonBasis::<ENotation>::new(2, Sort::Lexicographical);
//! let generator = BchSeriesGenerator::<ENotation>::new(basis, 5);
//! let bch_series: LieSeries<ENotation, Ratio<i64>> = generator.generate_lie_series();
//! ```
//!
//! ## Crate Features
//!
//! - `progress`: displays progress bars (via `indicatif`) for long-running
//!   computations.
//! - `tracing`: instruments the potentially expensive operations (basis
//!   factorization, Goldberg and BCH coefficient generation, Lie series
//!   construction, and commutator evaluation) with `tracing` spans, so
//!   timings flow into any subscriber configured by the consuming binary.
//!
//! To see the instrumentation in test output, run e.g.
//! `cargo test -p lie-rs --features tracing -- --nocapture`
//! (respects `RUST_LOG`, defaults to `lie_rs=debug`).
//!
//! # Environment diagnostics
//!
//! A few opt-in environment switches exist for tuning and ops work; all
//! are off by default and do not change numeric results:
//!
//! - `FOLD_ENGINE=0`: routes the class-ordered fold through the per-level
//!   dispatch instead of the default one-dispatch pack-slot walk (A/B
//!   comparison; results are bit-identical either way).
//! - `SLOT_POLICY_DEBUG=1`: one stderr line per walk dispatch with the
//!   slot policy's inputs and decision.
//! - `SIG_NO_COHORT=1` (consumed by `signature-rs`): disables the
//!   4-lane SIMD-across-folds cohort engine, forcing the scalar batch
//!   engine (bit-identical results; the kill switch behind
//!   `cohort_supported`-gated dispatch).

mod bch;
pub mod bch_series_generator;
mod constants;
mod feasible_decompositions;
pub mod lie_series;
pub mod rooted_tree;

use std::ops::{Div, MulAssign};

// Re-export main types at crate root
pub use bch_series_generator::{BchSeriesGenerator, Matrix2x2, MatrixTree};
pub use feasible_decompositions::FeasibleDecompositions;
pub use lie_series::{
    BlockShape, ClassBatchStage, ClassOrder, ClassOrderedCommutation, GatingCache, KernelJob,
    LieSeries, raw_add_assign, raw_add_assign_ptr, raw_mul, commutator_coefficients_batch,
    commutator_coefficients_batch_with_cache, commutator_coefficients_class_batch_with_cache,
    commutator_coefficients_class_fold_with_cache, plan_class_sweep_stages, run_class_batch,
    run_class_batch_with_work, planned_sweep_entries, work_adaptive_slots, IDENTITY_SHIFTS,
};
pub use lie_series::cohort::{add_mul4, supported as cohort_supported, LANES as COHORT_LANES};
pub use rooted_tree::{EdgePartitions, GraphPartitionTable, RootedTree};

use num_traits::{FromPrimitive, One, Zero};

/// Trait for types that can generate Lie series.
///
/// A Lie series is a formal power series in non-commuting variables that
/// respects the Lie algebra structure. This trait provides a uniform
/// interface for generating such series from various sources.
pub trait LieSeriesGenerator<T, U> {
    /// Generates a Lie series from this generator.
    fn generate_lie_series(&self) -> LieSeries<T, U>;
}

/// Computes the binomial coefficient "n choose k" for generic numeric types.
///
/// This function efficiently calculates C(n,k) = n!/(k!(n-k)!) using
/// the multiplicative formula to avoid computing large factorials directly.
pub(crate) fn binomial<U: One + Zero + FromPrimitive + MulAssign + Div<Output = U>>(
    n: usize,
    k: usize,
) -> U {
    if k == 0 {
        return U::one();
    }

    if k > n {
        return U::zero();
    }

    let k = k.min(n - k);
    let mut num = U::one();
    let mut den = U::one();

    for i in 0..k {
        num *= U::from_usize(n - i).expect("Failed to convert from usize");
        den *= U::from_usize(i + 1).expect("Failed to convert from usize");
    }

    num / den
}

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::bch_series_generator::BchSeriesGenerator;
    pub use crate::{LieSeries, LieSeriesGenerator};
}

#[cfg(test)]
mod test {
    use num_rational::Ratio;
    use ordered_float::NotNan;
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(0, 0, 1)]
    #[case(5, 0, 1)]
    #[case(100, 0, 1)]
    #[case(1, 1, 1)]
    #[case(2, 2, 1)]
    #[case(3, 3, 1)]
    #[case(4, 4, 1)]
    #[case(5, 5, 1)]
    #[case(5, 6, 0)]
    #[case(4, 1, 4)]
    #[case(4, 2, 6)]
    #[case(4, 3, 4)]
    #[case(4, 4, 1)]
    #[case(30, 15, 155_117_520)]
    fn test_binomial_integer(#[case] n: usize, #[case] k: usize, #[case] expected_result: i128) {
        assert_eq!(binomial::<i128>(n, k), expected_result);
    }

    #[rstest]
    #[case(0, 0, NotNan::from(1))]
    #[case(5, 0, NotNan::from(1))]
    #[case(100, 0, NotNan::from(1))]
    #[case(1, 1, NotNan::from(1))]
    #[case(2, 2, NotNan::from(1))]
    #[case(3, 3, NotNan::from(1))]
    #[case(4, 4, NotNan::from(1))]
    #[case(5, 5, NotNan::from(1))]
    #[case(4, 1, NotNan::from(4))]
    #[case(4, 2, NotNan::from(6))]
    #[case(4, 3, NotNan::from(4))]
    #[case(4, 4, NotNan::from(1))]
    #[case(30, 15, NotNan::from(155_117_520))]
    fn test_binomial_float(
        #[case] n: usize,
        #[case] k: usize,
        #[case] expected_result: NotNan<f64>,
    ) {
        assert_eq!(binomial::<NotNan<f64>>(n, k), expected_result);
    }

    #[rstest]
    #[case(0, 0, Ratio::from_integer(1))]
    #[case(5, 0, Ratio::from_integer(1))]
    #[case(100, 0, Ratio::from_integer(1))]
    #[case(1, 1, Ratio::from_integer(1))]
    #[case(2, 2, Ratio::from_integer(1))]
    #[case(3, 3, Ratio::from_integer(1))]
    #[case(4, 4, Ratio::from_integer(1))]
    #[case(5, 5, Ratio::from_integer(1))]
    #[case(5, 6, Ratio::from_integer(0))]
    #[case(4, 1, Ratio::from_integer(4))]
    #[case(4, 2, Ratio::from_integer(6))]
    #[case(4, 3, Ratio::from_integer(4))]
    #[case(4, 4, Ratio::from_integer(1))]
    #[case(30, 15, Ratio::from_integer(155_117_520))]
    fn test_binomial_ratio(
        #[case] n: usize,
        #[case] k: usize,
        #[case] expected_result: Ratio<i128>,
    ) {
        assert_eq!(binomial::<Ratio<i128>>(n, k), expected_result);
    }
}
