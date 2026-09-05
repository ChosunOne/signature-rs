//! # signature-rs
//!
//! A Rust library for computing log signatures from path data.
//!
//! ## Quick Start
//!
//! ```rust
//! use signature_rs::prelude::*;
//! use ndarray::array;
//! use ordered_float::NotNan;
//!
//! let path = array![[0.0, 0.0], [1.0, 0.5], [2.0, 1.0]];
//! let path = path.mapv(|v| NotNan::new(v).expect("value to be a number"));
//! let builder = LogSignatureBuilder::<ENotation>::new()
//!     .with_num_dimensions(2)
//!     .with_max_degree(3);
//! let log_sig = builder.build_from_path(&path.view());
//! ```
//!
//! # Plans and caching
//!
//! The Lyndon basis, the structure-constant table, the BCH series and the
//! commutator DAG are pure functions of the builder configuration, so they
//! are built once per process and shared across every subsequent
//! [`LogSignatureBuilder::build`] / [`LogSignatureBuilder::build_from_path`]
//! call: the first build for a configuration pays the expensive table
//! construction, later builds are cheap clones of the cached plan (this is
//! why the coefficient and generator types must be `'static` — the caches
//! key on [`core::any::TypeId`]).
//!
//! # Environment diagnostics
//!
//! `SIG_NO_COHORT=1` disables the 4-lane SIMD-across-folds cohort engine
//! for debugging/A-B (scalar batch engine instead; bit-identical results,
//! see `lie-rs`'s crate docs for the full list of diagnostic switches).

pub(crate) mod commutator_dag;
pub mod log_sig;

// Re-export main types for convenience
pub use log_sig::{LogSignature, LogSignatureBuilder};

// Re-export dependencies that are part of the public API
pub use commutator_rs::{Commutator, CommutatorTerm};
pub use lie_rs::{LieSeries, LieSeriesGenerator};
pub use lyndon_rs::generators::{ENotation, Generator};
pub use lyndon_rs::lyndon::{LyndonBasis, LyndonWord, Sort};

/// Prelude module for imports
pub mod prelude {
    pub use crate::log_sig::{LogSignature, LogSignatureBuilder};
    pub use commutator_rs::Commutator;
    pub use lie_rs::{LieSeries, LieSeriesGenerator};
    pub use lyndon_rs::generators::{ENotation, Generator};
    pub use lyndon_rs::lyndon::{LyndonBasis, LyndonWord, Sort};
}
