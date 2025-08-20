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
