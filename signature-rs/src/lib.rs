//! # Signature-RS
//!
//! A high-performance Rust library for computing log signatures and related algebraic operations
//! on paths and time series data.
//!
//! ## Overview
//!
//! The log signature is a powerful mathematical tool from rough path theory that provides
//! a compact representation of sequential data by capturing its geometric properties through
//! iterated integrals. This library provides efficient implementations of:
//!
//! - Log signature computation from paths
//! - Lyndon word basis for the free Lie algebra
//! - Baker-Campbell-Hausdorff formula computations
//! - Commutator operations in non-commutative algebras
//!
//! ## Quick Start
//!
//! ```rust,ignore
//! use signature_rs::prelude::*;
//! use ndarray::array;
//!
//! // Create a path in 2D space
//! let path = array![[0.0, 0.0], [1.0, 0.5], [2.0, 1.0]];
//!
//! // Build and compute log signature up to degree 3
//! let builder = LogSignatureBuilder::new()
//!     .with_num_dimensions(2)
//!     .with_max_degree(3);
//!
//! let log_sig = builder.build_from_slice(&path.view());
//! ```
//!
//! ## Features
//!
//! - **Performance**: Optimized for speed with parallel computation support
//! - **Flexibility**: Generic over numeric types supporting various precision requirements
//! - **Mathematical Rigor**: Based on established algebraic structures from rough path theory
//! - **Comprehensive**: Includes full support for Lyndon basis, BCH formula, and Lie series
//!
//! ## Main Components
//!
//! - [`LogSignature`]: The main type for log signature computation
//! - [`LogSignatureBuilder`]: Builder pattern for configuring log signatures
//! - [`LyndonBasis`]: Basis for the free Lie algebra using Lyndon words
//! - [`LieSeries`]: Formal Lie series with algebraic operations
//!
//! ## Mathematical Background
//!
//! The log signature transforms a path into an element of the free Lie algebra,
//! providing a graded summary of the path that is invariant under time reparametrization
//! and captures essential geometric information about the path's evolution.

pub mod log_sig;

// Re-export main types for convenience
pub use log_sig::{LogSignature, LogSignatureBuilder};

// Re-export dependencies that are part of the public API
pub use lie_rs::{LieSeries, LieSeriesGenerator};
pub use lyndon_rs::lyndon::{LyndonBasis, LyndonWord, Sort};
pub use lyndon_rs::generators::{Generator, ENotation};
pub use commutator_rs::{Commutator, CommutatorTerm};

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::log_sig::{LogSignature, LogSignatureBuilder};
    pub use lie_rs::{LieSeries, LieSeriesGenerator};
    pub use lyndon_rs::lyndon::{LyndonBasis, LyndonWord, Sort};
    pub use lyndon_rs::generators::Generator;
    pub use commutator_rs::Commutator;
}
