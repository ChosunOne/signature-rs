//! # Lyndon-RS
//!
//! A Rust library for working with Lyndon words and the Lyndon basis of the free Lie algebra.
//!
//! ## Overview
//!
//! Lyndon words are a fundamental concept in combinatorics and algebra, providing a canonical
//! basis for the free Lie algebra. This library offers efficient algorithms for:
//!
//! - Generating Lyndon words up to a given length
//! - Constructing the Lyndon basis for free Lie algebras
//! - Working with different generator types (numeric, symbolic, custom)
//!
//! ## Quick Start
//!
//! ```rust,ignore
//! use lyndon_rs::lyndon::{LyndonBasis, Sort};
//! use lyndon_rs::generators::ENotation;
//!
//! // Create a Lyndon basis with 3 generators up to degree 4
//! let basis = LyndonBasis::<ENotation> {
//!     alphabet_size: 3,
//!     max_degree: 4,
//!     sort: Sort::Shortlex,
//! };
//!
//! // Generate all Lyndon words
//! let words = basis.generate_lyndon_words();
//! ```
//!
//! ## Main Components
//!
//! - [`LyndonWord`](lyndon::LyndonWord): Represents a single Lyndon word
//! - [`LyndonBasis`](lyndon::LyndonBasis): Configuration and generation of Lyndon word bases
//! - [`Generator`](generators::Generator): Trait for types that can serve as generators
//! - [`ENotation`](generators::ENotation): Standard e-notation for generators (e1, e2, ...)

pub mod generators;
pub mod lyndon;

// Re-export main types at crate root
pub use generators::{Generator, ENotation};
pub use lyndon::{LyndonBasis, LyndonWord, LyndonWordError, Sort, moebius_mu};

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::generators::{Generator, ENotation};
    pub use crate::lyndon::{LyndonBasis, LyndonWord, Sort};
}
