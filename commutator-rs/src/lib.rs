//! # Commutator-RS
//!
//! A Rust library for computing commutators in non-commutative algebras.
//!
//! ## Overview
//!
//! This library provides tools for working with commutators, which are fundamental
//! operations in non-commutative algebra defined as [a, b] = ab - ba. Features include:
//!
//! - Generic commutator trait for custom types
//! - Formal indeterminates for symbolic computation
//! - Efficient representation of commutator expressions
//!
//! ## Quick Start
//!
//! ```rust,ignore
//! use commutator_rs::Commutator;
//! use commutator_rs::formal_indeterminate::FormalIndeterminate;
//!
//! // Create formal indeterminates
//! let x = FormalIndeterminate::new("x", 1.0);
//! let y = FormalIndeterminate::new("y", 1.0);
//!
//! // Compute commutator [x, y] using the trait
//! let result = x.commutator(&y);
//! ```
//!
//! ## Main Components
//!
//! - [`Commutator`]: Trait for types supporting commutator operations
//! - [`CommutatorTerm`]: Represents terms in commutator expressions
//! - [`FormalIndeterminate`](formal_indeterminate::FormalIndeterminate): Symbolic variables for formal computations
//!
//! ## Mathematical Background
//!
//! The commutator [a, b] = ab - ba measures how much two elements fail to commute.
//! It is zero if and only if the elements commute (ab = ba). Commutators satisfy
//! the Jacobi identity and form the basis of Lie algebra structure.

pub mod commutator;
pub mod formal_indeterminate;

// Re-export main types at crate root
pub use commutator::{Commutator, CommutatorTerm};
pub use formal_indeterminate::FormalIndeterminate;

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::commutator::{Commutator, CommutatorTerm};
    pub use crate::formal_indeterminate::FormalIndeterminate;
}
