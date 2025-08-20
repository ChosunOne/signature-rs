//! # commutator-rs
//!
//! A Rust library for commutator operations.
//!
//! ## Quick Start
//!
//! ```rust
//! use commutator_rs::{prelude::*, comm};
//!
//! let x = CommutatorTerm::Atom { coefficient: 1, atom: 'x' };
//! let y = CommutatorTerm::Atom { coefficient: 1, atom: 'y' };
//!
//! // Using the trait method
//! let result1 = x.commutator(&y);
//!
//! // Using the comm! macro
//! let result2 = comm![x, y];
//! assert_eq!(result1, result2);
//! ```

pub mod commutator;
pub mod formal_indeterminate;

// Re-export main types at crate root
pub use commutator::{Commutator, CommutatorTerm};
pub use formal_indeterminate::FormalIndeterminate;

/// Prelude module for imports
pub mod prelude {
    pub use crate::commutator::{Commutator, CommutatorTerm};
    pub use crate::formal_indeterminate::FormalIndeterminate;
}
