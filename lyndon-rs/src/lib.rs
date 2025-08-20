//! # lyndon-rs
//!
//! A Rust library for working with Lyndon words.
//!
//! ## Quick Start
//!
//! ```rust
//! use lyndon_rs::prelude::*;
//!
//! let basis = LyndonBasis::<ENotation>::new(3, Sort::Lexicographical);
//! let words = basis.generate_basis(4);
//! ```

pub mod generators;
pub mod lyndon;

// Re-export main types at crate root
pub use generators::{ENotation, Generator};
pub use lyndon::{LyndonBasis, LyndonWord, LyndonWordError, Sort, moebius_mu};

/// Prelude module for imports
pub mod prelude {
    pub use crate::generators::{ENotation, Generator};
    pub use crate::lyndon::{LyndonBasis, LyndonWord, Sort};
}
