//! Python bindings for commutator operations and formal indeterminates.
//!
//! # Classes
//!
//! - **CommutatorTerm**: Computation of nested commutators with automatic simplification
//! - **FormalIndeterminate**: Symbolic manipulation of non-commutative expressions

pub mod commutator;
pub mod formal_indeterminate;
use crate::{commutator::CommutatorTermPy, formal_indeterminate::FormalIndeterminatePy};
use pyo3::prelude::*;

/// Python module providing fast commutator operations and formal indeterminates.
///
/// This module exports:
/// - `CommutatorTerm`: Commutator expressions with nested bracket structure
/// - `FormalIndeterminate`: Symbolic representation for non-commutative algebra
#[pymodule]
fn commutator_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<CommutatorTermPy>()?;
    m.add_class::<FormalIndeterminatePy>()?;
    Ok(())
}
