//! Python bindings for path signature computation.
//!
//! # Classes
//!
//! - **LogSignature**: Path signature computation in log coordinates  
//! - **LogSignatureBuilder**: Factory for creating signatures with configurable parameters

pub mod log_sig;
use pyo3::{prelude::*, sync::GILOnceCell};

use crate::log_sig::{LogSignatureBuilderPy, LogSignaturePy};

pub(crate) static LIE_SERIES_CLASS: GILOnceCell<Py<PyAny>> = GILOnceCell::new();
pub(crate) static LYNDON_WORD_CLASS: GILOnceCell<Py<PyAny>> = GILOnceCell::new();

/// Python module for high-performance path signature computation.
///
/// This module exports:
/// - `LogSignature`: Computed path signatures with efficient operations
/// - `LogSignatureBuilder`: Factory for signature computation with configurable parameters
#[pymodule]
fn signature_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<LogSignatureBuilderPy>()?;
    m.add_class::<LogSignaturePy>()?;
    Ok(())
}
