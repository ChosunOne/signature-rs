pub mod log_sig;
use pyo3::{prelude::*, sync::GILOnceCell};

use crate::log_sig::{LogSignatureBuilderPy, LogSignaturePy};

pub(crate) static LIE_SERIES_CLASS: GILOnceCell<Py<PyAny>> = GILOnceCell::new();
pub(crate) static LYNDON_WORD_CLASS: GILOnceCell<Py<PyAny>> = GILOnceCell::new();

#[pymodule]
fn signature_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<LogSignatureBuilderPy>()?;
    m.add_class::<LogSignaturePy>()?;
    Ok(())
}
