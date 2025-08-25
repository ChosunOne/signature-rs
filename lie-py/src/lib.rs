pub mod bch_series_generator;
pub mod lie_series;
use crate::bch_series_generator::BchSeriesGeneratorPy;
use crate::lie_series::LieSeriesPy;
use pyo3::{prelude::*, sync::GILOnceCell};

pub(crate) static LYNDON_WORD_CLASS: GILOnceCell<Py<PyAny>> = GILOnceCell::new();
pub(crate) static COMMUTATOR_TERM_CLASS: GILOnceCell<Py<PyAny>> = GILOnceCell::new();

#[pymodule]
fn lie_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<BchSeriesGeneratorPy>()?;
    m.add_class::<LieSeriesPy>()?;
    Ok(())
}
