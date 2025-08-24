pub mod bch_series_generator;
pub mod lie_series;
use crate::bch_series_generator::BchSeriesGeneratorPy;
use pyo3::prelude::*;

#[pymodule]
fn lie_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<BchSeriesGeneratorPy>()?;
    Ok(())
}
