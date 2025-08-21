pub mod commutator;
use crate::commutator::CommutatorTermPy;
use pyo3::prelude::*;

#[pymodule]
fn commutator_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<CommutatorTermPy>()?;
    Ok(())
}
