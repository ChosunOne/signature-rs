pub mod commutator;
pub mod formal_indeterminate;
use crate::{commutator::CommutatorTermPy, formal_indeterminate::FormalIndeterminatePy};
use pyo3::prelude::*;

#[pymodule]
fn commutator_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<CommutatorTermPy>()?;
    m.add_class::<FormalIndeterminatePy>()?;
    Ok(())
}
