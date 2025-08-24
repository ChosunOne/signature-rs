use lie_rs::{BchSeriesGenerator, LieSeriesGenerator};
use lyndon_rs::{LyndonBasis, Sort};
use pyo3::{exceptions::PyValueError, prelude::*};

use crate::lie_series::LieSeriesPy;

#[pyclass(name = "BCHSeriesGenerator")]
pub struct BchSeriesGeneratorPy {
    pub inner: BchSeriesGenerator<u8>,
}

#[pymethods]
impl BchSeriesGeneratorPy {
    #[new]
    pub fn new(basis: Bound<'_, PyAny>, max_degree: usize) -> PyResult<Self> {
        let alphabet_size = basis.getattr("alphabet_size")?.extract::<usize>()?;

        let sort = match basis.getattr("sort")?.getattr("value")?.extract::<u8>()? {
            0 => Sort::Lexicographical,
            1 => Sort::Topological,
            _ => return Err(PyValueError::new_err("Invalid Sort")),
        };
        let basis = LyndonBasis::<u8>::new(alphabet_size, sort);
        Ok(Self {
            inner: BchSeriesGenerator::new(basis, max_degree),
        })
    }

    #[must_use]
    pub fn generate_goldberg_coefficient_numerators(&self) -> Vec<i128> {
        self.inner.generate_goldberg_coefficient_numerators()
    }

    #[must_use]
    pub fn generate_lie_series(&self) -> LieSeriesPy {
        LieSeriesPy {
            inner: self.inner.generate_lie_series(),
        }
    }
}
