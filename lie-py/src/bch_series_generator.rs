use lie_rs::{BchSeriesGenerator, LieSeriesGenerator};
use lyndon_rs::{LyndonBasis, Sort};
use pyo3::{exceptions::PyValueError, prelude::*, types::PyList};

use crate::{LYNDON_WORD_CLASS, lie_series::LieSeriesPy};

#[pyclass(name = "BCHSeriesGenerator")]
#[derive(Clone, Debug)]
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

    #[getter]
    #[must_use]
    pub fn get_alphabet_size(&self) -> usize {
        self.inner.alphabet_size
    }

    #[getter]
    #[must_use]
    pub fn get_max_degree(&self) -> usize {
        self.inner.max_degree
    }

    #[getter]
    #[must_use]
    pub fn get_left_factor(&self) -> &[usize] {
        &self.inner.left_factor
    }

    #[getter]
    #[must_use]
    pub fn get_right_factor(&self) -> &[usize] {
        &self.inner.right_factor
    }

    #[getter]
    #[must_use]
    pub fn get_word_lengths(&self) -> &[usize] {
        &self.inner.word_lengths
    }

    #[getter]
    #[must_use]
    pub fn get_index_of_degree(&self) -> &[usize] {
        &self.inner.index_of_degree
    }

    #[getter]
    #[must_use]
    pub fn multi_degree(&self) -> &[usize] {
        &self.inner.multi_degree
    }

    #[getter]
    pub fn get_basis(&self, py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
        let lyndon_word_class = LYNDON_WORD_CLASS.get_or_try_init(py, || {
            let module = py.import("lyndon_py")?;
            let class = module.getattr("LyndonWord")?;
            PyResult::Ok(class.unbind())
        })?;
        self.inner
            .basis
            .iter()
            .map(|w| PyList::new(py, &w.letters))
            .collect::<PyResult<Vec<_>>>()?
            .into_iter()
            .map(|l| lyndon_word_class.call1(py, (l,)))
            .collect::<PyResult<Vec<_>>>()
    }
}
