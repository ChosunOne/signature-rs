use numpy::PyReadonlyArray2;
use ordered_float::NotNan;
use pyo3::{exceptions::PyValueError, prelude::*, types::PyList};
use signature_rs::{LogSignature, LogSignatureBuilder};

use crate::{LIE_SERIES_CLASS, LYNDON_WORD_CLASS};

#[pyclass(name = "LogSignatureBuilder")]
#[derive(Default, Copy, Clone, Debug)]
pub struct LogSignatureBuilderPy {
    inner: LogSignatureBuilder<u8>,
}

#[pymethods]
impl LogSignatureBuilderPy {
    #[new]
    #[must_use]
    pub fn new(max_degree: Option<usize>, num_dimensions: Option<usize>) -> Self {
        let inner = LogSignatureBuilder::new()
            .with_max_degree(max_degree.unwrap_or_default())
            .with_num_dimensions(num_dimensions.unwrap_or_default());
        Self { inner }
    }

    #[getter]
    #[must_use]
    pub fn get_max_degree(&self) -> usize {
        self.inner.max_degree
    }

    #[setter]
    pub fn set_max_degree(&mut self, max_degree: usize) {
        self.inner.max_degree = max_degree;
    }

    #[getter]
    #[must_use]
    pub fn get_num_dimensions(&self) -> usize {
        self.inner.num_dimensions()
    }

    #[setter]
    pub fn set_num_dimensions(&mut self, num_dimensions: usize) {
        self.inner.lyndon_basis.alphabet_size = num_dimensions;
    }

    #[must_use]
    pub fn build(&self) -> LogSignaturePy {
        LogSignaturePy {
            inner: self.inner.build(),
        }
    }

    pub fn build_from_path(&self, path: PyReadonlyArray2<'_, f32>) -> PyResult<LogSignaturePy> {
        if path.as_array().is_any_nan() {
            return Err(PyValueError::new_err("Received NaN float value."));
        }
        let not_nan_path = path.as_array().mapv(|x| NotNan::try_from(x).unwrap());
        Ok(LogSignaturePy {
            inner: self.inner.build_from_path(&not_nan_path.view()),
        })
    }
}

#[pyclass(name = "LogSignature")]
#[derive(Debug, Clone)]
pub struct LogSignaturePy {
    pub inner: LogSignature<u8, NotNan<f32>>,
}

#[pymethods]
impl LogSignaturePy {
    #[must_use]
    pub fn __getitem__(&self, idx: usize) -> NotNan<f32> {
        self.inner[idx]
    }

    pub fn __setitem__(&mut self, idx: usize, coefficient: NotNan<f32>) {
        self.inner[idx] = coefficient;
    }

    #[getter]
    pub fn get_series(&self, py: Python) -> PyResult<Py<PyAny>> {
        let lie_series_class = LIE_SERIES_CLASS.get_or_try_init(py, || {
            let module = py.import("lie_py")?;
            let class = module.getattr("LieSeries")?;
            PyResult::Ok(class.unbind())
        })?;
        let lyndon_word_class = LYNDON_WORD_CLASS.get_or_try_init(py, || {
            let module = py.import("lyndon_py")?;
            let class = module.getattr("LyndonWord")?;
            PyResult::Ok(class.unbind())
        })?;
        let basis = self
            .inner
            .series
            .basis
            .iter()
            .map(|w| PyList::new(py, &w.letters))
            .collect::<PyResult<Vec<_>>>()?
            .into_iter()
            .map(|l| lyndon_word_class.call1(py, (l,)))
            .collect::<PyResult<Vec<_>>>()?;
        lie_series_class.call1(py, (basis, self.inner.series.coefficients.clone()))
    }

    #[getter]
    pub fn get_bch_series(&self, py: Python) -> PyResult<Py<PyAny>> {
        let lie_series_class = LIE_SERIES_CLASS.get_or_try_init(py, || {
            let module = py.import("lie_py")?;
            let class = module.getattr("LieSeries")?;
            PyResult::Ok(class.unbind())
        })?;
        let lyndon_word_class = LYNDON_WORD_CLASS.get_or_try_init(py, || {
            let module = py.import("lyndon_py")?;
            let class = module.getattr("LyndonWord")?;
            PyResult::Ok(class.unbind())
        })?;
        let basis = self
            .inner
            .bch_series
            .basis
            .iter()
            .map(|w| PyList::new(py, &w.letters))
            .collect::<PyResult<Vec<_>>>()?
            .into_iter()
            .map(|l| lyndon_word_class.call1(py, (l,)))
            .collect::<PyResult<Vec<_>>>()?;
        lie_series_class.call1(py, (basis, self.inner.bch_series.coefficients.clone()))
    }

    #[must_use]
    pub fn concatenate(&self, rhs: &Self) -> Self {
        Self {
            inner: self.inner.concatenate(&rhs.inner),
        }
    }

    pub fn concatenate_assign(&mut self, rhs: &Self) {
        self.inner.concatenate_assign(&rhs.inner);
    }
}
